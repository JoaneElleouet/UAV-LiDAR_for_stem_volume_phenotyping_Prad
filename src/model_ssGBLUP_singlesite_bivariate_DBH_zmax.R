#'@name model_ssGBLUP_singlesite_bivariate_DBH_zmax.R
#'@author Joane Elleouet

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(ASRgenomics)
library(AGHmatrix)
library(asreml)
library(corrplot)
library(reshape2)
library(purrr)

asreml.options(threads=20)

source("src/ped.extract.R")

data_dir <- "data/"
out_dir <- "out/"
trials <- c("trial_1","trial_2")
pheno_vars <-  c(dbh = "dbh.m", zmax = "zmax.m")

ped <- readRDS(
  paste0(data_dir,"pedigree_anonymised.rds")
) %>%
  mutate(across(everything(), as.factor))

dat0 <- readRDS(
  paste0(data_dir,"phenotypes_anonymised.rds")
) %>%
  filter(trial %in% trials) %>%
  mutate(trial = as.factor(trial))

dat <- dat0%>%
  filter(!grepl("ML|DL|FK|BW|BL|BT|DF|LL|LR|TO", tree_form_code))

dat_list <- split(dat, dat$trial)

#Load G matrix -----------------------------------------------------------------
G_align <- readRDS(
  paste0(data_dir,"Gmatrix_anonymised.rds")
)

# prepare variance components --------------------------------------------------
# uses varcomp output from model_ssGBLUP_singlesite.R
#load most recent model results

open_varcomp_res <- function(x){
  varcom_filename <- paste0(out_dir, x, "_H0.5_asr2_varcomp.csv")
  read.csv(varcom_filename)
}

vc_list <- lapply(pheno_vars, open_varcomp_res)
names(vc_list) <- pheno_vars

cor_pheno_list <- readRDS("out/phenotypic_correlations_single_trials.rds")
cor_pheno_list_sel <- lapply(cor_pheno_list,
                             function(mat){
                               mat_sub <- mat[c("dbh.cm","zmax.dm"), c("dbh.cm","zmax.dm")]
                               mat_sub[lower.tri(mat_sub, diag = FALSE)]})

# function to run model(s)
run_asr <- function(trial_dat,
                    Gmatrix = G_align,
                    pheno_vars = pheno_vars,
                    vc_list,
                    cor_pheno = cor_pheno_list_sel){

  trial_name <- trial_dat$trial[1]

  trial_dat <- trial_dat %>%
    mutate(across(c(GenotypeName, grid_id), droplevels))

  # genetic relationship matrix --------------------------------------------------
  #prepare A matrix
  ped_trial <- ped.extract(x = trial_dat[["GenotypeName"]], reference = ped)[[1]] %>%
    mutate(across(everything(), as.factor))
  A_trial <- AGHmatrix::Amatrix(data = ped_trial, ploidy = 2, dominance = FALSE)

  #prepare G matrix
  G_trial <- G_align[which(rownames(G_align) %in% ped_trial$Tree),
                which(colnames(G_align) %in% ped_trial$Tree)]
  G_trial_inv <- G.inverse(G = G_trial, sparseform = FALSE)$Ginv
  #prepare H matrix
  H_trial <- H.inverse(A = A_trial, G = G_trial_inv, lambda = 0.5, sparseform = TRUE)
  rm(A_trial, G_trial, G_trial_inv)

  #starting values
  vc_list_trial <- lapply(vc_list,
                          function(vc){vc %>% filter(site_year == trial_name)})

  cor_trial <- cor_pheno[[trial_name]]

  init_ga <- c(cor_trial, sapply(vc_list_trial, function(x) x$component[2]))
  init_ga2 <- init_ga;init_ga2[2] <- init_ga[2]*100;init_ga2[3] <- init_ga[3]*10;
  init_gna <- c(cor_trial, sapply(vc_list_trial, function(x) x$component[3]))
  init_gna2 <- init_gna;init_gna2[2] <- init_gna[2]*100;init_gna2[3] <- init_gna[3]*10;
  init_block <- sapply(vc_list_trial, function(x) x$component[1])
  init_block2 <- init_block;init_block[1] <- init_block[1]*100;init_block2[2] <- init_block2[2]*10;
  init_r <- c(cor_trial, sapply(vc_list_trial, function(x) x$component[4]))
  init_r2 <- init_r;init_r2[2] <- init_r[2]*100;init_r2[3] <- init_r[3]*10;

  #response formula
  lhs <- paste0("cbind(", paste(pheno_vars, collapse = ", "), ")")

  fix_f <- as.formula(paste(lhs, "~ trait"))

  random_f <- as.formula(
    bquote(
      ~ idh(trait, init = .(init_block2)):id(grid_id) +
        corgh(trait, init = .(init_ga2)):vm(GenotypeName, H_trial) +
        corgh(trait, init = .(init_gna2)):id(GenotypeName)
    )
  )

  resid_f <- as.formula(
    bquote(
      ~ id(units):corgh(trait, init = .(as.numeric(init_r)))
    )
  )

  #model
  m1 <- asreml(fixed = fix_f,
               random    = random_f,
               residual = resid_f,
               data      = trial_dat,
               maxit = 100,
               workspace = "4gb")

  saveRDS(m1,
          paste0(out_dir,
                 "bivariate_dbh_zmax_H0.5_",
                 trial_name,
                 ".rds")
          )
  m1
}

safe_run <- possibly(run_asr, otherwise = NULL)

asr_out_list <- map(
  dat_list, ~ safe_run(trial_dat = .x,
                       Gmatrix = G_align,
                       pheno_vars = pheno_vars,
                       vc_list = vc_list,
                       cor_pheno = cor_pheno_list_sel))



saveRDS(asr_out_list,
        paste0(out_dir,
               "bivariate_dbh_zmax_H0.5_list.rds")
)


sapply(asr_out_list, is.null)
asr_out_list2 <- asr_out_list[!vapply(asr_out_list, is.null, logical(1))]
lapply(asr_out_list2, function(x){x$converge})
lapply(asr_out_list2, function(x){summary(x)$varcomp})

#extract genetic correlations
gencor_list <- lapply(asr_out_list2, function(x){summary(x)$varcomp["trait:vm(GenotypeName, H_trial)!trait!zmax.m:!trait!dbh.m.cor",]})
gencor_df <- do.call(rbind,gencor_list) %>%
  rownames_to_column(var="trial")

