#'@name model_ssGBLUP_singlesite.R
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

asreml.options(threads=20)
source("src/ped.extract.R")


#load data ---------------------------------------------------------------------
data_dir <- "data/"
out_dir <- "out/"
phenotype <- "zmax.m"
trials <- c("trial_1","trial_2")

dat0 <- readRDS(
  paste0(data_dir,"phenotypes_anonymised.rds")
) %>%
  filter(trial %in% trials) %>%
  mutate(trial = as.factor(trial))

dat <- dat0%>%
  filter(!grepl("ML|DL|FK|BW|BL|BT|DF|LL|LR|TO", tree_form_code))


#Make inverse numerator matrix -------------------------------------------------
ped <- readRDS(
  paste0(data_dir,"pedigree_anonymised.rds")
) %>%
  mutate(across(everything(), as.factor))

#A
A <- AGHmatrix::Amatrix(data = ped, ploidy = 2, dominance = FALSE)

#G
G_align <- readRDS(
  paste0(data_dir,"Gmatrix_anonymised.rds")
)


# model  - single sites -------------------------------------------------------

#prepare data
dat_list <- split(dat, dat$trial)

#model
run_asr <- function(x){
  x <- x %>% mutate(across(c(grid_id,GenotypeName), droplevels))
  #prepare A matrix
  pedx <- ped.extract(x = x[["GenotypeName"]], reference = ped)[[1]] %>%
    mutate(across(everything(), as.factor))
  Ax <- AGHmatrix::Amatrix(data = pedx, ploidy = 2, dominance = FALSE)

  #prepare G matrix
  Gx <- G_align[which(rownames(G_align) %in% pedx$Tree),
                which(colnames(G_align) %in% pedx$Tree)]
  Gx_inv <- G.inverse(G = Gx, sparseform = FALSE)$Ginv
  #prepare H matrix
  Hx <- H.inverse(A = Ax, G = Gx_inv, lambda = 0.5, sparseform = TRUE)
  m1 <- asreml(fixed = as.formula(paste(phenotype, "~ 1")),
               random = ~ grid_id + vm(GenotypeName,Hx) + ide(GenotypeName),
               data = x)
  m1 <- update(m1)
  m1
  }

asr2_list <- lapply(dat_list, run_asr)

sapply(asr2_list, function(x){x$converge})
lapply(asr2_list, function(x){summary(x)$varcomp})

#get log likelihoods
loglik_list <- lapply(asr2_list, function(x){summary(x)$loglik})
loglik_df <- data.frame(trial = names(loglik_list),
                        loglik = do.call(c,loglik_list))
write.csv(loglik_df,
          paste0(out_dir,
                 phenotype,
                 "_H0.5_asr2_loglik.csv"
          )
)

#get heritabilities
#narrow-sense

asr2_h2_list <- lapply(asr2_list, function(x) {
  vpredict(x, h2 ~ V2/(V2 + V3 + V4))})
asr2_h2 <- do.call(bind_rows, asr2_h2_list) %>%
  mutate(trial = names(asr2_list))
#broad-sense
asr2_H2_list <- lapply(asr2_list, function(x) {
  vpredict(x, H2 ~ (V2 + V3)/(V2 + V3 + V4))})
asr2_H2 <- do.call(bind_rows, asr2_H2_list) %>%
  mutate(trial = names(asr2_list))


#for output
asr2_h2_out <- asr2_h2 %>%
  mutate(trait = phenotype) %>%
  select(trial, Trait = trait, h2_Estimate = Estimate, h2_SE = SE) %>%
  left_join((asr2_H2 %>% select(trial, H2 = Estimate, H2_SE = SE))) %>%
  mutate(across(c(h2_Estimate,h2_SE,H2,H2_SE), function(x){round(x,2)}))
write.csv(asr2_h2_out,
          paste0(out_dir,
                 phenotype,
                 "H0.5_asr2_h2.csv"
          ),
          row.names=FALSE
)

# checks
test_normality <- lapply(asr2_list, function(x){shapiro.test(x$residuals)})
asr2_plot <- lapply(asr2_list, plot)
lapply(asr2_plot, function(x){x$qq})

#check outliers
lapply(asr2_plot, function(x){x$fitted})

#output varcomp from all single models
asr2_varcomp_list <- lapply(asr2_list, function(x){summary(x)$varcomp %>% rownames_to_column(var = "variance_component")})
asr2_varcomp_list2 <- lapply(names(asr2_list), function(x) {
  df <- asr2_varcomp_list[[x]]
  df$site_year <- x  # Add new column with the list element's name
  return(df)
})

varcomp_out <- bind_rows(asr2_varcomp_list2) %>%
  select(7,1,2,3,4,5,6) %>%
  mutate(across(c(component, std.error, z.ratio), function(x){format(x, scientific = TRUE, digits = 2)}))

write.csv(varcomp_out,
          paste0(out_dir,
                 phenotype,
                 "_H0.5_asr2_varcomp.csv"
          ),
          row.names = FALSE
)
