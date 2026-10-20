#'@name model_ssGBLUP_multisite.R
#'@author Joane Elleouet

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(ASRgenomics)
library(AGHmatrix)
library(asreml)
library(ASRtriala)
library(corrplot)

asreml.options(threads=20)

#load data ---------------------------------------------------------------------
data_dir <- "data/"
out_dir <- "out/"
phenotype <- "vol.m3" # vol.m3 | dbh.m | convexhull3D_area | zmax.m
trials <- c("trial_1","trial_2")

dat0 <- readRDS(
  paste0(data_dir,"phenotypes_anonymised.rds")
  ) %>%
  filter(trial %in% trials) %>%
  mutate(trial = as.factor(trial))

dat <- dat0%>%
  filter(!grepl("ML|DL|FK|BW|BL|BT|DF|LL|LR|TO", tree_form_code)) %>%
  mutate(convexhull3D_area_rescale = convexhull3D_area/1000)


#Make inverse numerator matrix -------------------------------------------------
ped <- readRDS(
  paste0(data_dir,"pedigree_anonymised.rds")
  ) %>%
  mutate(across(everything(), as.factor))

#A
A <- AGHmatrix::Amatrix(data = ped, ploidy = 2, dominance = FALSE)

#G
Gsel <- readRDS(
  paste0(data_dir,"Gmatrix_anonymised.rds")
)

G_inv <- G.inverse(G = Gsel, sparseform = FALSE)$Ginv

#prepare H matrix
H <- H.inverse(A = A, G = G_inv, lambda = 0.5, sparseform = TRUE)
rm(G_align, Gsel, A)

# ------------------------------------------------------------------------------

if(phenotype == "convexhull3D_area"){
  phenotype2 <- "convexhull3D_area_rescale"
} else{
    phenotype2 <- phenotype
}

# with additive genetic variance of genotypes only
################################################
met1 <- asreml(fixed = as.formula(paste(phenotype2, "~ trial")),
               random = ~ at(trial):grid_id + corgh(trial):vm(GenotypeName,H),
               residual = ~ dsum(~units | trial),
               data = dat,
               workspace = "4gb")

saveRDS(met1,
        paste0(out_dir,
               phenotype,
               "_H0.5_met1.rds"
               )
)

update_thres <- 3
count <- 0
while (!met1$converge & count<update_thres) {
  met1 <- update(met1)
  count <- count + 1
  saveRDS(met1,
          paste0(out_dir,
                 phenotype,
                 "_H0.5_met1",
                 "_update",
                 count,
                 ".rds"
                 )
  )
}

#extract variance components
met1_vc <- summary(met1)$varcomp %>%
  rownames_to_column(var="Component") %>%
  select(Component, Value=component, Constraint = bound)


# with nonadditive genetic variance of genotypes
################################################

#initialise parameter table
start2 <- asreml(fixed = as.formula(paste(phenotype2, "~ trial")),
                 random = ~ at(trial):grid_id + corgh(trial):vm(GenotypeName,H) +
                   corgh(trial):ide(GenotypeName),
                 residual = ~ dsum(~units | trial),
                 data = dat,
                 start.values = TRUE,
                 workspace = "1gb")

start2.df <- start2$vparameters.table

#initialise table
start2.df <- start2.df %>%
  rows_update(met1_vc, by="Component")

#run model
met2 <- asreml(fixed = as.formula(paste(phenotype2, "~ trial")),
               random = ~ at(trial):grid_id +
                 corgh(trial):vm(GenotypeName, H) +
                 corgh(trial):ide(GenotypeName),
               residual = ~ dsum(~units | trial),
               data = dat,
               G.param = start2.df,
               R.param = start2.df,
               workspace = "4gb")
saveRDS(met2,
        paste0(out_dir,
               phenotype,
               "_H0.5_met2.rds"
        )
)

update_thres <- 3
count <- 0
while (!met2$converge & count<update_thres) {
  met2 <- update(met2)
  count <- count + 1
  saveRDS(met2,
          paste0(out_dir,
                 phenotype,
                 "_H0.5_met1",
                 "_update",
                 count,
                 ".rds"
          )
  )
}

