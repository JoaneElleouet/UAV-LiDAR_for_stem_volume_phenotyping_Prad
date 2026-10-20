#'@name calculate_ebv_function.R
#'@author Joane Elleouet
#'@description Calculate Estimated Breeding Values (EBVs) and deployment values
#'and their associated environmental values using estimated trial means
#'
#' @param dat_pheno Data frame containing phenotypic observations.
#' @param phenotype Character string specifying the phenotype to analyse.
#' @param trial_series Character. trial series.
#' @param trial_names Optional character vector of trial names to include.
#'
#' @return data frame containing EBVs (and associated statistics, if fitted)
#'   for each genotype.

library(dplyr)
library(tidyr)
library(tibble)
library(asreml)
library(ggplot2)

calculate_ebv <- function(dat_pheno, phenotype, trial_series, trial_names=NULL, out_dir) {

  model <- "met2"

  trials <- if (is.null(trial_names)) {
    if (trial_series == "BC52-55-59") {
      c("BC52_1","BC52_3","BC52_5","BC52_6","BC55_2","BC55_4","BC59_3")
    } else{
      c("BC66_1","BC66_2","BC66_4","BC66_5")}
  } else {trial_names}

  source("src/ped.extract.R")

  # load data --------------------------------------------------------------------
  dat_pheno <- dat_pheno %>%
    filter(trial %in% trials)


  #load most recent model results ------------------------------------------------
  prefix <- paste0(phenotype,"_H0.5_", model)
  asr_out <- readRDS(paste0(out_dir,phenotype,"_H0.5_", model,".rds"))
  if (!asr_out$converge) {
    stop(paste0("the model object in ",latest, " did not converge"))
  }


  #variance components -----------------------------------------------------------
  Vc <- summary(asr_out)$varcomp %>%
    as.data.frame() %>%
    rownames_to_column(var = "variance_component")

  #trial-level additive genetic variance
  Va_trial <- Vc %>%
    filter(grepl("trial:vm\\(GenotypeName, H\\)!trial_", variance_component)) %>%
    mutate(trial = gsub("trial:vm\\(GenotypeName, H\\)!trial_", "", variance_component)) %>%
    select(trial,additive_variance = component)


  #trial means -------------------------------------------------------------------
  trial_effects <- asr_out$coefficients$fixed %>% as.data.frame() %>%
    slice(-1) %>%
    rownames_to_column(var="trial") %>%
    mutate(trial = gsub("trial_","",trial))
  intercept <- asr_out$coefficients$fixed %>% as.data.frame() %>%
    slice(1)

  trial_intercepts <- trial_effects %>%
    mutate(trial_means = effect + intercept[1,1])

  write.csv(trial_intercepts,
            paste0(out_dir, phenotype, "_H0.5_", model, "_trial_intercepts.csv"),
            row.names=FALSE)

  # EBVs -------------------------------------------------------------------------

  #extract random effects
  re <- summary(asr_out, coef=TRUE)$coef.random %>%
    as.data.frame() %>%
    rownames_to_column(var = "effect_name")

  #extract blups with standard deviations
  blup_df <- re%>%
    filter(grepl("vm",re$effect_name)) %>%
    separate(effect_name, into = c("trial", "genotype"), sep = ":") %>%
    mutate(trial = gsub("trial_","",trial),
           genotype = gsub("vm\\(GenotypeName, H\\)_","",genotype)) %>%
    select(-z.ratio) %>%
    rename(add_effect = solution, add_effect_se = std.error) %>%
    mutate(add_effect_se_rel = add_effect_se/add_effect)

  #extract edv
  nonadd_df <- re %>%
    filter(grepl("ide",re$effect_name)) %>%
    separate(effect_name, into = c("trial", "genotype"), sep = ":") %>%
    mutate(trial = gsub("trial_","",trial),
           genotype = gsub("ide\\(GenotypeName\\)_","",genotype)) %>%
    select(-z.ratio) %>%
    rename(nonadd_effect = solution, nonadd_effect_se = std.error)

  #put together relevant effects
  perf_df <- blup_df %>%
    left_join((dat_pheno %>% select(genotype = GenotypeName) %>% distinct())) %>%
    left_join(nonadd_df) %>%
    left_join(trial_intercepts %>% select(-effect))

  #calculate env-BV and env_DV
  perf_df <- perf_df %>%
    mutate(genetic_effect = add_effect + nonadd_effect) %>%
    mutate(rel_add_effect = add_effect/trial_means,
           rel_add_effect_se = add_effect_se/trial_means,
           rel_genetic_effect = genetic_effect/trial_means,
           env_bv = add_effect + trial_means,
           env_dv = genetic_effect + trial_means,
           rel_env_bv = env_bv/trial_means,
           rel_env_dv = env_dv/trial_means)


  # accuracy ---------------------------------------------------------------------
  perf_df <- perf_df %>%
    left_join(Va_trial) %>%
    mutate(PEV = add_effect_se^2,
           bv_accuracy = ifelse(PEV/additive_variance<=1, sqrt(1 - PEV/additive_variance), 0)) %>%
    select(-additive_variance)


  #identify parents
  ped <- ped.extract(x = dat_pheno[["GenotypeName"]], reference = ped)[[1]] %>%
    mutate(across(everything(), as.factor))

  perf_df <- perf_df %>%
    mutate(status = ifelse(genotype %in% c(ped$Mom,ped$Dad), "parent", "clone"))

  #save
  write.csv(perf_df, paste0("out/",phenotype,"_H0.5_", model, "_ebvs.csv"), row.names=FALSE)

  perf_df2 <- perf_df %>%
    mutate(trial_series = trial_series,
           phenotype = phenotype)

  return(perf_df2)
}
