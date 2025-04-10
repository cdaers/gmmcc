library(mclust)
source("R/utils.R")
source("R/model.R")

# ! ------------------------- run demo ------------------------- ! #

# Real data analysis part
# demo 1:  Binary classification. LM-GMM achieves better result than ML.
dat_nm <- "brcancer_prog"
dat <- read.csv(sprintf("inst/extdata/real/%s.csv", dat_nm), header = T)
model_result <- model_run(dat, dat_nm = dat_nm,
                          features_ignored_gmm = "Tumor.size",
                          gmm_args = list(num_chains = 3),
                          model_saved = TRUE)
# model stability test
test_result <- stab_test_run(dat, dat_nm, model_result)
test_result$template_chosen # template chosen
quantile(test_result$result$Accuracy) # Accuracy
quantile(test_result$result$AUC) # AUC

# - demo 2 : Binary classification. Both ML and LM-GMM achieve good results.
dat_nm <- "brcancer_diag"
dat <- read.csv(sprintf("inst/extdata/real/%s.csv", dat_nm), header = T)
model_result <- model_run(dat, dat_nm = dat_nm,
                          lr_params = list(penalty = TRUE, min_features = 5),
                          model_saved = TRUE, opt_metric = "AUC")
# model stability test
test_result <- stab_test_run(dat, dat_nm, model_result)
test_result$template_chosen # template chosen
quantile(test_result$result$Accuracy) # Accuracy
quantile(test_result$result$AUC) # AUC

# - demo 3 : Multiclass classification. Both ML and LM-GMM and ML achieve good results.
dat_nm <- "Leukemia"
dat <- read.csv(sprintf("inst/extdata/real/%s.csv", dat_nm), header = T)
model_result <- model_run(dat, dat_nm = dat_nm,
                          features_ignored_gmm = c("Monocytes_percent", "Monocytes_G_L"),
                          model_saved = TRUE, opt_metric = "AUC")
# model stability test
test_result <- stab_test_run(dat, dat_nm, model_result)
test_result$template_chosen # template chosen
quantile(test_result$result$Accuracy) # Accuracy
quantile(test_result$result$AUC) # AUC

# - demo 4 : Binary classification. Both ML and LM-GMM and ML achieve good results.
dat_nm <- "bank"
dat <- read.csv(sprintf("inst/extdata/real/%s.csv", dat_nm), header = T)
model_result <- model_run(dat, dat_nm = dat_nm,
                          lr_params = list(penalty = TRUE),
                          gmm_args = list(num_chains = 3),
                          model_saved = TRUE)
# model stability test
test_result <- stab_test_run(dat, dat_nm, model_result)
test_result$template_chosen # template chosen
quantile(test_result$result$Accuracy) # Accuracy
quantile(test_result$result$AUC) # AUC

# - demo 5 : Multiclass classification. LM-GMM achieves better result than ML.
dat_nm <- "fog_presence"
dat <- read.csv(sprintf("inst/extdata/real/%s.csv", dat_nm), header = T)
model_result <- model_run(dat, dat_nm = dat_nm, gmm_args = list(num_chains = 3), model_saved = TRUE)
# model stability test
test_result <- stab_test_run(dat, dat_nm, model_result)
test_result$template_chosen # template chosen
quantile(test_result$result$Accuracy) # Accuracy
quantile(test_result$result$AUC) # AUC

# - demo 6 : Multiclass classification. LM-GMM achieves better result than ML.
dat_nm <- "sale"
dat <- read.csv(sprintf("inst/extdata/real/%s.csv", dat_nm), header = T)
model_result <- model_run(dat, dat_nm = dat_nm,
                          gmm_args = list(num_chains = 3),
                          model_saved = TRUE, opt_metric = "AUC")
# model stability test
test_result <- stab_test_run(dat, dat_nm, model_result)
test_result$template_chosen # template chosen
quantile(test_result$result$Accuracy) # Accuracy
quantile(test_result$result$AUC) # AUC

# Simulation part
# *** It will take hours to run OctoCore model on all simulated datasets.
# *** Please check the results in getwd()/result/synthetic directory.
simu_types <- c("heavytailed", "moon", "multimodal", "normal", "spiral")
for (simu_type in simu_types) {
  print(simu_type)
  ptm <- proc.time()
  dat_nms <- dir(sprintf("inst/extdata/synthetic/%s", simu_type)) # all of the synthetic data
  simu_result <- data.frame() # record stability test results
  for (dat_nm in dat_nms) {
    print(dat_nm)
    dat_nm <- gsub("dat_|\\.csv", "", dat_nm)
    temp <- strsplit(dat_nm, "_")[[1]]
    Distribution <- temp[1]
    Categories <- gsub("K", "", temp[2])
    Samples <- gsub("N", "", temp[3])
    Variables <- gsub("d", "", temp[4])
    dat <- read.csv(sprintf("inst/extdata/synthetic/%s/dat_%s.csv", simu_type, dat_nm), header = T)
    model_result <- model_run(dat, dat_nm = gsub("dat_", "", dat_nm), model_saved = TRUE,
                              result_dir = sprintf("%s/result/synthetic/%s", getwd(), simu_type))
    # model stability test
    test_result <- stab_test_run(dat, dat_nm, model_result)
    simu_temp <- data.frame(test_result$result, Categories = Categories, Samples = Samples, Variables = Variables,
                            Distribution = Distribution)
    simu_result <- rbind(simu_result, simu_temp)
  }
  write.csv(simu_result, sprintf("result/synthetic/simu_result_%s.csv", simu_type), row.names = F)
  print(sprintf("run [%s] costs %s seconds.", simu_type, proc.time()[3] - ptm[3]))
  Sys.sleep(0.2)
}

# ! ------------------------- run demo ------------------------- ! #
