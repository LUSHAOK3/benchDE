# poisson sim3 assessment.

library(parallel)
library(pROC)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
UTILS_DIR <- file.path(BASE_DIR, "utils")
DATASETS_DIR <- file.path(BASE_DIR, "datasets")
WORK_DIR <- file.path(BASE_DIR, "poisson_sim3")

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

source(file.path(UTILS_DIR, "poisson_sim3_function.R"))

dirs <- c("./counts", "./datp", "./auc")
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("BENCHDE_B", unset = "100"))
mc.cores <- as.integer(Sys.getenv("BENCHDE_CORES", unset = "10"))
SEED <- as.integer(Sys.getenv("BENCHDE_SEED", unset = "8848"))

nsample <- c(3, 5, 10, 30, 50)
Noise_levels <- c(0, 0.2, 0.4, 0.6, 0.8)

files <- c(file.path(DATASETS_DIR, "DGElist_E-ENAD-34.rds"))

for (f in files) {
  dge <- readRDS(f) 
  param <- get_param(dge)
  for (Noise in Noise_levels) {
    for (s in nsample) {
      res <- parallel::mclapply(
        X = 1:B,
        s = s,
        FUN = simulateOneSplit,
        rseed = SEED,
        Noise = Noise,
        param = param,
        mc.cores = mc.cores
      )
      noise_label <- format(Noise, trim = TRUE, scientific = FALSE)
      saveRDS(res, file = sprintf("./counts/%s-%s-%s", noise_label, s, gsub("DGElist_", "", basename(f))))
    }
  }
}

source(file.path(UTILS_DIR, "runDE.R"))
source(file.path(UTILS_DIR, "sim_DEanalysis_auc.R"))

files <- list.files("./counts", full.names = TRUE, pattern = ".rds$")
lapply(files, getp)

files <- list.files("./datp", full.names = TRUE, pattern = ".rds$")
parallel::mclapply(files, getAUC, mc.cores = mc.cores)

