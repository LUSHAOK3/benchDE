# Negative-binomial outlier simulation.

library(parallel)
library(pROC)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
UTILS_DIR <- file.path(BASE_DIR, "utils")
DATASETS_DIR <- file.path(BASE_DIR, "datasets")
WORK_DIR <- file.path(BASE_DIR, "sim2")

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

source(file.path(UTILS_DIR, "sim2_function.R"))

lapply(c("counts", "datp", "auc"), dir.create, recursive = TRUE,
       showWarnings = FALSE)

B <- as.integer(Sys.getenv("BENCHDE_B", unset = "100"))
CORES <- as.integer(Sys.getenv("BENCHDE_CORES", unset = "10"))
SEED <- as.integer(Sys.getenv("BENCHDE_SEED", unset = "8848"))
NSAMPLE <- c(3, 5, 10, 30, 50)
OUTLIER_LEVELS <- c(0.1, 0.5, 1, 2, 5)
files <- file.path(
  DATASETS_DIR,
  c("DGElist_E-ENAD-34.rds", "DGElist_GSE150910.rds", "DGElist_LUSC.rds")
)

for (f in files) {
  dge <- readRDS(f)
  param <- get_param(dge)
  for (level in OUTLIER_LEVELS) {
    for (s in NSAMPLE) {
      res <- parallel::mclapply(
        X = seq_len(B), FUN = simulateOneSplit, s = s, rseed = SEED,
        RO.prop = level, param = param, mc.cores = CORES
      )
      saveRDS(res, file.path("counts", sprintf("%s-%s-%s", level, s,
        sub("DGElist_", "", basename(f)))))
    }
  }
}

source(file.path(UTILS_DIR, "runDE.R"))
source(file.path(UTILS_DIR, "sim_DEanalysis_auc.R"))

lapply(list.files("counts", full.names = TRUE, pattern = "\\.rds$"), getp)
parallel::mclapply(
  list.files("datp", full.names = TRUE, pattern = "\\.rds$"),
  getAUC,
  mc.cores = CORES
)
