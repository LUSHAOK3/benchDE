# poisson sim2 assessment resume safe.

library(parallel)
library(pROC)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
UTILS_DIR <- file.path(BASE_DIR, "utils")
DATASETS_DIR <- file.path(BASE_DIR, "datasets")
WORK_DIR <- file.path(BASE_DIR, "poisson_sim2")

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

source(file.path(UTILS_DIR, "poisson_sim2_function.R"))
source(file.path(UTILS_DIR, "runDE.R"))

dirs <- c("./counts", "./datp", "./auc")
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("BENCHDE_B", unset = "100"))
mc.cores <- as.integer(Sys.getenv("BENCHDE_CORES", unset = "4"))
SEED <- as.integer(Sys.getenv("BENCHDE_SEED", unset = "8848"))

nsample <- c(3, 5, 10, 30, 50)
RO.props <- c(0.1, 0.5, 1, 2, 5)

files <- c(
  file.path(DATASETS_DIR, "DGElist_E-ENAD-34.rds"),
  file.path(DATASETS_DIR, "DGElist_GSE150910.rds"),
  file.path(DATASETS_DIR, "DGElist_LUSC.rds")
)

for (f in files) {
  dge <- readRDS(f)
  param <- get_param(dge)

  for (RO.prop in RO.props) {
    for (s in nsample) {
      ro_label <- format(RO.prop, trim = TRUE, scientific = FALSE)
      out_file <- sprintf("./counts/%s-%s-%s", ro_label, s, gsub("DGElist_", "", basename(f)))

      if (file.exists(out_file)) {
        next
      }

      res <- parallel::mclapply(
        X = 1:B,
        s = s,
        FUN = simulateOneSplit,
        rseed = SEED,
        RO.prop = RO.prop,
        param = param,
        mc.cores = mc.cores
      )
      saveRDS(res, file = out_file)
      rm(res)
      gc()
    }
  }
  rm(dge, param)
  gc()
}

getDE_split_safe <- function(y) {
  label <- as.numeric(y$truth)
  dge <- y$dge
  keep <- rowSums(dge$counts) > 0
  label <- label[keep]
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  plis <- lapply(algos, function(f) f(dge))
  gene <- rownames(plis[[1]])
  plis <- lapply(plis, function(x) x[gene, , drop = FALSE])
  p <- data.frame(lapply(plis, '[[', "PValue"), row.names = gene)
  p[is.na(p)] <- 1
  list("p" = p, "label" = label)
}

getp_safe <- function(file) {
  out_file <- sprintf("./datp/%s.rds", gsub(".rds", "", basename(file)))
  if (file.exists(out_file)) {
    return(invisible(out_file))
  }
  dat <- readRDS(file)
  datp <- parallel::mclapply(
    seq_along(dat),
    function(n) getDE_split_safe(dat[[n]]),
    mc.cores = mc.cores
  )
  saveRDS(datp, file = out_file)
  rm(dat, datp)
  gc()
  invisible(out_file)
}

count_files <- list.files("./counts", full.names = TRUE, pattern = ".rds$")
for (cf in count_files) {
  getp_safe(cf)
}

getAUC_safe <- function(file) {
  out_file <- gsub("datp", "auc", file)
  if (file.exists(out_file)) {
    return(invisible(out_file))
  }
  dat <- readRDS(file)
  AUC <- lapply(dat, function(x) {
    p <- x$p
    label <- x$label
    p[is.na(p)] <- 1
    auc <- lapply(p, function(col) {
      pROC::roc(label, col, levels = c(0, 1), direction = ">", quiet = TRUE)$auc
    })
    auc <- unlist(auc)
    data.frame(auc = auc, methods = names(auc))
  })
  saveRDS(AUC, file = out_file)
  rm(dat, AUC)
  gc()
  invisible(out_file)
}

datp_files <- list.files("./datp", full.names = TRUE, pattern = ".rds$")
for (df in datp_files) {
  getAUC_safe(df)
}

