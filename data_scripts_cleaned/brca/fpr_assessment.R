# fpr assessment.

library(parallel)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
WORK_DIR <- file.path(BASE_DIR, "FPR")
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

dirs <- c("./counts", "./datp")
lapply(dirs, dir.create)

B <- 100
nsample <- c(3,5,10,30,50)

dge <- readRDS(file.path(BASE_DIR, "datasets", "DGElist_BRCA.rds"))
counts <- dge$counts

for (s in nsample) {
  set.seed(123)
  myseed <- sample(1:2000, B)
  counts.list <- list()
  
  for (i in 1:length(myseed)) {
    set.seed(myseed[i])
    inds <- sample(colnames(counts), 2 * s, replace = FALSE)
    group <- rep(c(1, 2), each = s)
    dat <- counts[, inds]
    d <- DGEList(counts = dat, group = group)
    counts.list[[i]] <- d
  }
  
  saveRDS(counts.list, sprintf("./counts/%s-BRCA.rds", s))
}

source(file.path(BASE_DIR, "utils", "runDE.R"))

getDE_split <- function(y) {
  keep <- rowSums(y$counts >= 5) >= 2 
  y <- y[keep, , keep.lib.sizes = FALSE]
  plis <- lapply(algos, function(f) f(y))
  gene <- rownames(plis[[1]])
  plis <- lapply(plis, function(x) x[gene, , drop = FALSE])
  p <- data.frame(lapply(plis, '[[', "PValue"), row.names = gene)
  p[is.na(p)] <- 1
  return(p)
}

getp <- function(file) {
  dat <- readRDS(file)
  datp <- mclapply(1:length(dat), function(n) getDE_split(dat[[n]]), mc.cores = 10)
  saveRDS(datp, file = sprintf("./datp/%s.rds", gsub(".rds", "", basename(file))))
}

files <- list.files("./counts", full.names = TRUE, pattern = ".rds")
lapply(files, getp)
