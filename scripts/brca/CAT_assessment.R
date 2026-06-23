# CAT assessment.

library(edgeR)
library(tidyverse)
library(parallel)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
WORK_DIR <- file.path(BASE_DIR, "CAT")
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

dirs <- c("./counts", "./datp", "./conc")
lapply(dirs, dir.create)

B <- 10
nSample <- c(3, 5, 10, 30, 50)

dge <- readRDS(file.path(BASE_DIR, "datasets", "DGElist_BRCA.rds"))
set.seed(123)
myseed <- sample(1:99999, 10)

for (n in nSample) {
  
  subset1 <- lapply(myseed, function(seed) {
    set.seed(seed)
    index.a <- sample(which(dge$samples$group == "Cancer"), 2 * n)
    index.b <- sample(which(dge$samples$group == "Normal"), 2 * n)
    counts <- dge$counts[, c(head(index.a, n = n), head(index.b, n = n))]
    group <- c(rep("a", n), rep("b", n))
    d <- DGEList(counts = counts, group = group)
    return(d)
  })
  saveRDS(subset1, file = sprintf("./counts/%s-subset1.rds", n))
  
  subset2 <- lapply(myseed, function(seed) {
    set.seed(seed)
    index.a <- sample(which(dge$samples$group == "Cancer"), 2 * n)
    index.b <- sample(which(dge$samples$group == "Normal"), 2 * n)
    counts <- dge$counts[, c(tail(index.a, n = n), tail(index.b, n = n))]
    group <- c(rep("a", n), rep("b", n))
    d <- DGEList(counts = counts, group = group)
    return(d)
  })
  saveRDS(subset2, file = sprintf("./counts/%s-subset2.rds", n))
  
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

source(file.path(BASE_DIR, "utils", "get_conc_function.R"))

files <- list.files("./datp", full.names = TRUE, pattern = ".rds")
lapply(files, get_conc)

lapply(nSample, get_conc_w)
