library(parallel)

getDE_split <- function(y) {
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
  res <- list("p" = p, "label" = label)
}

getp <- function(file) {
  dat <- readRDS(file)
  datp <- mclapply(1:length(dat), function(n) getDE_split(dat[[n]]), mc.cores = 10)
  saveRDS(datp, file = sprintf("./datp/%s.rds", gsub(".rds", "", basename(file))))
}

getAUC <- function(file) {
  dat <- readRDS(file)
  AUC <- lapply(dat, function(x) {
    p <- x$p
    label <- x$label
    p[is.na(p)] <- 1
    auc <- lapply(p, function(col) {
      auc <- roc(label, col, levels = c(0, 1), direction = ">")$auc
    })
    
    auc <- unlist(auc)
    auc <- data.frame(auc = auc, methods = names(auc))
  })
  saveRDS(AUC, file = gsub("datp", "auc", file))
}
