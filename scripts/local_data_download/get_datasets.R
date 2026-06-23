library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(data.table)
library(edgeR)
library(GEOquery)
library(GenomicRanges)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
DATASETS_DIR <- file.path(BASE_DIR, "datasets")
dir.create(DATASETS_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(BASE_DIR)

## Download and preprocess E-ENAD-34 dataset from EBI repository for further analysis.

acc <- c("E-ENAD-34")
exps <- ExpressionAtlas::getAtlasData(acc)
eset <- exps[[acc]]$rnaseq
se <- exps[[1]]$rnaseq
counts <- assays(se)$counts 
coldata <- colData(se)
counts <- counts[, coldata$AtlasAssayGroup %in% c("g2", "g5")]
group <- coldata[coldata$AtlasAssayGroup %in% c("g2", "g5"), "AtlasAssayGroup"]
group <- dplyr::case_when(
  group == "g2" ~ "neutrophil-non somker",
  group == "g5" ~ "neutrophil-somker"
)  
dge <- DGEList(counts = counts, group = group)
saveRDS(dge, file = file.path(DATASETS_DIR, sprintf("DGElist_%s.rds", acc)))

## Download and preprocess GSE150910 dataset from GEO repository for further analysis.
counts <- fread(file.path(DATASETS_DIR, "GSE150910_gene-level_count_file.csv.gz"))
counts <- as.data.frame(counts)
rownames(counts) <- counts$symbol
counts$symbol <- NULL
group <- sapply(strsplit(colnames(counts), "_"), "[[", 1)
keep <- group != "chp"
counts <- counts[, keep, drop = FALSE]
group <- group[keep]
dge <- DGEList(counts = counts, group = group)
saveRDS(dge, file = file.path(DATASETS_DIR, "DGElist_GSE150910.rds"))

## Download and preprocess TCGA-LUSC and TCGA-BRCA datasets from GDC repository for further analysis.

cancer_list <- c("TCGA-LUSC", "TCGA-BRCA")
for (name in cancer_list) {
  proj <- TCGAbiolinks:::getProjectSummary(name)
  # query
  query <- GDCquery(
    project       = name,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # download
  GDCdownload(
    query           = query,
    method          = "api",
    files.per.chunk = 5
  )
  
  # prepare
  dat <- GDCprepare(
    query         = query,
    save          = FALSE,
    remove.files.prepared = FALSE
  )
  
  df <- as.data.frame(rowRanges(dat))
  gene_expr <- as.data.frame(assay(dat, i = "unstranded"))
  rt <- limma::avereps(gene_expr, ID = df$gene_name)
  rt <- as.data.frame(rt)
  
  dimnames <- list(rownames(rt), colnames(rt))
  rt <- matrix(as.numeric(as.matrix(rt)), nrow = nrow(rt), dimnames = dimnames)
  rt <- as.data.frame(rt)
  rt_N <- rt[, substr(colnames(rt), 14, 15) > 10, drop = F]
  rt_T <- rt[, substr(colnames(rt), 14, 15) < 10, drop = F]
  rt <- cbind(rt_N, rt_T)
  column_names <- data.frame(sample = names(rt))
  column_names$group <- ifelse(substr(column_names$sample, 14, 15) > 10,  "Normal", "Cancer")
  dge <- DGEList(counts = rt, group = column_names$group)
  saveRDS(dge, file = file.path(DATASETS_DIR,
    sprintf("DGElist_%s.rds", gsub("TCGA-", "", name))))
}

## Filter out outlier samples based on t-SNE distances using a 3 IQR threshold.
files <- list.files(DATASETS_DIR, pattern = "^DGElist_.*\\.rds$", full.names = TRUE)
library(Rtsne)
for (file in files) {
  d <- readRDS(file)
  lsize <- colSums(d$counts)
  datas <- as.data.frame(log2(t(t(d$counts) / lsize * 1e6) + 1))
  set.seed(1234)
  tsne <- Rtsne(t(datas), dims = 2, perplexity = 30, check_duplicates = FALSE)
  
  tsne_data <- data.frame(
    Sample = colnames(datas), 
    group = d$samples$group,
    TSNE1 = tsne$Y[, 1], 
    TSNE2 = tsne$Y[, 2]
  )
  
  rownames(tsne_data) <- tsne_data$Sample
  cons <- unique(tsne_data$group)
  outlier_samples <- lapply(cons, function(con) {
    dist_matrix <- dist(tsne_data[tsne_data$group == con, c("TSNE1", "TSNE2")])
    avg_distance <- apply(as.matrix(dist_matrix), 1, mean)
    threshold <- quantile(avg_distance, 0.75) + 3 * IQR(avg_distance) 
    outlier_samples <- names(which(avg_distance > threshold))
  })
  badreps <- do.call(c, outlier_samples)
  
  counts <- d$counts[, colnames(d$counts) != badreps]
  group <- d$samples$group[colnames(d$counts) != badreps]
  dge <- DGEList(counts = counts, group = group)
  saveRDS(dge, file = file)
}
