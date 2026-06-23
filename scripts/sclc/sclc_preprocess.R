# Prepare the SCLC count data for the FPR and CAT benchmarks.

library(edgeR)
library(data.table)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())
COUNTS_FILE <- Sys.getenv(
  "SCLC_COUNTS",
  unset = file.path(BASE_DIR, "datasets", "count_ensg_STAR.txt.gz")
)
CLINICAL_FILE <- Sys.getenv(
  "SCLC_CLINICAL",
  unset = file.path(BASE_DIR, "datasets", "Cell2023_SCLC_data_clinical.txt")
)
DATASET_DIR <- file.path(BASE_DIR, "datasets")

B_FPR <- as.integer(Sys.getenv("SCLC_FPR_B", unset = "100"))
B_CAT <- as.integer(Sys.getenv("SCLC_CAT_B", unset = "10"))
FPR_MODE <- Sys.getenv("SCLC_FPR_MODE", unset = "normal_only")

raw_counts <- data.table::fread(
  COUNTS_FILE,
  sep = "auto",
  data.table = FALSE,
  check.names = FALSE
)

gene_id <- as.character(raw_counts[[1]])
raw_counts[[1]] <- NULL
counts <- as.matrix(raw_counts)
storage.mode(counts) <- "numeric"
rownames(counts) <- gene_id

normal_samples <- colnames(counts)[grepl("^N", colnames(counts))]
cancer_samples <- colnames(counts)[grepl("^T", colnames(counts))]
counts <- counts[, c(normal_samples, cancer_samples), drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

if (anyDuplicated(rownames(counts))) {
  counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
}

group <- ifelse(grepl("^N", colnames(counts)), "Normal", "Cancer")
dge <- edgeR::DGEList(counts = counts, group = group)

sample_info <- data.frame(
  sample = colnames(counts),
  group = group,
  patient_id = sub("^[NT]", "", colnames(counts)),
  stringsAsFactors = FALSE
)

clinical <- data.table::fread(
  CLINICAL_FILE,
  sep = "\t",
  data.table = FALSE,
  check.names = FALSE
)

if ("Sample.ID" %in% names(clinical)) {
  clinical$Sample.ID <- as.character(clinical$Sample.ID)
  sample_info <- merge(
    sample_info,
    clinical,
    by.x = "patient_id",
    by.y = "Sample.ID",
    all.x = TRUE,
    sort = FALSE
  )
  rownames(sample_info) <- sample_info$sample
  sample_info <- sample_info[colnames(counts), , drop = FALSE]
}

attr(dge, "sample_info") <- sample_info
dir.create(DATASET_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(dge, file.path(DATASET_DIR, "DGElist_SCLC_from_newdata.rds"))
