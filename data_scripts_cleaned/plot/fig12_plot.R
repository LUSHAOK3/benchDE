# Figure 12: BRCA false-positive-rate control.

set.seed(8848)

BASE_DIR   <- "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE"
FPR_SUBDIR <- "FPR"
FPR_DIR    <- file.path(BASE_DIR, FPR_SUBDIR)

DATP_DIR_CANDIDATES <- c(
  file.path(FPR_DIR, "datp")
)
DATP_DIR <- DATP_DIR_CANDIDATES[file.exists(DATP_DIR_CANDIDATES)][1]

if (is.na(DATP_DIR) || length(DATP_DIR) == 0) {
  stop("No valid DATP_DIR found under: ", FPR_DIR,
       "\nExpected one of: ", paste(DATP_DIR_CANDIDATES, collapse = " ; "))
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR   <- file.path(BASE_DIR, "figures", paste0("Fig12_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
})

NSAMPLE_ORDER <- c(3, 5, 10, 30, 50)
FDR_CUTOFFS   <- c(0.01, 0.05, 0.10)

METHOD_ORDER <- c(
  "ABSSeq", "DESeq", "NOISeq", "ROTS", "T.test", "DSS",
  "Wilcoxon", "voom", "DESeq2", "edgeR.lrt", "edgeR.qlf", "NBPSeq"
)

METHOD_COLS <- c(
  "ABSSeq"    = "#666666",
  "DESeq"     = "#6A3D9A",
  "NOISeq"    = "#cc79a7",
  "ROTS"      = "#E31A1C",
  "T.test"    = "#FB9A99",
  "DSS"       = "#e69f00",
  "Wilcoxon"  = "#d55e00",
  "voom"      = "#11264f",
  "DESeq2"    = "#CAB2D6",
  "edgeR.lrt" = "#A6CEE3",
  "edgeR.qlf" = "#56b4e9",
  "NBPSeq"    = "#009e73"
)

plot_format <- theme(
  axis.line        = element_line(color = "black", linewidth = 0.25),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.25),
  axis.title       = element_text(size = 8),
  axis.text        = element_text(size = 7, color = "black"),
  plot.title       = element_text(size = 8, hjust = 0.5),
  legend.title     = element_text(size = 8),
  legend.text      = element_text(size = 7),
  plot.margin      = unit(c(0.12, 0.12, 0.12, 0.12), "cm")
)

canonical_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x
}

parse_meta_from_filename <- function(f) {
  bn <- basename(f)
  m <- regexec("^([0-9]+)-(.+)\\.rds$", bn)
  reg <- regmatches(bn, m)[[1]]
  if (length(reg) == 3) {
    return(list(nSample = as.integer(reg[2]), Dataset = reg[3]))
  }
  list(nSample = NA_integer_, Dataset = NA_character_)
}

contains_try_error <- function(x) {
  if (inherits(x, "try-error")) return(TRUE)
  if (inherits(x, "error")) return(TRUE)
  if (is.list(x)) return(any(vapply(x, contains_try_error, logical(1))))
  FALSE
}

valid_datp <- function(x) {
  if (inherits(x, "error")) return(FALSE)
  if (!is.list(x) || length(x) == 0) return(FALSE)
  if (contains_try_error(x)) return(FALSE)
  all(vapply(x, is.data.frame, logical(1)))
}

read_rds_safe <- function(f) {
  tryCatch(readRDS(f), error = function(e) e)
}

extract_fpr_from_one_file <- function(file) {
  meta <- parse_meta_from_filename(file)
  if (is.na(meta$nSample) || is.na(meta$Dataset)) {
    warning("Skip file with unrecognized filename format: ", file)
    return(data.frame())
  }

  datp <- read_rds_safe(file)
  if (!valid_datp(datp)) {
    warning("Invalid datp skipped: ", file)
    return(data.frame())
  }

  rows <- list()
  idx  <- 1L
  for (rep_i in seq_along(datp)) {
    pmat <- datp[[rep_i]]
    colnames(pmat) <- canonical_method(colnames(pmat))
    methods_here <- intersect(METHOD_ORDER, colnames(pmat))
    if (length(methods_here) == 0) next

    for (m in methods_here) {
      pv <- suppressWarnings(as.numeric(pmat[[m]]))
      ok <- is.finite(pv) & !is.na(pv) & pv >= 0 & pv <= 1
      if (!any(ok)) next
      padj <- p.adjust(pv[ok], method = "BH")

      for (cutoff in FDR_CUTOFFS) {
        rows[[idx]] <- data.frame(
          Dataset    = meta$Dataset,
          nSample    = meta$nSample,
          replicate  = rep_i,
          Method     = m,
          cutoff     = as.character(cutoff),
          cutoff_num = cutoff,
          FPR        = mean(padj < cutoff),
          nGenes     = sum(ok),
          source_file = file,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

all_files <- list.files(DATP_DIR, pattern = "\\.rds$", full.names = TRUE)
all_fpr <- do.call(rbind, lapply(all_files, extract_fpr_from_one_file))

if (is.null(all_fpr) || nrow(all_fpr) == 0) {
  stop("No FPR data extracted. Please check DATP_DIR and filename format.")
}

all_fpr$Method <- canonical_method(all_fpr$Method)
all_fpr <- all_fpr[all_fpr$Method %in% METHOD_ORDER, , drop = FALSE]
all_fpr$nSample <- as.integer(all_fpr$nSample)
all_fpr <- all_fpr[all_fpr$nSample %in% NSAMPLE_ORDER, , drop = FALSE]
all_fpr$Method <- factor(all_fpr$Method, levels = METHOD_ORDER)
all_fpr$cutoff <- factor(all_fpr$cutoff, levels = as.character(FDR_CUTOFFS))

write.csv(all_fpr, file.path(OUT_DIR, "Fig12_data.csv"), row.names = FALSE)

fpr_summary <- aggregate(FPR ~ Dataset + nSample + Method + cutoff + cutoff_num, all_fpr, function(x) {
  c(
    n      = length(x),
    min    = min(x),
    q25    = unname(quantile(x, 0.25)),
    median = median(x),
    mean   = mean(x),
    q75    = unname(quantile(x, 0.75)),
    max    = max(x),
    sd     = sd(x),
    n_zero = sum(x == 0)
  )
})
fpr_summary <- do.call(data.frame, fpr_summary)
colnames(fpr_summary) <- gsub("FPR\\.", "", colnames(fpr_summary))
write.csv(fpr_summary, file.path(OUT_DIR, "Fig12_summary.csv"), row.names = FALSE)

make_one_panel <- function(s) {
  dat <- all_fpr[all_fpr$nSample == s, , drop = FALSE]
  dat$Method <- droplevels(dat$Method)
  cols_here <- METHOD_COLS[levels(dat$Method)]

  p <- ggplot(dat, aes(x = cutoff, y = FPR, color = Method)) +
    geom_boxplot(
      aes(group = interaction(cutoff, Method)),
      position = position_dodge2(width = 0.82, preserve = "single", padding = 0.18),
      width = 0.55,
      linewidth = 0.30,
      outlier.size = 0.45,
      outlier.alpha = 0.65
    ) +
    coord_trans(y = "sqrt") +
    annotate("segment", x = 0.60, xend = 1.40, y = 0.01, yend = 0.01,
             color = "red", linewidth = 0.28, linetype = 2) +
    annotate("segment", x = 1.60, xend = 2.40, y = 0.05, yend = 0.05,
             color = "red", linewidth = 0.28, linetype = 2) +
    annotate("segment", x = 2.60, xend = 3.40, y = 0.10, yend = 0.10,
             color = "red", linewidth = 0.28, linetype = 2) +
    scale_color_manual(values = cols_here, drop = FALSE) +
    xlab("Nominal False Discovery Rate") +
    ylab("False Positive Rate") +
    ggtitle(paste0("n = ", s)) +
    plot_format +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 7, angle = 0, vjust = 0.9))
  p
}

plotlist <- lapply(NSAMPLE_ORDER, make_one_panel)

legend_data <- all_fpr[all_fpr$nSample == NSAMPLE_ORDER[1], , drop = FALSE]
legend_data$Method <- droplevels(legend_data$Method)
legend_cols <- METHOD_COLS[levels(legend_data$Method)]

g_legend <- ggplot(legend_data, aes(x = cutoff, y = FPR, color = Method)) +
  geom_boxplot(
    aes(group = interaction(cutoff, Method)),
    position = position_dodge2(width = 0.82, preserve = "single", padding = 0.18),
    width = 0.55
  ) +
  scale_color_manual(values = legend_cols, drop = FALSE) +
  plot_format +
  theme(legend.position = "right", legend.box = "vertical") +
  guides(color = guide_legend(ncol = 3))

legend_grob <- ggpubr::get_legend(g_legend)
plotlist[[length(plotlist) + 1]] <- legend_grob

fig <- cowplot::plot_grid(
  plotlist = plotlist,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "")
)

pdf_file  <- file.path(OUT_DIR, "Fig12.pdf")
tiff_file <- file.path(OUT_DIR, "Fig12.tiff")
png_file  <- file.path(OUT_DIR, "Fig12.png")

ggsave(pdf_file,  fig, width = 7.6, height = 8.6, units = "in")
ggsave(tiff_file, fig, width = 7.6, height = 8.6, units = "in", dpi = 600, compression = "lzw")
ggsave(png_file,  fig, width = 7.6, height = 8.6, units = "in", dpi = 300)

