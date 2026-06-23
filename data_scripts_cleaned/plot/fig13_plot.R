# Figure 13: SCLC false-positive-rate control.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
DATASET <- Sys.getenv("SCLC_DATASET_LABEL", unset = "SCLC")
NSAMPLE <- as.integer(strsplit(Sys.getenv("SCLC_PLOT_NSAMPLE", unset = "3,5,10,30,50"), ",")[[1]])
FDR_CUTOFFS <- as.numeric(strsplit(Sys.getenv("SCLC_FPR_CUTOFFS", unset = "0.01,0.05,0.1"), ",")[[1]])
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig13_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
})

METHODS <- c("ABSSeq", "DESeq", "NOISeq", "ROTS", "T.test", "DSS",
             "Wilcoxon", "voom", "DESeq2", "edgeR.lrt", "edgeR.qlf", "NBPSeq")
METHODS_ORIGINAL <- c("DESeq2", "DESeq", "edgeR.lrt", "edgeR.qlf", "DSS",
                      "voom", "NBPSeq", "Wilcoxon", "T.test", "NOISeq", "ROTS", "ABSSeq")
COLS <- c("#666666", "#6A3D9A", "#cc79a7", "#E31A1C", "#FB9A99", "#e69f00",
          "#d55e00", "#11264f", "#CAB2D6", "#A6CEE3", "#56b4e9", "#009e73")
names(COLS) <- METHODS

canonical_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x
}

choose_file <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NA_character_)
  paths[1]
}

plot_format <- theme(
  axis.line = element_line(color = "black", linewidth = 0.25),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6, color = "black"),
  plot.title = element_text(size = 7, hjust = 0.5),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 6),
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
)

get_fpr_long_one_file <- function(file, nsample) {
  datp <- readRDS(file)
  rows <- list()
  idx <- 1L

  for (sim_i in seq_along(datp)) {
    p <- datp[[sim_i]]
    colnames(p) <- canonical_method(colnames(p))
    available <- intersect(METHODS_ORIGINAL, colnames(p))

    if (length(available) == 0) {
      warning("No known methods in file: ", file)
      next
    }

    p <- p[, available, drop = FALSE]
    for (method in available) {
      pv <- suppressWarnings(as.numeric(p[[method]]))
      ok <- is.finite(pv) & !is.na(pv) & pv >= 0 & pv <= 1
      if (!any(ok)) next
      padj <- p.adjust(pv[ok], method = "BH")

      for (cutoff in FDR_CUTOFFS) {
        rows[[idx]] <- data.frame(
          Dataset = DATASET,
          nSample = nsample,
          sim = sim_i,
          Method = canonical_method(method),
          cutoff = as.character(cutoff),
          cutoff_num = cutoff,
          FPR = mean(padj < cutoff),
          nGenes = sum(ok),
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

all_fpr <- data.frame()

for (s in NSAMPLE) {
  file <- choose_file(c(
    file.path(BASE_DIR, "FPR_SCLC", "datp", sprintf("%s-SCLC.rds", s)),
    file.path(BASE_DIR, "FPR_SCLC", "datp_method_combined", sprintf("%s-SCLC.rds", s))
  ))
  if (is.na(file)) {
    warning("Missing FPR datp for n=", s)
    next
  }
  all_fpr <- rbind(all_fpr, get_fpr_long_one_file(file, s))
}

if (nrow(all_fpr) == 0) stop("No FPR data generated.")

all_fpr$Method <- canonical_method(all_fpr$Method)
method_levels <- METHODS[METHODS %in% unique(all_fpr$Method)]
all_fpr$Method <- factor(all_fpr$Method, levels = method_levels)
all_fpr$cutoff <- factor(all_fpr$cutoff, levels = as.character(FDR_CUTOFFS))

fpr_summary <- aggregate(FPR ~ Dataset + nSample + Method + cutoff + cutoff_num, all_fpr, function(x) {
  c(mean = mean(x), median = median(x), sd = sd(x), q25 = unname(quantile(x, 0.25)), q75 = unname(quantile(x, 0.75)))
})
fpr_summary <- do.call(data.frame, fpr_summary)
colnames(fpr_summary) <- gsub("FPR\\.", "", colnames(fpr_summary))

write.csv(all_fpr, file.path(OUT_DIR, "Fig13_data.csv"), row.names = FALSE)
write.csv(fpr_summary, file.path(OUT_DIR, "Fig13_summary.csv"), row.names = FALSE)

get_fpr_plot <- function(s) {
  dat <- all_fpr[all_fpr$nSample == s, , drop = FALSE]
  dat$Method <- droplevels(dat$Method)

  p <- ggplot(dat, aes(color = Method)) +
    geom_boxplot(data = dat[dat$cutoff == as.character(FDR_CUTOFFS[1]), ],
                 aes(x = as.character(FDR_CUTOFFS[1]), y = FPR),
                 outlier.size = 0.7, linewidth = 0.25, width = 0.65) +
    geom_boxplot(data = dat[dat$cutoff == as.character(FDR_CUTOFFS[2]), ],
                 aes(x = as.character(FDR_CUTOFFS[2]), y = FPR),
                 outlier.size = 0.7, linewidth = 0.25, width = 0.65) +
    geom_boxplot(data = dat[dat$cutoff == as.character(FDR_CUTOFFS[3]), ],
                 aes(x = as.character(FDR_CUTOFFS[3]), y = FPR),
                 outlier.size = 0.7, linewidth = 0.25, width = 0.65) +
    coord_trans(y = "sqrt") +
    geom_segment(aes(x = 1 - 0.5, xend = 1 + 0.5, y = FDR_CUTOFFS[1], yend = FDR_CUTOFFS[1]),
                 color = "red", linewidth = 0.25, lty = 2) +
    geom_segment(aes(x = 2 - 0.5, xend = 2 + 0.5, y = FDR_CUTOFFS[2], yend = FDR_CUTOFFS[2]),
                 color = "red", linewidth = 0.25, lty = 2) +
    geom_segment(aes(x = 3 - 0.5, xend = 3 + 0.5, y = FDR_CUTOFFS[3], yend = FDR_CUTOFFS[3]),
                 color = "red", linewidth = 0.25, lty = 2) +
    xlab("Nominal False Discovery Rate") +
    ylab("False Positive Rate") +
    ggtitle(paste0("n = ", s)) +
    scale_color_manual(values = COLS, drop = FALSE) +
    plot_format +
    theme(legend.position = "none")
  p
}

plotlist <- lapply(NSAMPLE, get_fpr_plot)

legend_data <- all_fpr[all_fpr$nSample == NSAMPLE[1], , drop = FALSE]
g_legend <- ggplot(legend_data, aes(x = cutoff, y = FPR, color = Method)) +
  geom_boxplot() +
  scale_color_manual(values = COLS, drop = FALSE) +
  plot_format +
  theme(legend.position = "right", legend.box = "vertical") +
  guides(color = guide_legend(ncol = 3))
legend <- ggpubr::get_legend(g_legend)
plotlist[[length(plotlist) + 1]] <- legend

pdf_file <- file.path(OUT_DIR, "Fig13.pdf")
tiff_file <- file.path(OUT_DIR, "Fig13.tiff")
ggsave(pdf_file, cowplot::plot_grid(plotlist = plotlist, ncol = 2, labels = c("A", "B", "C", "D", "E", "")),
       width = 6.69, height = 7.88, units = "in")
ggsave(tiff_file, cowplot::plot_grid(plotlist = plotlist, ncol = 2, labels = c("A", "B", "C", "D", "E", "")),
       width = 6.69, height = 7.88, units = "in", dpi = 600, compression = "lzw")

