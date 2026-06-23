# Figure 10: BRCA CAT stability.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
CAT_DIR  <- Sys.getenv("ORIGINAL_CAT_DIR", unset = file.path(BASE_DIR, "CAT"))
CONC_DIR <- Sys.getenv("ORIGINAL_CAT_CONC_DIR", unset = file.path(CAT_DIR, "conc"))

NSAMPLE <- as.integer(strsplit(Sys.getenv("ORIGINAL_CAT_NSAMPLE", unset = "3,5,10,30,50"), ",")[[1]])

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig10_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(CONC_DIR)) {
  stop("CONC_DIR does not exist: ", CONC_DIR,
       "\nPlease check /public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE/CAT/conc.")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Package cowplot is required. Please install/load it in the R environment.")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  stop("Package ggpubr is required. Please install/load it in the R environment.")
}

METHOD_ORDER <- c(
  "ABSSeq", "DESeq", "DESeq2", "DSS", "edgeR.lrt", "edgeR.qlf",
  "NBPSeq", "NOISeq", "ROTS", "T.test", "voom", "Wilcoxon"
)

METHOD_COLS <- c(
  "ABSSeq"    = "#666666",
  "DESeq"     = "#6A3D9A",
  "DESeq2"    = "#CAB2D6",
  "DSS"       = "#E69F00",
  "edgeR.lrt" = "#A6CEE3",
  "edgeR.qlf" = "#56B4E9",
  "NBPSeq"    = "#009E73",
  "NOISeq"    = "#CC79A7",
  "ROTS"      = "#E31A1C",
  "T.test"    = "#FB9A99",
  "voom"      = "#11264F",
  "Wilcoxon"  = "#D55E00"
)

canonical_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x
}

read_conc_one_n <- function(n) {
  possible_files <- c(
    file.path(CONC_DIR, sprintf("%s-1vs2.rds", n)),
    file.path(CONC_DIR, sprintf("%s-1_vs_2.rds", n)),
    file.path(CONC_DIR, sprintf("%s-subset1_vs_subset2.rds", n)),
    file.path(CONC_DIR, sprintf("%s-conc.rds", n))
  )
  f <- possible_files[file.exists(possible_files)][1]

  if (is.na(f) || length(f) == 0) {
    warning("Missing conc file for n=", n, ". Tried: ", paste(possible_files, collapse = " ; "))
    return(data.frame())
  }

  x <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(x, "error")) {
    warning("Failed to read conc file: ", f, " | ", conditionMessage(x))
    return(data.frame())
  }

  if (is.data.frame(x)) {
    dat <- x
    if (!"iteration" %in% colnames(dat)) dat$iteration <- 1L
  } else if (is.list(x)) {
    parts <- list()
    for (i in seq_along(x)) {
      if (is.data.frame(x[[i]]) && nrow(x[[i]]) > 0) {
        tmp <- x[[i]]
        tmp$iteration <- i
        parts[[length(parts) + 1L]] <- tmp
      }
    }
    if (length(parts) == 0) {
      warning("No data.frame elements found in conc object: ", f)
      return(data.frame())
    }
    dat <- do.call(rbind, parts)
  } else {
    warning("Unsupported conc object type in: ", f)
    return(data.frame())
  }

  if (!"conc" %in% colnames(dat)) {
    alt <- intersect(c("Concordance", "concordance", "CAT", "cat"), colnames(dat))
    if (length(alt) > 0) dat$conc <- dat[[alt[1]]]
  }

  if (!"method1" %in% colnames(dat)) {
    alt <- intersect(c("Method", "method", "methods", "method_1"), colnames(dat))
    if (length(alt) > 0) dat$method1 <- dat[[alt[1]]]
  }

  if (!"method2" %in% colnames(dat)) {
    alt <- intersect(c("method_2"), colnames(dat))
    if (length(alt) > 0) {
      dat$method2 <- dat[[alt[1]]]
    } else {
      dat$method2 <- dat$method1
    }
  }

  if (!"rank" %in% colnames(dat)) dat$rank <- NA_integer_

  if (!all(c("method1", "method2", "conc") %in% colnames(dat))) {
    warning("Conc object lacks required columns method1/method2/conc in: ", f)
    warning("Available columns: ", paste(colnames(dat), collapse = ", "))
    return(data.frame())
  }

  dat$method1 <- canonical_method(dat$method1)
  dat$method2 <- canonical_method(dat$method2)
  dat$conc <- suppressWarnings(as.numeric(dat$conc))
  dat$nSample <- n
  dat$source_file <- f

  dat <- dat[dat$method1 == dat$method2, , drop = FALSE]
  dat <- dat[is.finite(dat$conc) & !is.na(dat$conc) & dat$conc >= 0 & dat$conc <= 1, , drop = FALSE]
  dat <- dat[dat$method1 %in% METHOD_ORDER, , drop = FALSE]

  dat
}

all_conc <- do.call(rbind, lapply(NSAMPLE, read_conc_one_n))

if (is.null(all_conc) || nrow(all_conc) == 0) {
  stop("No valid CAT concordance data were extracted from: ", CONC_DIR)
}

all_conc$method1 <- factor(canonical_method(all_conc$method1), levels = METHOD_ORDER)
all_conc$nSample <- factor(all_conc$nSample, levels = NSAMPLE)

write.csv(all_conc, file.path(OUT_DIR, "Fig10_data.csv"), row.names = FALSE)

cat_summary <- aggregate(conc ~ nSample + method1, all_conc, function(x) {
  c(
    n = length(x),
    min = min(x, na.rm = TRUE),
    q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q75 = unname(quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
})
cat_summary <- do.call(data.frame, cat_summary)
colnames(cat_summary) <- gsub("conc\\.", "", colnames(cat_summary))
write.csv(cat_summary, file.path(OUT_DIR, "Fig10_summary.csv"), row.names = FALSE)

count_table <- as.data.frame(table(all_conc$nSample, all_conc$method1))
colnames(count_table) <- c("nSample", "Method", "nRows")
write.csv(count_table, file.path(OUT_DIR, "Fig10_counts.csv"), row.names = FALSE)

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
  plot.margin = unit(c(0.10, 0.10, 0.10, 0.10), "cm")
)

make_one_panel <- function(n) {
  dat <- all_conc[as.character(all_conc$nSample) == as.character(n), , drop = FALSE]
  if (nrow(dat) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste0("n = ", n)))
  }

  dat$method_chr <- as.character(dat$method1)

  med <- aggregate(conc ~ method_chr, dat, median, na.rm = TRUE)
  med <- med[order(med$conc, decreasing = FALSE), , drop = FALSE]
  dat$method_plot <- factor(dat$method_chr, levels = med$method_chr)

  ggplot(dat, aes(x = method_plot, y = conc, color = method_chr)) +
    geom_boxplot(
      outlier.size = 0.65,
      outlier.alpha = 0.80,
      linewidth = 0.28,
      width = 0.65
    ) +
    coord_flip() +
    scale_color_manual(values = METHOD_COLS, limits = METHOD_ORDER, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    xlab("Method") +
    ylab("Concordance") +
    ggtitle(paste0("n = ", n)) +
    theme_minimal(base_size = 7) +
    plot_format +
    theme(legend.position = "none")
}

plotlist <- lapply(NSAMPLE, make_one_panel)

legend_dat <- all_conc[as.character(all_conc$nSample) == as.character(NSAMPLE[1]), , drop = FALSE]
if (nrow(legend_dat) == 0) legend_dat <- all_conc
legend_dat$method_chr <- as.character(legend_dat$method1)

g_legend <- ggplot(legend_dat, aes(x = method_chr, y = conc, color = method_chr)) +
  geom_boxplot(linewidth = 0.28) +
  scale_color_manual(values = METHOD_COLS, limits = METHOD_ORDER, drop = FALSE) +
  labs(color = "Method") +
  theme_minimal(base_size = 7) +
  plot_format +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(color = guide_legend(nrow = 4))

legend_grob <- ggpubr::get_legend(g_legend)
plotlist[[length(plotlist) + 1]] <- legend_grob

fig <- cowplot::plot_grid(
  plotlist = plotlist,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "")
)

pdf_file  <- file.path(OUT_DIR, "Fig10.pdf")
tiff_file <- file.path(OUT_DIR, "Fig10.tiff")
png_file  <- file.path(OUT_DIR, "Fig10.png")

ggplot2::ggsave(pdf_file, fig, width = 6.69, height = 7.88, units = "in")
ggplot2::ggsave(tiff_file, fig, width = 6.69, height = 7.88, units = "in", dpi = 600, compression = "lzw")
ggplot2::ggsave(png_file, fig, width = 6.69, height = 7.88, units = "in", dpi = 300)

