# Figure 2: negative-binomial standard-simulation AUC.

set.seed(8848)

BASE_DIR <- "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE"
AUC_DIR  <- "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE/sim1/auc"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig2_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(AUC_DIR)) stop("AUC_DIR does not exist: ", AUC_DIR)

suppressPackageStartupMessages(library(ggplot2))
if (!requireNamespace("cowplot", quietly = TRUE)) stop("Package cowplot is required.")

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

NSAMPLE_ORDER <- c(3, 5, 10, 30, 50)
DATASET_ORDER <- c("E-ENAD-34", "GSE150910", "LUSC")

NORMAL_POINT_SIZE  <- 0.18
NORMAL_POINT_ALPHA <- 0.40
OUTLIER_POINT_SIZE <- 0.72
OUTLIER_POINT_ALPHA <- 0.95
POINT_JITTER_WIDTH <- 0.055

canonical_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x
}

parse_meta <- function(f) {
  stem <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(stem, "-", fixed = TRUE)[[1]]

  if (length(parts) < 2) {
    warning("Cannot parse filename: ", basename(f))
    return(NULL)
  }

  nSample <- suppressWarnings(as.integer(parts[1]))
  if (!is.finite(nSample)) {
    warning("Cannot parse nSample from filename: ", basename(f))
    return(NULL)
  }

  dataset <- paste(parts[-1], collapse = "-")
  list(nSample = nSample, Dataset = dataset)
}

tidy_auc_object <- function(x, file, meta) {
  handle_df <- function(df, replicate_default = NA_integer_) {
    cn <- colnames(df)
    method_col <- intersect(c("methods", "method", "Method", "METHOD"), cn)
    auc_col <- intersect(c("auc", "AUC", "Auc"), cn)
    rep_col <- intersect(c("replicate", "rep", "sim", "iteration", "iter"), cn)

    if (length(method_col) >= 1 && length(auc_col) >= 1) {
      reps <- if (length(rep_col) >= 1) df[[rep_col[1]]] else seq_len(nrow(df))
      return(data.frame(
        Dataset = meta$Dataset,
        nSample = meta$nSample,
        replicate = reps,
        Method = canonical_method(df[[method_col[1]]]),
        AUC = suppressWarnings(as.numeric(df[[auc_col[1]]])),
        source_file = file,
        stringsAsFactors = FALSE
      ))
    }

    numeric_cols <- cn[vapply(df, is.numeric, logical(1))]
    method_cols <- numeric_cols[canonical_method(numeric_cols) %in% METHOD_ORDER]

    if (length(method_cols) > 0) {
      out <- list()
      k <- 1L
      for (cc in method_cols) {
        out[[k]] <- data.frame(
          Dataset = meta$Dataset,
          nSample = meta$nSample,
          replicate = seq_len(nrow(df)),
          Method = canonical_method(cc),
          AUC = suppressWarnings(as.numeric(df[[cc]])),
          source_file = file,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
      return(do.call(rbind, out))
    }

    rn <- canonical_method(rownames(df))
    if (length(intersect(rn, METHOD_ORDER)) > 0 && ncol(df) >= 1) {
      return(data.frame(
        Dataset = meta$Dataset,
        nSample = meta$nSample,
        replicate = replicate_default,
        Method = rn,
        AUC = suppressWarnings(as.numeric(df[[1]])),
        source_file = file,
        stringsAsFactors = FALSE
      ))
    }

    data.frame()
  }

  if (is.data.frame(x)) return(handle_df(x))
  if (is.matrix(x)) return(handle_df(as.data.frame(x)))

  if (is.numeric(x) && !is.null(names(x))) {
    return(data.frame(
      Dataset = meta$Dataset,
      nSample = meta$nSample,
      replicate = 1L,
      Method = canonical_method(names(x)),
      AUC = suppressWarnings(as.numeric(x)),
      source_file = file,
      stringsAsFactors = FALSE
    ))
  }

  if (is.list(x)) {
    out <- list()
    k <- 1L

    for (i in seq_along(x)) {
      xi <- x[[i]]
      tmp <- data.frame()

      if (is.data.frame(xi)) {
        tmp <- handle_df(xi, replicate_default = i)
      } else if (is.matrix(xi)) {
        tmp <- handle_df(as.data.frame(xi), replicate_default = i)
      } else if (is.numeric(xi) && !is.null(names(xi))) {
        tmp <- data.frame(
          Dataset = meta$Dataset,
          nSample = meta$nSample,
          replicate = i,
          Method = canonical_method(names(xi)),
          AUC = suppressWarnings(as.numeric(xi)),
          source_file = file,
          stringsAsFactors = FALSE
        )
      }

      if (nrow(tmp) > 0) {
        if (!"replicate" %in% colnames(tmp) || all(is.na(tmp$replicate))) tmp$replicate <- i
        out[[k]] <- tmp
        k <- k + 1L
      }
    }

    if (length(out) > 0) return(do.call(rbind, out))
  }

  warning("Unsupported AUC object structure in: ", file)
  data.frame()
}

all_files <- list.files(AUC_DIR, full.names = TRUE, recursive = FALSE)
files <- all_files[tolower(tools::file_ext(all_files)) == "rds"]
files <- sort(files)

print(basename(files))

if (length(files) == 0) {
  stop("No RDS files found by file_ext() in AUC_DIR: ", AUC_DIR)
}

auc_long <- data.frame()
structure_log <- character()

for (f in files) {
  meta <- parse_meta(f)
  if (is.null(meta)) next

  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error")) {
    warning("Cannot read RDS: ", f, " | ", conditionMessage(obj))
    next
  }

  structure_log <- c(
    structure_log,
    paste0("\n\n================ ", basename(f), " ================"),
    capture.output(str(obj, max.level = 2))
  )

  tmp <- tidy_auc_object(obj, f, meta)
  if (nrow(tmp) > 0) auc_long <- rbind(auc_long, tmp)
}

if (nrow(auc_long) == 0) {
  stop("No AUC rows extracted. Check the input RDS structure.")
}

auc_long$Method <- canonical_method(auc_long$Method)
auc_long <- auc_long[auc_long$Method %in% METHOD_ORDER, , drop = FALSE]
auc_long$nSample <- as.integer(auc_long$nSample)
auc_long <- auc_long[auc_long$nSample %in% NSAMPLE_ORDER, , drop = FALSE]

auc_long$AUC <- suppressWarnings(as.numeric(auc_long$AUC))
auc_long <- auc_long[is.finite(auc_long$AUC) & !is.na(auc_long$AUC) & auc_long$AUC >= 0 & auc_long$AUC <= 1, , drop = FALSE]

if (nrow(auc_long) == 0) {
  stop("All extracted AUC rows were filtered out. Check AUC values and method names.")
}

auc_long$Method <- factor(auc_long$Method, levels = METHOD_ORDER)
auc_long$nSample <- factor(auc_long$nSample, levels = NSAMPLE_ORDER)

dataset_levels <- unique(c(
  DATASET_ORDER[DATASET_ORDER %in% unique(auc_long$Dataset)],
  sort(setdiff(unique(auc_long$Dataset), DATASET_ORDER))
))
auc_long$Dataset <- factor(auc_long$Dataset, levels = dataset_levels)

write.csv(auc_long, file.path(OUT_DIR, "Fig2_data.csv"), row.names = FALSE)

auc_summary <- aggregate(AUC ~ Dataset + nSample + Method, auc_long, function(x) {
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
auc_summary <- do.call(data.frame, auc_summary)
colnames(auc_summary) <- sub("AUC.", "", colnames(auc_summary), fixed = TRUE)
write.csv(auc_summary, file.path(OUT_DIR, "Fig2_summary.csv"), row.names = FALSE)

calc_outlier_flag <- function(df) {
  q1 <- unname(quantile(df$AUC, 0.25, na.rm = TRUE))
  q3 <- unname(quantile(df$AUC, 0.75, na.rm = TRUE))
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr

  df$q1 <- q1
  df$q3 <- q3
  df$iqr <- iqr
  df$lower_whisker_rule <- lower
  df$upper_whisker_rule <- upper
  df$is_outlier <- df$AUC < lower | df$AUC > upper
  df
}

split_key <- interaction(auc_long$Dataset, auc_long$nSample, auc_long$Method, drop = TRUE)
auc_with_outlier <- do.call(rbind, lapply(split(auc_long, split_key), calc_outlier_flag))

write.csv(auc_with_outlier, file.path(OUT_DIR, "Fig2_data_with_outlier_flag.csv"), row.names = FALSE)
write.csv(auc_with_outlier[auc_with_outlier$is_outlier, ],
          file.path(OUT_DIR, "Fig2_outliers.csv"),
          row.names = FALSE)

get_breaks <- function(v) {
  br <- pretty(range(v, na.rm = TRUE), n = 5)
  br[br >= 0 & br <= 1]
}

plot_format <- theme(
  axis.line = element_line(color = "black", linewidth = 0.25),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 5.2, color = "black"),
  plot.title = element_text(size = 7.2, hjust = 0.5, margin = margin(b = 1)),
  plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")
)

make_one_panel <- function(s, ds, row_i, col_j, n_rows) {
  dat <- auc_with_outlier[
    as.character(auc_with_outlier$nSample) == as.character(s) &
      as.character(auc_with_outlier$Dataset) == as.character(ds),
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) return(ggplot() + theme_void())

  y_breaks <- get_breaks(dat$AUC)
  if (length(y_breaks) < 2) y_breaks <- pretty(c(0, 1), n = 5)

  ylab_text <- if (col_j == 1) "AUC" else NULL
  subtitle_text <- if (col_j == 1) paste0("n = ", s) else NULL

  ggplot(dat, aes(x = Method, y = AUC, color = Method)) +
    geom_point(
      data = dat[!dat$is_outlier, , drop = FALSE],
      size = NORMAL_POINT_SIZE,
      alpha = NORMAL_POINT_ALPHA,
      position = position_jitter(width = POINT_JITTER_WIDTH, height = 0, seed = 8848)
    ) +
    geom_boxplot(
      outlier.shape = NA,
      linewidth = 0.28,
      width = 0.55,
      fill = NA
    ) +
    geom_point(
      data = dat[dat$is_outlier, , drop = FALSE],
      size = OUTLIER_POINT_SIZE,
      alpha = OUTLIER_POINT_ALPHA,
      position = position_jitter(width = POINT_JITTER_WIDTH, height = 0, seed = 8848)
    ) +
    scale_color_manual(values = METHOD_COLS, limits = METHOD_ORDER, drop = FALSE) +
    scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.03, 0.05))) +
    labs(x = NULL, y = ylab_text, title = as.character(ds), subtitle = subtitle_text) +
    plot_format +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(size = 6.8, face = "bold", hjust = 0, margin = margin(b = 1)),
      axis.text.x = if (row_i == n_rows) element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.7) else element_blank(),
      axis.ticks.x = if (row_i == n_rows) element_line(color = "black", linewidth = 0.25) else element_blank(),
      axis.title.y = if (col_j == 1) element_text(size = 7) else element_blank(),
      axis.text.y = element_text(size = 5.2)
    )
}

plotlist <- list()
labels <- character()

for (i in seq_along(NSAMPLE_ORDER)) {
  n_now <- NSAMPLE_ORDER[i]
  for (j in seq_along(dataset_levels)) {
    ds_now <- dataset_levels[j]
    plotlist[[length(plotlist) + 1]] <- make_one_panel(n_now, ds_now, i, j, length(NSAMPLE_ORDER))
    labels <- c(labels, ifelse(j == 1, LETTERS[i], ""))
  }
}

fig <- cowplot::plot_grid(
  plotlist = plotlist,
  ncol = length(dataset_levels),
  labels = labels,
  label_size = 11,
  align = "hv"
)

title_grob <- cowplot::ggdraw() +
  cowplot::draw_label("Negative-binomial-based standard simulation AUC results", fontface = "bold", size = 11, x = 0.02, hjust = 0)

fig_final <- cowplot::plot_grid(title_grob, fig, ncol = 1, rel_heights = c(0.045, 1))

pdf_file  <- file.path(OUT_DIR, "Fig2.pdf")
tiff_file <- file.path(OUT_DIR, "Fig2.tiff")
png_file  <- file.path(OUT_DIR, "Fig2.png")

ggsave(pdf_file,  fig_final, width = 6.4, height = 7.0, units = "in")
ggsave(tiff_file, fig_final, width = 6.4, height = 7.0, units = "in", dpi = 600, compression = "lzw")
ggsave(png_file,  fig_final, width = 6.4, height = 7.0, units = "in", dpi = 300)

ggsave(file.path(OUT_DIR, "Fig2_small.pdf"),
       fig_final, width = 5.8, height = 6.4, units = "in")
ggsave(file.path(OUT_DIR, "Fig2_small.tiff"),
       fig_final, width = 5.8, height = 6.4, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(OUT_DIR, "Fig2_small.png"),
       fig_final, width = 5.8, height = 6.4, units = "in", dpi = 300)

