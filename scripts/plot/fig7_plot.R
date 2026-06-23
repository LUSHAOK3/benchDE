# Figure 7: Poisson outlier robustness.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
AUC_DIR <- Sys.getenv("POISSON_SIM2_AUC_DIR", unset = "")

if (!nzchar(AUC_DIR)) {
  candidate_auc_dirs <- c(
    file.path(BASE_DIR, "poisson_sim2", "auc")
  )
  candidate_auc_dirs <- candidate_auc_dirs[dir.exists(candidate_auc_dirs)]
  if (length(candidate_auc_dirs) == 0) {
    stop(
      "No Poisson outlier AUC directory found. Please set POISSON_SIM2_AUC_DIR manually."
    )
  }
  AUC_DIR <- candidate_auc_dirs[1]
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig7_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(AUC_DIR)) stop("AUC_DIR does not exist: ", AUC_DIR)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

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

canonical_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x
}

tidy_auc_object <- function(x, file, meta) {
  handle_df <- function(df, rep_default = NA_integer_) {
    cn <- colnames(df)
    method_col <- intersect(c("methods", "method", "Method", "METHOD"), cn)
    auc_col <- intersect(c("auc", "AUC", "Auc"), cn)
    rep_col <- intersect(c("replicate", "rep", "sim", "iteration", "iter"), cn)

    if (length(method_col) >= 1 && length(auc_col) >= 1) {
      reps <- if (length(rep_col) >= 1) df[[rep_col[1]]] else seq_len(nrow(df))
      return(data.frame(
        replicate = reps,
        Method = canonical_method(df[[method_col[1]]]),
        AUC = suppressWarnings(as.numeric(df[[auc_col[1]]])),
        source_file = basename(file),
        stringsAsFactors = FALSE
      ))
    }

    numeric_cols <- cn[vapply(df, is.numeric, logical(1))]
    method_cols <- numeric_cols[canonical_method(numeric_cols) %in% METHOD_ORDER]
    if (length(method_cols) > 0) {
      out <- lapply(method_cols, function(cc) {
        data.frame(
          replicate = seq_len(nrow(df)),
          Method = canonical_method(cc),
          AUC = suppressWarnings(as.numeric(df[[cc]])),
          source_file = basename(file),
          stringsAsFactors = FALSE
        )
      })
      return(do.call(rbind, out))
    }

    rn <- canonical_method(rownames(df))
    if (length(intersect(rn, METHOD_ORDER)) > 0 && ncol(df) >= 1) {
      return(data.frame(
        replicate = rep_default,
        Method = rn,
        AUC = suppressWarnings(as.numeric(df[[1]])),
        source_file = basename(file),
        stringsAsFactors = FALSE
      ))
    }

    data.frame()
  }

  if (is.data.frame(x)) return(handle_df(x))
  if (is.matrix(x)) return(handle_df(as.data.frame(x)))

  if (is.numeric(x) && !is.null(names(x))) {
    return(data.frame(
      replicate = 1L,
      Method = canonical_method(names(x)),
      AUC = suppressWarnings(as.numeric(x)),
      source_file = basename(file),
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
        tmp <- handle_df(xi, rep_default = i)
      } else if (is.matrix(xi)) {
        tmp <- handle_df(as.data.frame(xi), rep_default = i)
      } else if (is.numeric(xi) && !is.null(names(xi))) {
        tmp <- data.frame(
          replicate = i,
          Method = canonical_method(names(xi)),
          AUC = suppressWarnings(as.numeric(xi)),
          source_file = basename(file),
          stringsAsFactors = FALSE
        )
      }
      if (nrow(tmp) > 0) {
        out[[k]] <- tmp
        k <- k + 1L
      }
    }
    if (length(out) > 0) return(do.call(rbind, out))
  }

  data.frame()
}

DATASET_ORDER <- c("E-ENAD-34", "GSE150910", "LUSC")
ROW_LETTERS <- LETTERS[seq_along(NSAMPLE_ORDER)]

parse_meta <- function(f) {
  stem <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(stem, "-", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  prop <- suppressWarnings(as.numeric(parts[1]))
  nSample <- suppressWarnings(as.integer(parts[2]))
  dataset <- paste(parts[-c(1, 2)], collapse = "-")
  if (!is.finite(prop) || !is.finite(nSample)) return(NULL)
  list(prop = prop, nSample = nSample, Dataset = dataset)
}

files_all <- list.files(AUC_DIR, full.names = TRUE, recursive = FALSE)
files <- sort(files_all[tolower(tools::file_ext(files_all)) == "rds"])
if (length(files) == 0) stop("No .rds files found in AUC_DIR: ", AUC_DIR)

auc_long <- data.frame()
structure_log <- character()

for (f in files) {
  meta <- parse_meta(f)
  if (is.null(meta)) next
  if (!meta$nSample %in% NSAMPLE_ORDER) next
  if (!meta$Dataset %in% DATASET_ORDER) next

  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error")) next

  structure_log <- c(
    structure_log,
    paste0("\n\n================ ", basename(f), " ================"),
    capture.output(str(obj, max.level = 2))
  )

  tmp <- tidy_auc_object(obj, f, meta)
  if (nrow(tmp) > 0) {
    tmp$prop <- meta$prop
    tmp$Dataset <- meta$Dataset
    tmp$nSample <- meta$nSample
    auc_long <- rbind(auc_long, tmp)
  }
}

auc_long <- auc_long[
  auc_long$Method %in% METHOD_ORDER &
    auc_long$nSample %in% NSAMPLE_ORDER &
    auc_long$Dataset %in% DATASET_ORDER &
    is.finite(auc_long$prop) &
    is.finite(auc_long$AUC) &
    !is.na(auc_long$AUC),
  ,
  drop = FALSE
]
if (nrow(auc_long) == 0) stop("No usable AUC rows extracted.")

auc_long$Method <- factor(auc_long$Method, levels = METHOD_ORDER)
auc_long$Dataset <- factor(auc_long$Dataset, levels = DATASET_ORDER)
write.csv(auc_long, file.path(OUT_DIR, "Fig7_data.csv"), row.names = FALSE)

auc_summary <- aggregate(AUC ~ prop + Dataset + nSample + Method, auc_long, function(x) {
  c(n = length(x), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
})
auc_summary <- do.call(data.frame, auc_summary)
colnames(auc_summary) <- sub("AUC.", "", colnames(auc_summary), fixed = TRUE)
auc_summary$Method <- factor(auc_summary$Method, levels = METHOD_ORDER)
auc_summary$Dataset <- factor(auc_summary$Dataset, levels = DATASET_ORDER)
write.csv(auc_summary, file.path(OUT_DIR, "Fig7_summary.csv"), row.names = FALSE)

split_key <- interaction(auc_summary$nSample, auc_summary$Dataset, auc_summary$Method, drop = TRUE)
slope_list <- lapply(split(auc_summary, split_key), function(df) {
  if (nrow(df) < 2) return(NULL)
  fit <- lm(mean ~ prop, data = df)
  data.frame(
    nSample = df$nSample[1],
    Dataset = as.character(df$Dataset[1]),
    Method = as.character(df$Method[1]),
    Slope = unname(coef(fit)[["prop"]]),
    stringsAsFactors = FALSE
  )
})
slope_df <- do.call(rbind, slope_list)
if (is.null(slope_df) || nrow(slope_df) == 0) stop("Failed to compute slopes.")

slope_df$Method <- factor(slope_df$Method, levels = METHOD_ORDER)
slope_df$Dataset <- factor(slope_df$Dataset, levels = DATASET_ORDER)
slope_df$nSample <- as.integer(slope_df$nSample)
write.csv(slope_df, file.path(OUT_DIR, "Fig7_slope.csv"), row.names = FALSE)

panel_theme <- theme_bw(base_size = 7.2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line = element_blank(),
    axis.text = element_text(color = "black", size = 5.2),
    axis.title = element_text(size = 7.0),
    plot.title = element_text(size = 6.2, hjust = 0.5, margin = margin(b = 0.35)),
    plot.margin = margin(t = 0.05, r = 2.0, b = 0.05, l = 2.0),
    legend.position = "none",
    aspect.ratio = 0.52
  )

make_panel <- function(n_now, dataset_now, row_i, col_i) {
  dat <- slope_df[
    slope_df$nSample == n_now &
      as.character(slope_df$Dataset) == dataset_now,
    ,
    drop = FALSE
  ]
  if (nrow(dat) == 0) return(ggplot() + theme_void())

  y_rng <- range(dat$Slope, na.rm = TRUE)
  pad <- max(0.0008, diff(y_rng) * 0.12)
  y_lim <- c(y_rng[1] - pad, y_rng[2] + pad)

  ggplot(dat, aes(x = Method, y = Slope, color = Method)) +
    geom_errorbar(aes(ymin = Slope, ymax = Slope), width = 0.54, linewidth = 0.9, show.legend = FALSE) +
    scale_color_manual(values = METHOD_COLS, limits = METHOD_ORDER, drop = FALSE) +
    scale_y_continuous(limits = y_lim, breaks = pretty(y_lim, n = 4), expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      title = dataset_now,
      x = NULL,
      y = if (col_i == 1) "Slope" else NULL
    ) +
    panel_theme +
    theme(
      axis.title.y = if (col_i == 1) element_text(size = 7.0) else element_blank(),
      axis.text.x = if (row_i == length(NSAMPLE_ORDER)) {
        element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.0)
      } else {
        element_blank()
      },
      axis.ticks.x = element_line(linewidth = 0.25)
    )
}

make_label_panel <- function(letter, n_now) {
  ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    annotate("text", x = 0.03, y = 0.90, label = letter, fontface = "bold", size = 6.4, hjust = 0) +
    annotate("text", x = 0.45, y = 0.90, label = paste0("n = ", n_now), size = 3.3, hjust = 0) +
    theme_void() +
    theme(plot.margin = margin(t = -1.5, r = 0.3, b = -1.5, l = 5))
}

panel_list <- vector("list", length(NSAMPLE_ORDER) * length(DATASET_ORDER))
panel_index <- 1
for (i in seq_along(NSAMPLE_ORDER)) {
  for (j in seq_along(DATASET_ORDER)) {
    panel_list[[panel_index]] <- make_panel(NSAMPLE_ORDER[i], DATASET_ORDER[j], i, j)
    panel_index <- panel_index + 1
  }
}

aligned_panels <- cowplot::align_plots(plotlist = panel_list, align = "hv", axis = "tblr")

row_plots <- vector("list", length(NSAMPLE_ORDER))
panel_index <- 1
for (i in seq_along(NSAMPLE_ORDER)) {
  lab <- make_label_panel(ROW_LETTERS[i], NSAMPLE_ORDER[i])
  p1 <- aligned_panels[[panel_index]]; panel_index <- panel_index + 1
  p2 <- aligned_panels[[panel_index]]; panel_index <- panel_index + 1
  p3 <- aligned_panels[[panel_index]]; panel_index <- panel_index + 1

  row_plots[[i]] <- plot_grid(
    lab, p1, p2, p3,
    ncol = 4,
    rel_widths = c(0.38, 1, 1, 1),
    align = "hv",
    axis = "tblr"
  )
}

fig <- plot_grid(
  plotlist = row_plots,
  ncol = 1,
  rel_heights = rep(1, length(row_plots)),
  align = "v",
  axis = "lr"
)

pdf_file <- file.path(OUT_DIR, "Fig7.pdf")

grDevices::pdf(
  file = pdf_file,
  width = 6.68,
  height = 7.20,
  onefile = TRUE,
  useDingbats = FALSE,
  paper = "special",
  family = "Helvetica"
)
print(fig)
grDevices::dev.off()

SAVE_TIFF <- toupper(Sys.getenv("SAVE_TIFF", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
if (SAVE_TIFF) {
  if (!requireNamespace("ragg", quietly = TRUE)) {
    warning("SAVE_TIFF=TRUE but package 'ragg' is not installed. TIFF was not generated.")
  } else {
    tiff_file <- file.path(OUT_DIR, "Fig7.tiff")
    ragg::agg_tiff(
      filename = tiff_file,
      width = 6.68,
      height = 7.20,
      units = "in",
      res = 600,
      compression = "lzw",
      background = "white",
      scaling = 1
    )
    print(fig)
    grDevices::dev.off()
  }
}

