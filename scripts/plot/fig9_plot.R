# Figure 9: Poisson noise robustness.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
AUC_DIR <- Sys.getenv(
  "POISSON_SIM3_AUC_DIR",
  unset = file.path(BASE_DIR, "poisson_sim3", "auc")
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig9_", timestamp))
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

DATASET_KEEP <- "E-ENAD-34"
NOISE_LEVELS_FOR_SLOPE <- c(0.2, 0.4, 0.6, 0.8)
ROW_LETTERS <- LETTERS[seq_along(NSAMPLE_ORDER)]

parse_meta <- function(f) {
  stem <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(stem, "-", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  noise <- suppressWarnings(as.numeric(parts[1]))
  nSample <- suppressWarnings(as.integer(parts[2]))
  dataset <- paste(parts[-c(1, 2)], collapse = "-")
  if (!is.finite(noise) || !is.finite(nSample)) return(NULL)
  list(noise = noise, nSample = nSample, Dataset = dataset)
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
  if (meta$Dataset != DATASET_KEEP) next

  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error")) next

  structure_log <- c(
    structure_log,
    paste0("\n\n================ ", basename(f), " ================"),
    capture.output(str(obj, max.level = 2))
  )

  tmp <- tidy_auc_object(obj, f, meta)
  if (nrow(tmp) > 0) {
    tmp$noise <- meta$noise
    tmp$Dataset <- meta$Dataset
    tmp$nSample <- meta$nSample
    auc_long <- rbind(auc_long, tmp)
  }
}

auc_long <- auc_long[
  auc_long$Method %in% METHOD_ORDER &
    auc_long$nSample %in% NSAMPLE_ORDER &
    auc_long$Dataset == DATASET_KEEP &
    is.finite(auc_long$noise) &
    is.finite(auc_long$AUC) &
    !is.na(auc_long$AUC),
  ,
  drop = FALSE
]
if (nrow(auc_long) == 0) stop("No usable AUC rows extracted.")

auc_long$Method <- factor(auc_long$Method, levels = METHOD_ORDER)
write.csv(auc_long, file.path(OUT_DIR, "Fig9_data.csv"), row.names = FALSE)

auc_summary <- aggregate(AUC ~ noise + Dataset + nSample + Method, auc_long, function(x) {
  c(n = length(x), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
})
auc_summary <- do.call(data.frame, auc_summary)
colnames(auc_summary) <- sub("AUC.", "", colnames(auc_summary), fixed = TRUE)
auc_summary$Method <- factor(auc_summary$Method, levels = METHOD_ORDER)
write.csv(auc_summary, file.path(OUT_DIR, "Fig9_summary.csv"), row.names = FALSE)

slope_input <- auc_summary[auc_summary$noise %in% NOISE_LEVELS_FOR_SLOPE, , drop = FALSE]
if (nrow(slope_input) == 0) {
  stop("No rows for NOISE_LEVELS_FOR_SLOPE = ", paste(NOISE_LEVELS_FOR_SLOPE, collapse = ", "))
}

split_key <- interaction(slope_input$nSample, slope_input$Method, drop = TRUE)
slope_list <- lapply(split(slope_input, split_key), function(df) {
  if (nrow(df) < 2) return(NULL)
  fit <- lm(mean ~ noise, data = df)
  data.frame(
    nSample = df$nSample[1],
    Dataset = DATASET_KEEP,
    Method = as.character(df$Method[1]),
    Slope = unname(coef(fit)[["noise"]]),
    stringsAsFactors = FALSE
  )
})
slope_df <- do.call(rbind, slope_list)
if (is.null(slope_df) || nrow(slope_df) == 0) stop("Failed to compute slopes.")

slope_df$Method <- factor(slope_df$Method, levels = METHOD_ORDER)
slope_df$nSample <- as.integer(slope_df$nSample)
write.csv(slope_df, file.path(OUT_DIR, "Fig9_slope.csv"), row.names = FALSE)

panel_theme <- theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line = element_blank(),
    axis.text = element_text(color = "black", size = 5.8),
    axis.title = element_text(size = 8),
    plot.margin = margin(t = 2, r = 6, b = 6, l = 7),
    legend.position = "none"
  )

make_panel <- function(n_now) {
  dat <- slope_df[slope_df$nSample == n_now, , drop = FALSE]
  if (nrow(dat) == 0) return(ggplot() + theme_void())

  y_rng <- range(dat$Slope, na.rm = TRUE)
  pad <- max(0.004, diff(y_rng) * 0.12)
  y_lim <- c(y_rng[1] - pad, y_rng[2] + pad)

  ggplot(dat, aes(x = Method, y = Slope, color = Method)) +
    geom_errorbar(aes(ymin = Slope, ymax = Slope), width = 0.55, linewidth = 1.0, show.legend = FALSE) +
    scale_color_manual(values = METHOD_COLS, limits = METHOD_ORDER, drop = FALSE) +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0.02, 0.03))) +
    labs(x = NULL, y = "Slope") +
    panel_theme +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.6),
      axis.ticks.x = element_line(linewidth = 0.25)
    )
}

make_header <- function(letter, n_now) {
  ggdraw() +
    draw_label(letter, x = 0.000, y = 0.54, hjust = 0, vjust = 0.5, fontface = "bold", size = 16) +
    draw_label(paste0("n = ", n_now), x = 0.095, y = 0.54, hjust = 0, vjust = 0.5, size = 8.2) +
    draw_label(DATASET_KEEP, x = 0.54, y = 0.54, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 8)
}

make_cell <- function(letter, n_now) {
  header <- make_header(letter, n_now)
  panel <- make_panel(n_now)
  plot_grid(header, panel, ncol = 1, rel_heights = c(0.16, 1))
}

pA <- make_cell("A", 3)
pB <- make_cell("B", 5)
pC <- make_cell("C", 10)
pD <- make_cell("D", 30)
pE <- make_cell("E", 50)
blank <- ggplot() + theme_void()

row1 <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1), align = "hv", axis = "tblr")
row2 <- plot_grid(pC, pD, ncol = 2, rel_widths = c(1, 1), align = "hv", axis = "tblr")
row3 <- plot_grid(pE, blank, ncol = 2, rel_widths = c(1, 1), align = "hv", axis = "tblr")

fig <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 1, 1))

ggsave(file.path(OUT_DIR, "Fig9.pdf"),
       fig, width = 9.4, height = 8.4, units = "in", limitsize = FALSE)
ggsave(file.path(OUT_DIR, "Fig9.tiff"),
       fig, width = 9.4, height = 8.4, units = "in", dpi = 600, compression = "lzw", limitsize = FALSE)
ggsave(file.path(OUT_DIR, "Fig9.png"),
       fig, width = 9.4, height = 8.4, units = "in", dpi = 300, limitsize = FALSE)

