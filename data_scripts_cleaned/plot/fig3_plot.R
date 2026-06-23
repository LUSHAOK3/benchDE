# Figure 3: negative-binomial standard-simulation FDR-TPR curves.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
DATP_DIR <- Sys.getenv("NB_STANDARD_DATP_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE/sim1/datp")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig3_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(DATP_DIR)) stop("DATP_DIR does not exist: ", DATP_DIR)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
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
ROW_LABELS <- paste0(LETTERS[seq_along(NSAMPLE_ORDER)], "   n = ", NSAMPLE_ORDER)
DATASET_ORDER <- c("E-ENAD-34", "GSE150910", "LUSC")
FDR_GRID <- seq(0, 0.15, by = 0.01)
FDR_DASH <- c(0.01, 0.05, 0.10)

canonical_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x
}

label_to_binary <- function(label) {
  if (is.logical(label)) return(as.integer(label))
  if (is.numeric(label) || is.integer(label)) return(as.integer(label == 1))
  y <- toupper(trimws(as.character(label)))
  as.integer(y %in% c("1", "TRUE", "T", "DE", "DEG", "DIFF", "DIFFERENTIAL", "UP", "DOWN"))
}

parse_meta <- function(f) {
  stem <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(stem, "-", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NULL)
  nSample <- suppressWarnings(as.integer(parts[1]))
  dataset <- paste(parts[-1], collapse = "-")
  if (!is.finite(nSample)) return(NULL)
  list(nSample = nSample, Dataset = dataset)
}

extract_p_label <- function(obj) {
  if (inherits(obj, "try-error")) stop("Replicate is try-error.")

  if (is.list(obj) && !is.null(obj$p)) {
    p <- as.data.frame(obj$p, check.names = FALSE)
    label <- obj$label
    if (is.null(label)) label <- obj$truth
    if (is.null(label)) stop("Missing label/truth in list(p=...) object.")
    return(list(p = p, label = label))
  }

  if (is.data.frame(obj) && ("label" %in% colnames(obj) || "truth" %in% colnames(obj))) {
    label_col <- if ("label" %in% colnames(obj)) "label" else "truth"
    label <- obj[[label_col]]
    p <- obj[, setdiff(colnames(obj), label_col), drop = FALSE]
    return(list(p = p, label = label))
  }

  stop("Unsupported datp replicate structure.")
}

align_label_p <- function(pmat, label) {
  pmat <- as.data.frame(pmat, check.names = FALSE)

  if (!is.null(names(label)) && !is.null(rownames(pmat))) {
    common <- intersect(names(label), rownames(pmat))
    if (length(common) > 10) {
      return(list(pmat = pmat[common, , drop = FALSE], label = label[common]))
    }
  }

  L <- min(length(label), nrow(pmat))
  list(pmat = pmat[seq_len(L), , drop = FALSE], label = label[seq_len(L)])
}

calc_curve_one_method <- function(pvec, label, method, rep_id, meta, fdr_grid = FDR_GRID) {
  pvec <- suppressWarnings(as.numeric(pvec))
  label <- label_to_binary(label)

  ok <- is.finite(pvec) & !is.na(pvec) & !is.na(label)
  pvec <- pvec[ok]
  label <- label[ok]

  n_pos <- sum(label == 1)
  n_neg <- sum(label == 0)
  if (n_pos == 0 || n_neg == 0 || length(pvec) == 0) return(data.frame())

  ord <- order(pvec, decreasing = FALSE, na.last = NA)
  lab_ord <- label[ord]

  tp <- cumsum(lab_ord == 1)
  fp <- cumsum(lab_ord == 0)
  called <- seq_along(tp)

  empirical_fdr <- fp / pmax(called, 1)
  tpr <- tp / n_pos

  tpr_grid <- vapply(fdr_grid, function(g) {
    idx <- which(empirical_fdr <= g)
    if (length(idx) == 0) return(0)
    max(tpr[idx], na.rm = TRUE)
  }, numeric(1))

  data.frame(
    Dataset = meta$Dataset,
    nSample = meta$nSample,
    FDR = fdr_grid,
    TPR = tpr_grid,
    Method = method,
    replicate = rep_id,
    stringsAsFactors = FALSE
  )
}

read_all_curves <- function(datp_dir) {
  files_all <- list.files(datp_dir, full.names = TRUE, recursive = FALSE)
  files <- sort(files_all[tolower(tools::file_ext(files_all)) == "rds"])
  if (length(files) == 0) stop("No RDS files found in DATP_DIR: ", datp_dir)

  all_curve <- list()
  skipped <- list()
  structure_log <- character()
  k <- 1L
  s <- 1L

  for (f in files) {
    meta <- parse_meta(f)
    if (is.null(meta)) next
    if (!meta$nSample %in% NSAMPLE_ORDER) next
    if (!meta$Dataset %in% DATASET_ORDER) next

    datp <- tryCatch(readRDS(f), error = function(e) e)
    if (inherits(datp, "error")) {
      skipped[[s]] <- data.frame(file = basename(f), replicate = NA_integer_, reason = conditionMessage(datp))
      s <- s + 1L
      next
    }

    structure_log <- c(
      structure_log,
      paste0("\n\n================ ", basename(f), " ================"),
      capture.output(str(datp, max.level = 2))
    )

    for (i in seq_along(datp)) {
      ep <- tryCatch(extract_p_label(datp[[i]]), error = function(e) e)
      if (inherits(ep, "error")) {
        skipped[[s]] <- data.frame(file = basename(f), replicate = i, reason = conditionMessage(ep))
        s <- s + 1L
        next
      }

      aligned <- align_label_p(ep$p, ep$label)
      pmat <- aligned$pmat
      label <- aligned$label
      colnames(pmat) <- canonical_method(colnames(pmat))
      methods <- intersect(METHOD_ORDER, colnames(pmat))

      for (m in methods) {
        tmp <- calc_curve_one_method(pmat[[m]], label, m, i, meta)
        if (nrow(tmp) > 0) {
          all_curve[[k]] <- tmp
          k <- k + 1L
        }
      }
    }
  }

  if (length(skipped) > 0) {
    write.csv(do.call(rbind, skipped), file.path(OUT_DIR, "Fig3_skipped.csv"), row.names = FALSE)
  }

  if (length(all_curve) == 0) stop("No curves computed from datp files.")
  do.call(rbind, all_curve)
}

curve_long <- read_all_curves(DATP_DIR)

curve_long$Method <- factor(curve_long$Method, levels = METHOD_ORDER)
curve_long$Dataset <- factor(curve_long$Dataset, levels = DATASET_ORDER)
curve_long$nSample <- as.integer(curve_long$nSample)
curve_long$rowLabel <- factor(
  paste0(LETTERS[match(curve_long$nSample, NSAMPLE_ORDER)], "   n = ", curve_long$nSample),
  levels = ROW_LABELS
)

write.csv(curve_long, file.path(OUT_DIR, "Fig3_data.csv"), row.names = FALSE)

curve_summary <- aggregate(TPR ~ Dataset + nSample + rowLabel + FDR + Method, curve_long, function(x) {
  c(n = length(x), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
})
curve_summary <- do.call(data.frame, curve_summary)
colnames(curve_summary) <- sub("TPR.", "", colnames(curve_summary), fixed = TRUE)
curve_summary$Method <- factor(curve_summary$Method, levels = METHOD_ORDER)
curve_summary$Dataset <- factor(curve_summary$Dataset, levels = DATASET_ORDER)
curve_summary$rowLabel <- factor(curve_summary$rowLabel, levels = ROW_LABELS)

write.csv(curve_summary, file.path(OUT_DIR, "Fig3_summary.csv"), row.names = FALSE)

p <- ggplot(curve_summary, aes(x = FDR, y = mean, color = Method, group = Method)) +
  geom_vline(xintercept = FDR_DASH, linetype = "dashed", linewidth = 0.30, color = "black") +
  geom_line(linewidth = 0.45) +
  facet_grid(rows = vars(rowLabel), cols = vars(Dataset), scales = "free_y", switch = "y") +
  scale_color_manual(values = METHOD_COLS, limits = METHOD_ORDER, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0, 0.05, 0.10, 0.15),
    labels = c("0.00", "0.05", "0.10", "0.15"),
    limits = c(-0.010, 0.15),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "Negative-binomial-based standard simulation FDR-TPR results",
    x = "FDR",
    y = "TPR",
    color = "Methods"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1.0))) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),

    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 7, face = "bold", margin = margin(b = 3)),
    strip.text.y.left = element_text(size = 7, face = "bold", angle = 0, margin = margin(r = 6)),

    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 8),

    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.box = "horizontal",
    legend.margin = margin(t = 2, b = 2),

    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(t = 6, r = 8, b = 6, l = 6),

    panel.spacing.x = unit(0.55, "lines"),
    panel.spacing.y = unit(0.70, "lines")
  )

ggsave(file.path(OUT_DIR, "Fig3.pdf"), p, width = 8.2, height = 8.8, units = "in", limitsize = FALSE)
ggsave(file.path(OUT_DIR, "Fig3.tiff"), p, width = 8.2, height = 8.8, units = "in", dpi = 600, compression = "lzw", limitsize = FALSE)
ggsave(file.path(OUT_DIR, "Fig3.png"), p, width = 8.2, height = 8.8, units = "in", dpi = 300, limitsize = FALSE)

