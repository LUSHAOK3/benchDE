# Figure 14 source: overall-performance bubble panels.

set.seed(8848)

BASE_DIR <- Sys.getenv("MULTIDIMDE_BASE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
SAVE_TIFF <- toupper(Sys.getenv("SAVE_TIFF", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")

NB_STD_DATP_DIR   <- file.path(BASE_DIR, "sim1", "datp")
NB_STD_AUC_DIR    <- file.path(BASE_DIR, "sim1", "auc")
NB_OUT_AUC_DIR    <- file.path(BASE_DIR, "sim2", "auc")
NB_NOISE_AUC_DIR  <- file.path(BASE_DIR, "sim3_restored", "auc")

POI_STD_DATP_DIR  <- file.path(BASE_DIR, "poisson", "datp")
POI_STD_AUC_DIR   <- file.path(BASE_DIR, "poisson", "auc")
POI_OUT_AUC_DIR   <- file.path(BASE_DIR, "poisson_sim2", "auc")
POI_NOISE_AUC_DIR <- file.path(BASE_DIR, "poisson_sim3", "auc")

BRCA_FPR_DIR      <- file.path(BASE_DIR, "FPR", "datp")
BRCA_CAT_DIR      <- file.path(BASE_DIR, "CAT", "conc")
SCLC_FPR_DIR      <- file.path(BASE_DIR, "FPR_SCLC", "datp")
SCLC_CAT_DIR      <- file.path(BASE_DIR, "CAT_SCLC", "conc")

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(BASE_DIR, "figures", paste0("Fig14_", TIMESTAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUBBLE_LABEL_TYPE <- Sys.getenv("BUBBLE_LABEL_TYPE", unset = "score_raw")
if (!BUBBLE_LABEL_TYPE %in% c("score_raw", "raw_value", "score_plot", "none")) {
  stop("BUBBLE_LABEL_TYPE must be 'score_raw', 'raw_value', 'score_plot', or 'none'.")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

METHOD_ORDER <- c(
  "DSS", "ROTS", "T.test", "Wilcoxon", "ABSSeq", "voom",
  "NOISeq", "DESeq", "edgeR.lrt", "edgeR.qlf", "DESeq2", "NBPSeq"
)

NSAMPLE_ORDER <- c(3, 5, 10, 30, 50)
LOW_FDR_MAX <- 0.10
FPR_CUTOFF <- 0.05
CAT_TOP_PROP <- 0.10

NORMALIZE_SCOPE <- Sys.getenv("NORMALIZE_SCOPE", unset = "panel_nsample")
if (!NORMALIZE_SCOPE %in% c("panel_nsample", "panel", "none")) {
  stop("NORMALIZE_SCOPE must be 'panel_nsample', 'panel', or 'none'.")
}

BUBBLE_COLORS <- c("#2047c7", "#7b62aa", "#c8b07d", "#e0ad55", "#d86e35", "#b80000")

clean_method <- function(x) {
  x <- as.character(x)
  x[x == "T-test"] <- "T.test"
  x[x == "T.test."] <- "T.test"
  x
}

minmax01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & !is.na(x)
  if (!any(ok)) return(out)
  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(out)
  if (abs(rng[2] - rng[1]) < .Machine$double.eps) {
    out[ok] <- 1
  } else {
    out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
  }
  out
}

list_rds_files <- function(d) {
  if (!dir.exists(d)) return(character())
  fs <- list.files(d, full.names = TRUE, recursive = FALSE)
  sort(fs[tolower(tools::file_ext(fs)) == "rds"])
}

safe_read_rds <- function(f) {
  tryCatch(readRDS(f), error = function(e) {
    warning("Failed to read: ", f, " | ", conditionMessage(e))
    NULL
  })
}

parse_file_general <- function(file) {
  b <- basename(file)
  b <- sub("\\.rds$", "", b)
  parts <- strsplit(b, "-", fixed = TRUE)[[1]]

  if (length(parts) >= 3 && suppressWarnings(!is.na(as.numeric(parts[1])))) {
    level <- suppressWarnings(as.numeric(parts[1]))
    n <- suppressWarnings(as.integer(parts[2]))
    dataset <- paste(parts[-c(1, 2)], collapse = "-")
    return(data.frame(Level = level, nSample = n, Dataset = dataset, stringsAsFactors = FALSE))
  }

  if (length(parts) >= 2) {
    n <- suppressWarnings(as.integer(parts[1]))
    dataset <- paste(parts[-1], collapse = "-")
    return(data.frame(Level = NA_real_, nSample = n, Dataset = dataset, stringsAsFactors = FALSE))
  }

  data.frame(Level = NA_real_, nSample = NA_integer_, Dataset = NA_character_, stringsAsFactors = FALSE)
}

read_auc_files <- function(auc_dir, source_name = "auc") {
  files <- list_rds_files(auc_dir)
  if (length(files) == 0) {
    warning("No AUC RDS files found in: ", auc_dir)
    return(NULL)
  }

  dat_list <- lapply(files, function(f) {
    x <- safe_read_rds(f)
    if (is.null(x)) return(NULL)

    if (is.list(x) && !is.data.frame(x)) {
      ok <- vapply(x, function(z) is.data.frame(z) || is.matrix(z), logical(1))
      if (any(ok)) {
        x <- do.call(rbind, lapply(x[ok], as.data.frame))
      } else if (!is.null(names(x)) && any(clean_method(names(x)) %in% METHOD_ORDER)) {
        x <- data.frame(methods = clean_method(names(x)), auc = as.numeric(x), stringsAsFactors = FALSE)
      }
    }

    if (is.matrix(x)) x <- as.data.frame(x)
    if (!is.data.frame(x)) return(NULL)

    x$file <- basename(f)
    x
  })

  dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
  if (length(dat_list) == 0) return(NULL)

  dat <- do.call(rbind, dat_list)
  dat$file_no_ext <- sub("\\.rds$", "", dat$file)

  if (!"methods" %in% colnames(dat)) {
    possible <- intersect(c("method", "Method", "methods", "algorithm", "Algorithm"), colnames(dat))
    if (length(possible) > 0) dat$methods <- dat[[possible[1]]]
  }

  if (!"auc" %in% colnames(dat)) {
    possible <- intersect(c("AUC", "auc", "AUROC", "auroc", "value", "Value"), colnames(dat))
    if (length(possible) > 0) dat$auc <- dat[[possible[1]]]
  }

  if (!all(c("methods", "auc") %in% colnames(dat))) {
    warning("AUC files parsed but required columns methods/auc were not found: ", auc_dir)
    return(NULL)
  }

  dat$methods <- clean_method(dat$methods)

  meta <- do.call(rbind, lapply(dat$file_no_ext, parse_file_general))
  dat$Level <- meta$Level
  dat$nSample <- meta$nSample
  dat$Dataset <- meta$Dataset
  dat$source <- source_name
  dat
}

make_auc_metric <- function(auc_dir, group_name, metric_name, panel_label) {
  dat <- read_auc_files(auc_dir, source_name = paste(group_name, metric_name))
  if (is.null(dat)) return(NULL)

  dat$auc <- suppressWarnings(as.numeric(dat$auc))
  dat <- dat[dat$methods %in% METHOD_ORDER & dat$nSample %in% NSAMPLE_ORDER &
               is.finite(dat$auc) & !is.na(dat$auc), , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  raw <- aggregate(auc ~ methods + nSample, dat, mean, na.rm = TRUE)
  data.frame(
    Group = group_name,
    Metric = metric_name,
    ScoreMetric = metric_name,
    Panel = panel_label,
    methods = raw$methods,
    nSample = raw$nSample,
    raw_value = raw$auc,
    oriented_value = raw$auc,
    direction = "higher is better",
    stringsAsFactors = FALSE
  )
}

make_robustness_metric <- function(auc_dir, group_name, metric_name, panel_label) {
  dat <- read_auc_files(auc_dir, source_name = paste(group_name, metric_name))
  if (is.null(dat)) return(NULL)

  dat$auc <- suppressWarnings(as.numeric(dat$auc))
  dat$Level <- suppressWarnings(as.numeric(dat$Level))
  dat <- dat[dat$methods %in% METHOD_ORDER & dat$nSample %in% NSAMPLE_ORDER &
               is.finite(dat$auc) & is.finite(dat$Level), , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  keys <- interaction(dat$methods, dat$nSample, dat$Dataset, drop = TRUE)
  slope_list <- lapply(split(dat, keys), function(d) {
    if (nrow(d) < 2 || length(unique(d$Level)) < 2) return(NULL)
    fit <- lm(auc ~ Level, data = d)
    data.frame(
      methods = d$methods[1],
      nSample = d$nSample[1],
      Dataset = d$Dataset[1],
      slope = unname(coef(fit)[2]),
      stringsAsFactors = FALSE
    )
  })
  slope_df <- do.call(rbind, slope_list)
  if (is.null(slope_df) || nrow(slope_df) == 0) return(NULL)

  raw <- aggregate(slope ~ methods + nSample, slope_df, mean, na.rm = TRUE)

  data.frame(
    Group = group_name,
    Metric = metric_name,
    ScoreMetric = metric_name,
    Panel = panel_label,
    methods = raw$methods,
    nSample = raw$nSample,
    raw_value = raw$slope,
    oriented_value = -abs(raw$slope),
    direction = "slope closer to zero is better: oriented_value = -abs(slope)",
    stringsAsFactors = FALSE
  )
}

extract_p_matrix <- function(rep_obj) {
  if (inherits(rep_obj, "try-error")) return(NULL)

  if (is.list(rep_obj) && !is.null(rep_obj$p)) {
    pmat <- as.data.frame(rep_obj$p, check.names = FALSE)
    colnames(pmat) <- clean_method(colnames(pmat))
    return(pmat)
  }

  if (is.data.frame(rep_obj)) {
    label_col <- intersect(c("label", "truth"), colnames(rep_obj))
    pcols <- setdiff(colnames(rep_obj), label_col)
    pmat <- as.data.frame(rep_obj[, pcols, drop = FALSE], check.names = FALSE)
    colnames(pmat) <- clean_method(colnames(pmat))
    return(pmat)
  }

  NULL
}

extract_label <- function(rep_obj) {
  if (inherits(rep_obj, "try-error")) return(NULL)

  if (is.list(rep_obj)) {
    lab <- rep_obj$label
    if (is.null(lab)) lab <- rep_obj$truth
    if (!is.null(lab)) return(as.numeric(lab))
  }

  if (is.data.frame(rep_obj)) {
    if ("label" %in% colnames(rep_obj)) return(as.numeric(rep_obj$label))
    if ("truth" %in% colnames(rep_obj)) return(as.numeric(rep_obj$truth))
  }

  NULL
}

label_to_binary <- function(label) {
  if (is.logical(label)) return(as.integer(label))
  if (is.numeric(label) || is.integer(label)) return(as.integer(label == 1))
  y <- toupper(trimws(as.character(label)))
  as.integer(y %in% c("1", "TRUE", "T", "DE", "DEG", "DIFF", "DIFFERENTIAL", "UP", "DOWN"))
}

calc_low_fdr_tpr_repeat <- function(pmat, label, low_fdr_max = LOW_FDR_MAX) {
  if (is.null(pmat) || is.null(label)) return(NULL)

  pmat <- as.data.frame(pmat, check.names = FALSE)
  colnames(pmat) <- clean_method(colnames(pmat))

  if (!is.null(names(label)) && !is.null(rownames(pmat))) {
    common <- intersect(names(label), rownames(pmat))
    if (length(common) > 10) {
      pmat <- pmat[common, , drop = FALSE]
      label <- label[common]
    }
  } else {
    L <- min(length(label), nrow(pmat))
    pmat <- pmat[seq_len(L), , drop = FALSE]
    label <- label[seq_len(L)]
  }

  label <- label_to_binary(label)
  n_pos <- sum(label == 1, na.rm = TRUE)
  n_neg <- sum(label == 0, na.rm = TRUE)
  if (n_pos == 0 || n_neg == 0) return(NULL)

  methods <- intersect(METHOD_ORDER, colnames(pmat))
  if (length(methods) == 0) return(NULL)

  fdr_grid <- seq(0, low_fdr_max, by = 0.01)

  out <- lapply(methods, function(m) {
    pv <- suppressWarnings(as.numeric(pmat[[m]]))
    ok <- is.finite(pv) & !is.na(pv) & !is.na(label)
    pv <- pv[ok]
    lab <- label[ok]

    n_pos_m <- sum(lab == 1, na.rm = TRUE)
    n_neg_m <- sum(lab == 0, na.rm = TRUE)
    if (length(pv) == 0 || n_pos_m == 0 || n_neg_m == 0) return(NULL)

    ord <- order(pv, decreasing = FALSE, na.last = NA)
    lab_ord <- lab[ord]

    tp <- cumsum(lab_ord == 1)
    fp <- cumsum(lab_ord == 0)
    called <- seq_along(tp)

    empirical_fdr <- fp / pmax(called, 1)
    tpr <- tp / n_pos_m

    tpr_values <- vapply(fdr_grid, function(g) {
      idx <- which(empirical_fdr <= g)
      if (length(idx) == 0) return(0)
      max(tpr[idx], na.rm = TRUE)
    }, numeric(1))

    data.frame(
      methods = m,
      low_fdr_tpr = mean(tpr_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

process_fdr_tpr_file <- function(f) {
  meta <- parse_file_general(basename(f))
  dat <- safe_read_rds(f)
  if (is.null(dat) || !is.list(dat)) return(NULL)

  reps <- lapply(seq_along(dat), function(i) {
    pmat <- extract_p_matrix(dat[[i]])
    lab <- extract_label(dat[[i]])
    tmp <- calc_low_fdr_tpr_repeat(pmat, lab, LOW_FDR_MAX)
    if (is.null(tmp)) return(NULL)
    tmp$replicate <- i
    tmp$nSample <- meta$nSample
    tmp$Dataset <- meta$Dataset
    tmp
  })
  reps <- reps[!vapply(reps, is.null, logical(1))]
  if (length(reps) == 0) return(NULL)
  do.call(rbind, reps)
}

make_fdr_tpr_metric <- function(datp_dir, group_name, metric_name, panel_label) {
  files <- list_rds_files(datp_dir)
  if (length(files) == 0) return(NULL)

  all <- lapply(files, process_fdr_tpr_file)
  all <- all[!vapply(all, is.null, logical(1))]
  if (length(all) == 0) return(NULL)

  dat <- do.call(rbind, all)
  dat <- dat[dat$methods %in% METHOD_ORDER & dat$nSample %in% NSAMPLE_ORDER, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  raw <- aggregate(low_fdr_tpr ~ methods + nSample, dat, mean, na.rm = TRUE)

  data.frame(
    Group = group_name,
    Metric = metric_name,
    ScoreMetric = metric_name,
    Panel = panel_label,
    methods = raw$methods,
    nSample = raw$nSample,
    raw_value = raw$low_fdr_tpr,
    oriented_value = raw$low_fdr_tpr,
    direction = paste0("higher mean max TPR under empirical FDR <= ", LOW_FDR_MAX, " is better"),
    stringsAsFactors = FALSE
  )
}

calc_fpr_repeat <- function(pmat, cutoff = FPR_CUTOFF) {
  if (is.null(pmat)) return(NULL)
  pmat <- as.data.frame(pmat, check.names = FALSE)
  colnames(pmat) <- clean_method(colnames(pmat))
  methods <- intersect(METHOD_ORDER, colnames(pmat))
  if (length(methods) == 0) return(NULL)

  out <- lapply(methods, function(m) {
    pv <- suppressWarnings(as.numeric(pmat[[m]]))
    pv <- pv[is.finite(pv) & !is.na(pv) & pv >= 0 & pv <= 1]
    if (length(pv) == 0) return(NULL)
    qv <- p.adjust(pv, method = "BH")
    data.frame(methods = m, fpr = mean(qv <= cutoff), stringsAsFactors = FALSE)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

process_fpr_file <- function(f) {
  meta <- parse_file_general(basename(f))
  dat <- safe_read_rds(f)
  if (is.null(dat)) return(NULL)

  if (is.data.frame(dat)) dat <- list(dat)
  if (!is.list(dat)) return(NULL)

  reps <- lapply(seq_along(dat), function(i) {
    pmat <- extract_p_matrix(dat[[i]])
    if (is.null(pmat) && is.data.frame(dat[[i]])) {
      pmat <- as.data.frame(dat[[i]], check.names = FALSE)
      colnames(pmat) <- clean_method(colnames(pmat))
    }
    tmp <- calc_fpr_repeat(pmat, FPR_CUTOFF)
    if (is.null(tmp)) return(NULL)
    tmp$replicate <- i
    tmp$nSample <- meta$nSample
    tmp$Dataset <- meta$Dataset
    tmp
  })

  reps <- reps[!vapply(reps, is.null, logical(1))]
  if (length(reps) == 0) return(NULL)
  do.call(rbind, reps)
}

make_fpr_metric <- function(fpr_dir, group_name, metric_name, panel_label) {
  files <- list_rds_files(fpr_dir)
  if (length(files) == 0) return(NULL)

  all <- lapply(files, process_fpr_file)
  all <- all[!vapply(all, is.null, logical(1))]
  if (length(all) == 0) return(NULL)

  dat <- do.call(rbind, all)
  dat <- dat[dat$methods %in% METHOD_ORDER & dat$nSample %in% NSAMPLE_ORDER, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  raw <- aggregate(fpr ~ methods + nSample, dat, mean, na.rm = TRUE)

  data.frame(
    Group = group_name,
    Metric = metric_name,
    ScoreMetric = "FPR-Control",
    Panel = panel_label,
    methods = raw$methods,
    nSample = raw$nSample,
    raw_value = raw$fpr,
    oriented_value = -raw$fpr,
    direction = paste0("lower FPR under BH-FDR <= ", FPR_CUTOFF, " is better: oriented_value = -FPR"),
    stringsAsFactors = FALSE
  )
}

choose_conc_at_top <- function(df, method_name, top_prop = CAT_TOP_PROP) {
  if (is.null(df) || nrow(df) == 0) return(NA_real_)

  required <- c("method1", "method2", "rank", "conc")
  if (!all(required %in% colnames(df))) return(NA_real_)

  d <- df
  d$method1 <- clean_method(d$method1)
  d$method2 <- clean_method(d$method2)
  d <- d[d$method1 == method_name & d$method2 == method_name, , drop = FALSE]
  if (nrow(d) == 0) return(NA_real_)

  d$rank_num <- suppressWarnings(as.numeric(d$rank))
  d$conc_num <- suppressWarnings(as.numeric(d$conc))
  d <- d[is.finite(d$rank_num) & is.finite(d$conc_num), , drop = FALSE]
  if (nrow(d) == 0) return(NA_real_)

  if ("nfeatures" %in% colnames(d)) {
    nf <- suppressWarnings(as.numeric(d$nfeatures))
    nf <- nf[is.finite(nf)]
    if (length(nf) > 0) {
      target <- max(1, median(nf, na.rm = TRUE) * top_prop)
    } else {
      target <- max(d$rank_num, na.rm = TRUE) * top_prop
    }
  } else {
    target <- max(d$rank_num, na.rm = TRUE) * top_prop
  }

  idx <- which.min(abs(d$rank_num - target))
  d$conc_num[idx]
}

process_cat_file <- function(f) {
  meta <- parse_file_general(basename(f))
  dat <- safe_read_rds(f)
  if (is.null(dat)) return(NULL)
  if (is.data.frame(dat)) dat <- list(dat)
  if (!is.list(dat)) return(NULL)

  reps <- lapply(seq_along(dat), function(i) {
    df <- dat[[i]]
    if (!is.data.frame(df)) return(NULL)

    vals <- lapply(METHOD_ORDER, function(m) {
      data.frame(
        methods = m,
        stability = choose_conc_at_top(df, m, CAT_TOP_PROP),
        replicate = i,
        nSample = meta$nSample,
        Dataset = meta$Dataset,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, vals)
  })

  reps <- reps[!vapply(reps, is.null, logical(1))]
  if (length(reps) == 0) return(NULL)
  do.call(rbind, reps)
}

make_stability_metric <- function(cat_dir, group_name, metric_name, panel_label) {
  files <- list_rds_files(cat_dir)
  if (length(files) == 0) return(NULL)

  all <- lapply(files, process_cat_file)
  all <- all[!vapply(all, is.null, logical(1))]
  if (length(all) == 0) return(NULL)

  dat <- do.call(rbind, all)
  dat <- dat[dat$methods %in% METHOD_ORDER & dat$nSample %in% NSAMPLE_ORDER &
               is.finite(dat$stability) & !is.na(dat$stability), , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  raw <- aggregate(stability ~ methods + nSample, dat, mean, na.rm = TRUE)

  data.frame(
    Group = group_name,
    Metric = metric_name,
    ScoreMetric = "Stability",
    Panel = panel_label,
    methods = raw$methods,
    nSample = raw$nSample,
    raw_value = raw$stability,
    oriented_value = raw$stability,
    direction = paste0("higher CAT concordance at top ", CAT_TOP_PROP * 100, "% is better"),
    stringsAsFactors = FALSE
  )
}

metric_list <- list(
  make_auc_metric(NB_STD_AUC_DIR, "NB simulation", "AUC", "NB\nAUC"),
  make_fdr_tpr_metric(NB_STD_DATP_DIR, "NB simulation", "FDR-TPR", "NB\nFDR-TPR"),
  make_robustness_metric(NB_OUT_AUC_DIR, "NB simulation", "O-Robustness", "NB\nO-Robustness"),
  make_robustness_metric(NB_NOISE_AUC_DIR, "NB simulation", "N-Robustness", "NB\nN-Robustness"),

  make_auc_metric(POI_STD_AUC_DIR, "Poisson simulation", "AUC", "Poisson\nAUC"),
  make_fdr_tpr_metric(POI_STD_DATP_DIR, "Poisson simulation", "FDR-TPR", "Poisson\nFDR-TPR"),
  make_robustness_metric(POI_OUT_AUC_DIR, "Poisson simulation", "O-Robustness", "Poisson\nO-Robustness"),
  make_robustness_metric(POI_NOISE_AUC_DIR, "Poisson simulation", "N-Robustness", "Poisson\nN-Robustness"),

  make_fpr_metric(BRCA_FPR_DIR, "Real data", "BRCA FPR", "BRCA\nFPR"),
  make_stability_metric(BRCA_CAT_DIR, "Real data", "BRCA Stability", "BRCA\nStability"),
  make_fpr_metric(SCLC_FPR_DIR, "Real data", "SCLC FPR", "SCLC\nFPR"),
  make_stability_metric(SCLC_CAT_DIR, "Real data", "SCLC Stability", "SCLC\nStability")
)

metric_list <- metric_list[!vapply(metric_list, is.null, logical(1))]
if (length(metric_list) == 0) stop("No metric data were extracted.")

raw_all <- do.call(rbind, metric_list)
raw_all$methods <- clean_method(raw_all$methods)
raw_all <- raw_all[
  raw_all$methods %in% METHOD_ORDER &
    raw_all$nSample %in% NSAMPLE_ORDER &
    is.finite(raw_all$oriented_value) &
    !is.na(raw_all$oriented_value),
  ,
  drop = FALSE
]

panel_order <- c(
  "NB\nAUC", "NB\nFDR-TPR", "NB\nO-Robustness", "NB\nN-Robustness",
  "Poisson\nAUC", "Poisson\nFDR-TPR", "Poisson\nO-Robustness", "Poisson\nN-Robustness",
  "BRCA\nFPR", "BRCA\nStability", "SCLC\nFPR", "SCLC\nStability"
)

scoremetric_order <- c("AUC", "FDR-TPR", "O-Robustness", "N-Robustness", "FPR-Control", "Stability")

raw_all$Panel <- factor(raw_all$Panel, levels = panel_order)
raw_all$ScoreMetric <- factor(raw_all$ScoreMetric, levels = scoremetric_order)
raw_all$methods <- factor(raw_all$methods, levels = rev(METHOD_ORDER))
raw_all$nSample <- factor(raw_all$nSample, levels = NSAMPLE_ORDER)

score_all <- raw_all
score_all$score <- NA_real_

for (sm in levels(score_all$ScoreMetric)) {
  idx <- which(score_all$ScoreMetric == sm)
  if (length(idx) > 0) {
    score_all$score[idx] <- minmax01(score_all$oriented_value[idx])
  }
}

score_all <- score_all[is.finite(score_all$score) & !is.na(score_all$score), , drop = FALSE]

score_all$score_raw <- score_all$score

if (NORMALIZE_SCOPE == "panel_nsample") {
  grp <- interaction(score_all$Panel, score_all$nSample, drop = TRUE)
  score_all$score_plot <- ave(score_all$score_raw, grp, FUN = minmax01)
} else if (NORMALIZE_SCOPE == "panel") {
  score_all$score_plot <- ave(score_all$score_raw, score_all$Panel, FUN = minmax01)
} else {
  score_all$score_plot <- score_all$score_raw
}

score_all$raw_display <- score_all$raw_value
robust_idx <- as.character(score_all$ScoreMetric) %in% c("O-Robustness", "N-Robustness")
score_all$raw_display[robust_idx] <- abs(score_all$raw_value[robust_idx])

fmt2 <- function(z) {
  ifelse(is.na(z), "", sprintf("%.2f", z))
}
fmt3 <- function(z) {
  ifelse(is.na(z), "", sprintf("%.3f", z))
}

score_all$label_score_raw <- fmt2(score_all$score_raw)
score_all$label_score_plot <- fmt2(score_all$score_plot)
score_all$label_raw_value <- ifelse(
  as.character(score_all$ScoreMetric) %in% c("O-Robustness", "N-Robustness"),
  fmt3(score_all$raw_display),
  fmt2(score_all$raw_display)
)

score_all$bubble_label <- switch(
  BUBBLE_LABEL_TYPE,
  "score_raw" = score_all$label_score_raw,
  "score_plot" = score_all$label_score_plot,
  "raw_value" = score_all$label_raw_value,
  "none" = ""
)

score_all$label_color <- ifelse(score_all$score_plot >= 0.68, "white", "black")

write.csv(raw_all, file.path(OUT_DIR, "Fig14_data.csv"), row.names = FALSE)
write.csv(score_all, file.path(OUT_DIR, "Fig14_scores.csv"), row.names = FALSE)

panel_status <- aggregate(score_plot ~ Panel, score_all, function(x) c(n = length(x), min = min(x), max = max(x)))
panel_status <- do.call(data.frame, panel_status)
write.csv(panel_status, file.path(OUT_DIR, "Fig14_panel_ranges.csv"), row.names = FALSE)

scoremetric_status <- aggregate(score_plot ~ ScoreMetric, score_all, function(x) c(n = length(x), min = min(x), max = max(x)))
scoremetric_status <- do.call(data.frame, scoremetric_status)
write.csv(scoremetric_status, file.path(OUT_DIR, "Fig14_metric_ranges.csv"), row.names = FALSE)

score_all$bubble_label <- ""
score_all$label_color <- "black"

pdf_file <- file.path(OUT_DIR, "Fig14.pdf")

make_main_panel_plot <- function(df_panel, panel_name) {
  ggplot(df_panel, aes(x = nSample, y = methods)) +
    geom_point(
      aes(size = score_plot, fill = score_plot),
      shape = 21, color = "grey35", stroke = 0.25, alpha = 0.98,
      show.legend = FALSE
    ) +
    scale_size_continuous(range = c(1.8, 7.2), limits = c(0, 1)) +
    scale_fill_gradientn(colours = BUBBLE_COLORS, limits = c(0, 1)) +
    labs(title = panel_name, x = NULL, y = NULL) +
    theme_bw(base_size = 8.8) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.text.x = element_text(size = 7.2, angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 7.8),
      axis.title = element_blank(),
      plot.title = element_text(size = 10.2, face = "bold", hjust = 0.5, margin = margin(b = 3)),
      plot.margin = margin(4, 2, 4, 4)
    )
}

make_scale_panel_plot <- function(panel_name, df_panel) {
  raw_min <- suppressWarnings(min(df_panel$score_raw, na.rm = TRUE))
  raw_max <- suppressWarnings(max(df_panel$score_raw, na.rm = TRUE))
  fmt2 <- function(z) ifelse(is.finite(z), sprintf("%.2f", z), "NA")

  bar_df <- data.frame(
    x = 1,
    y = seq(0, 1, length.out = 101)
  )
  bar_df$score_plot <- bar_df$y

  tick_df <- data.frame(
    y = c(1, 0.5, 0),
    label = c("1.00", "0.50", "0.00")
  )

  ggplot(bar_df, aes(x = x, y = y)) +
    annotate(
      "rect",
      xmin = 0.62, xmax = 1.66,
      ymin = -0.12, ymax = 1.20,
      fill = "white", color = "grey55", linewidth = 0.35
    ) +
    geom_tile(
      aes(fill = score_plot),
      width = 0.26,
      height = 0.018,
      show.legend = FALSE
    ) +
    geom_segment(
      data = tick_df,
      aes(x = 1.15, xend = 1.22, y = y, yend = y),
      inherit.aes = FALSE,
      color = "grey35",
      linewidth = 0.25
    ) +
    geom_text(
      data = tick_df,
      aes(x = 1.25, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.15,
      color = "black"
    ) +
    annotate(
      "text",
      x = 1.14, y = 1.12,
      label = "Score",
      fontface = "bold",
      size = 2.45
    ) +
    annotate(
      "text",
      x = 1.14, y = -0.22,
      label = paste0("raw: ", fmt2(raw_min), "–", fmt2(raw_max)),
      size = 1.90,
      color = "grey20"
    ) +
    scale_fill_gradientn(colours = BUBBLE_COLORS, limits = c(0, 1)) +
    coord_cartesian(xlim = c(0.55, 1.68), ylim = c(-0.28, 1.25), clip = "off") +
    theme_void(base_size = 8.5) +
    theme(
      plot.margin = margin(12, 8, 10, 0)
    )
}

panel_levels <- levels(score_all$Panel)
if (is.null(panel_levels) || length(panel_levels) == 0) {
  panel_levels <- unique(as.character(score_all$Panel))
}

main_plots <- list()
scale_plots <- list()
for (pn in panel_levels) {
  dfp <- subset(score_all, Panel == pn)
  main_plots[[pn]] <- make_main_panel_plot(dfp, pn)
  scale_plots[[pn]] <- make_scale_panel_plot(pn, dfp)
}

grDevices::pdf(pdf_file, width = 20, height = 11.2, onefile = TRUE, paper = "special")
grid::grid.newpage()
lay <- grid::grid.layout(
  nrow = 3, ncol = 8,
  widths = grid::unit(rep(c(6.25, 1.25), 4), "null"),
  heights = grid::unit(rep(1, 3), "null")
)
grid::pushViewport(grid::viewport(layout = lay))

for (i in seq_along(panel_levels)) {
  pn <- panel_levels[i]
  rr <- ((i - 1) %/% 4) + 1
  cc_main  <- (((i - 1) %% 4) * 2) + 1
  cc_scale <- cc_main + 1

  print(main_plots[[pn]],
        vp = grid::viewport(layout.pos.row = rr, layout.pos.col = cc_main))
  print(scale_plots[[pn]],
        vp = grid::viewport(layout.pos.row = rr, layout.pos.col = cc_scale))
}

grid::upViewport()
grDevices::dev.off()

