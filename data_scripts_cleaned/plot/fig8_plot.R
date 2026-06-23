# Figure 8: negative-binomial noise robustness.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = "/public/workspace/fangzhaoyuan/lsk/plos/MultiDimDE")
AUC_DIR <- Sys.getenv("NB_SIM3_AUC_DIR", unset = file.path(BASE_DIR, "sim3_restored", "auc"))
OUT_DIR <- file.path(BASE_DIR, "Fig8")
DATASET_TO_PLOT <- Sys.getenv("NB_SIM3_DATASET", unset = "E-ENAD-34")

NOISE_LEVELS_MAIN <- as.numeric(trimws(strsplit(
  Sys.getenv("NB_SIM3_NOISE_LEVELS", unset = "0.2,0.4,0.6,0.8"),
  ",",
  fixed = TRUE
)[[1]]))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

method_order <- c(
  "ABSSeq", "DESeq", "DESeq2", "DSS",
  "edgeR.lrt", "edgeR.qlf", "NBPSeq", "NOISeq",
  "ROTS", "T.test", "voom", "Wilcoxon"
)

sample_order <- c(3, 5, 10, 30, 50)

method_cols <- c(
  "ABSSeq"    = "#6B6B6B",
  "DESeq"     = "#6A3D9A",
  "DESeq2"    = "#CAB2D6",
  "DSS"       = "#FDBF6F",
  "edgeR.lrt" = "#A6CEE3",
  "edgeR.qlf" = "#56B4E9",
  "NBPSeq"    = "#009E73",
  "NOISeq"    = "#CC79A7",
  "ROTS"      = "#E31A1C",
  "T.test"    = "#FB9A99",
  "voom"      = "#1F355E",
  "Wilcoxon"  = "#D95F02"
)


clean_method <- function(x) {
  x <- as.character(x)
  x[x %in% c("T-test", "Ttest", "ttest", "t.test", "T test")] <- "T.test"
  x[grepl("^wilcox", x, ignore.case = TRUE)] <- "Wilcoxon"
  x
}

safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) {
    warning("Failed to readRDS: ", path, " | ", conditionMessage(e))
    NULL
  })
}

parse_file_general <- function(file_base) {
  x <- sub("\\.rds$", "", basename(file_base))
  parts <- strsplit(x, "-", fixed = TRUE)[[1]]
  nums <- suppressWarnings(as.numeric(parts))

  level <- NA_real_
  n_sample <- NA_real_
  dataset <- NA_character_

  if (length(parts) >= 3 && is.finite(nums[1]) && is.finite(nums[2])) {
    level <- nums[1]
    n_sample <- nums[2]
    dataset <- paste(parts[-c(1, 2)], collapse = "-")
  } else {
    idx <- which(is.finite(nums) & nums %in% sample_order)
    if (length(idx) > 0) {
      n_sample <- nums[idx[1]]
      if (idx[1] > 1 && is.finite(nums[idx[1] - 1])) level <- nums[idx[1] - 1]
      if (idx[1] < length(parts)) dataset <- paste(parts[(idx[1] + 1):length(parts)], collapse = "-")
    }
  }

  data.frame(
    file = basename(file_base),
    Level = level,
    nSample = n_sample,
    Dataset = dataset,
    stringsAsFactors = FALSE
  )
}

read_auc_files <- function(auc_dir) {
  files <- list.files(auc_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No auc .rds files found in: ", auc_dir)

  dat_list <- lapply(files, function(f) {
    x <- safe_read_rds(f)
    if (is.null(x)) return(NULL)

    if (is.list(x) && !is.data.frame(x)) {
      x <- do.call(rbind, x)
    }
    if (!is.data.frame(x)) {
      x <- as.data.frame(x)
    }

    x$file <- basename(f)
    x
  })

  dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
  if (length(dat_list) == 0) stop("No readable auc files in: ", auc_dir)

  dat <- do.call(rbind, dat_list)

  if (!"methods" %in% colnames(dat)) {
    possible <- intersect(c("methods", "method", "Method", "algorithm", "Algorithm"), colnames(dat))
    if (length(possible) > 0) dat$methods <- dat[[possible[1]]]
  }

  if (!"auc" %in% colnames(dat)) {
    possible <- intersect(c("auc", "AUC", "value", "Value", "score", "Score"), colnames(dat))
    if (length(possible) > 0) dat$auc <- dat[[possible[1]]]
  }

  if (!all(c("methods", "auc", "file") %in% colnames(dat))) {
    stop("AUC file structure not recognized. Required columns: methods and auc.")
  }

  dat$methods <- clean_method(dat$methods)
  dat$auc <- suppressWarnings(as.numeric(dat$auc))

  meta <- do.call(rbind, lapply(dat$file, parse_file_general))
  dat$Level <- meta$Level
  dat$nSample <- meta$nSample
  dat$Dataset <- meta$Dataset

  dat
}

make_auc_summary <- function(auc_dat) {
  out <- aggregate(
    auc ~ Dataset + nSample + Level + methods,
    data = auc_dat,
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  colnames(out)[colnames(out) == "auc"] <- "mean_auc"
  out
}

calc_slope_table <- function(auc_summary, levels_to_use, label) {
  dat <- auc_summary[
    auc_summary$Level %in% levels_to_use &
      auc_summary$Dataset == DATASET_TO_PLOT &
      auc_summary$nSample %in% sample_order &
      auc_summary$methods %in% method_order,
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) {
    stop("No AUC summary rows left for slope set: ", label)
  }

  groups <- unique(dat[, c("Dataset", "nSample", "methods")])
  slope_list <- vector("list", nrow(groups))

  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    sub <- dat[
      dat$Dataset == g$Dataset &
        dat$nSample == g$nSample &
        dat$methods == g$methods,
      ,
      drop = FALSE
    ]

    if (length(unique(sub$Level)) < 2) {
      slope <- NA_real_
      intercept <- NA_real_
    } else {
      fit <- lm(mean_auc ~ Level, data = sub)
      slope <- unname(coef(fit)[2])
      intercept <- unname(coef(fit)[1])
    }

    slope_list[[i]] <- data.frame(
      Dataset = g$Dataset,
      nSample = as.numeric(g$nSample),
      methods = clean_method(g$methods),
      slope = slope,
      intercept = intercept,
      n_levels = length(unique(sub$Level)),
      levels_used = paste(sort(unique(sub$Level)), collapse = ","),
      mean_auc = mean(sub$mean_auc, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  slope_df <- do.call(rbind, slope_list)

  slope_df <- slope_df[
    slope_df$methods %in% method_order &
      slope_df$nSample %in% sample_order &
      is.finite(slope_df$slope),
    ,
    drop = FALSE
  ]

  slope_df$methods <- factor(slope_df$methods, levels = method_order)
  slope_df$nSample <- factor(slope_df$nSample, levels = sample_order)
  slope_df$calculation <- label

  slope_df
}

nice_ylim <- function(y) {
  y <- y[is.finite(y)]
  if (length(y) == 0) return(c(-1, 1))

  yr <- range(y)
  pad <- max(0.010, diff(yr) * 0.18)
  if (!is.finite(pad) || pad == 0) {
    pad <- max(0.010, abs(yr[1]) * 0.15)
  }
  c(yr[1] - pad, yr[2] + pad)
}

plot_one_panel <- function(df_panel, panel_letter, n_value, dataset_name) {
  ylim <- nice_ylim(df_panel$slope)

  plot(
    NA,
    xlim = c(0.5, length(method_order) + 0.5),
    ylim = ylim,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "Slope",
    main = dataset_name,
    bty = "o",
    cex.main = 0.86,
    cex.lab = 0.82,
    mgp = c(1.55, 0.38, 0),
    xaxs = "i",
    yaxs = "i"
  )

  axis(2, las = 1, cex.axis = 0.72, tck = -0.025, mgp = c(1.35, 0.35, 0))
  axis(
    1,
    at = seq_along(method_order),
    labels = method_order,
    las = 2,
    cex.axis = 0.66,
    tck = -0.025,
    mgp = c(1.0, 0.35, 0)
  )

  for (m in method_order) {
    row <- df_panel[as.character(df_panel$methods) == m, , drop = FALSE]
    if (nrow(row) == 0) next

    x <- match(m, method_order)
    yy <- row$slope[1]

    segments(
      x - 0.26,
      yy,
      x + 0.26,
      yy,
      col = method_cols[m],
      lwd = 2.4,
      lend = "round"
    )
  }

  mtext(panel_letter, side = 3, line = 1.25, adj = -0.18, cex = 1.45, font = 2)
  mtext(paste0("n = ", n_value), side = 3, line = 1.18, adj = 0.28, cex = 0.72)
}

draw_fig <- function(slope_df, file_prefix, width = 8.4, height = 7.2, res = 300) {
  draw_device <- function(device, file) {
    if (device == "pdf") {
      grDevices::pdf(file, width = width, height = height)
    } else if (device == "png") {
      grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
      grDevices::tiff(file, width = width, height = height, units = "in", res = res, compression = "lzw")
    }

    old_par <- par(no.readonly = TRUE)
    on.exit({
      par(old_par)
      grDevices::dev.off()
    }, add = TRUE)

    layout(
      matrix(c(1, 2,
               3, 4,
               5, 0), nrow = 3, byrow = TRUE),
      widths = c(1, 1),
      heights = c(1, 1, 1)
    )

    par(oma = c(0.15, 0.25, 0.35, 0.15), xpd = NA)

    for (i in seq_along(sample_order)) {
      n <- sample_order[i]
      par(mar = c(4.45, 3.35, 2.15, 1.05))

      df_panel <- slope_df[
        as.numeric(as.character(slope_df$nSample)) == n,
        ,
        drop = FALSE
      ]

      plot_one_panel(
        df_panel = df_panel,
        panel_letter = LETTERS[i],
        n_value = n,
        dataset_name = DATASET_TO_PLOT
      )
    }
  }

  pdf_file <- file.path(OUT_DIR, paste0(file_prefix, ".pdf"))
  png_file <- file.path(OUT_DIR, paste0(file_prefix, ".png"))
  tiff_file <- file.path(OUT_DIR, paste0(file_prefix, ".tiff"))

  draw_device("pdf", pdf_file)

  tryCatch(
    draw_device("png", png_file),
    error = function(e) warning("PNG export failed: ", conditionMessage(e))
  )

  tryCatch(
    draw_device("tiff", tiff_file),
    error = function(e) warning("TIFF export failed: ", conditionMessage(e))
  )

}

auc_dat <- read_auc_files(AUC_DIR)

auc_dat <- auc_dat[
  auc_dat$methods %in% method_order &
    auc_dat$nSample %in% sample_order &
    is.finite(auc_dat$Level) &
    is.finite(auc_dat$auc),
  ,
  drop = FALSE
]

if (nrow(auc_dat) == 0) {
  stop("AUC table became empty after filtering. Check method names and file names.")
}

write.csv(auc_dat, file.path(OUT_DIR, "Fig8_data.csv"), row.names = FALSE)

level_table <- as.data.frame.matrix(table(auc_dat$Level, auc_dat$nSample))
level_table$noise_level <- rownames(level_table)
level_table <- level_table[, c("noise_level", setdiff(colnames(level_table), "noise_level"))]
write.csv(level_table, file.path(OUT_DIR, "Fig8_levels.csv"), row.names = FALSE)

auc_summary <- make_auc_summary(auc_dat)
write.csv(auc_summary, file.path(OUT_DIR, "Fig8_summary.csv"), row.names = FALSE)

all_levels_for_dataset <- sort(unique(auc_summary$Level[auc_summary$Dataset == DATASET_TO_PLOT]))

slope_exclude0 <- calc_slope_table(
  auc_summary = auc_summary,
  levels_to_use = NOISE_LEVELS_MAIN,
  label = "exclude0_0.2_0.4_0.6_0.8"
)

slope_with0 <- calc_slope_table(
  auc_summary = auc_summary,
  levels_to_use = all_levels_for_dataset,
  label = "with0_all_levels"
)

write.csv(slope_exclude0, file.path(OUT_DIR, "Fig8_slope.csv"), row.names = FALSE)
write.csv(slope_with0, file.path(OUT_DIR, "Fig8_slope_with_zero.csv"), row.names = FALSE)

merge_key <- c("Dataset", "nSample", "methods")
cmp <- merge(
  slope_exclude0[, c(merge_key, "slope", "levels_used")],
  slope_with0[, c(merge_key, "slope", "levels_used")],
  by = merge_key,
  suffixes = c("_exclude0", "_with0")
)
cmp$delta_exclude0_minus_with0 <- cmp$slope_exclude0 - cmp$slope_with0

write.csv(
  cmp,
  file.path(OUT_DIR, "Fig8_slope_comparison.csv"),
  row.names = FALSE
)

draw_fig(
  slope_df = slope_exclude0,
  file_prefix = "Fig8"
)

draw_fig(
  slope_df = slope_with0,
  file_prefix = "Fig8_with_zero"
)
