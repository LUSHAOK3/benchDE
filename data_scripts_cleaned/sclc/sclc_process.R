# run sclc all methods.

set.seed(8848)

BASE_DIR <- Sys.getenv("BENCHDE_DIR", unset = getwd())

SEED_VALUE <- 8848
WORKERS <- as.integer(Sys.getenv("SCLC_ALL_WORKERS", unset = "10"))
BATCH_SIZE <- as.integer(Sys.getenv("SCLC_ALL_BATCH_SIZE", unset = "10"))
NSAMPLE <- as.integer(strsplit(Sys.getenv("SCLC_ALL_NSAMPLE", unset = "3,5,10,30,50"), ",")[[1]])
RUN_MODULES <- toupper(gsub(" ", "", Sys.getenv("SCLC_ALL_MODULES", unset = "FPR,CAT")))
WRITE_ORIGINAL_DATP <- toupper(Sys.getenv("SCLC_WRITE_ORIGINAL_DATP", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
COMPUTE_CAT_CONC <- toupper(Sys.getenv("SCLC_COMPUTE_CAT_CONC", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")

Sys.setenv(
  BENCHDE_DIR = BASE_DIR,
  BENCHDE_SEED = as.character(SEED_VALUE),
  BENCHDE_CORES = as.character(WORKERS),
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1",
  SCLC_COUNTS = file.path(BASE_DIR, "datasets", "count_ensg_STAR.txt.gz"),
  SCLC_CLINICAL = file.path(BASE_DIR, "datasets", "Cell2023_SCLC_data_clinical.txt"),
  SCLC_NSAMPLE = paste(NSAMPLE, collapse = ","),
  SCLC_FPR_B = "100",
  SCLC_CAT_B = "10",
  SCLC_FPR_MODE = "normal_only",
  SCLC_RUN_FPR = "FALSE",
  SCLC_RUN_CAT = "FALSE",
  SCLC_OUT_SUFFIX = ""
)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::register(BiocParallel::SerialParam())
if (requireNamespace("dplyr", quietly = TRUE)) suppressPackageStartupMessages(library(dplyr))

UTILS_DIR <- file.path(BASE_DIR, "utils")
MAIN_SCRIPT <- file.path(UTILS_DIR, "sclc_preprocess.R")
RUNDE_FILE <- file.path(UTILS_DIR, "runDE.R")
CONC_FILE <- file.path(UTILS_DIR, "get_conc_function.R")

required_files <- c(MAIN_SCRIPT, RUNDE_FILE, CONC_FILE, Sys.getenv("SCLC_COUNTS"), Sys.getenv("SCLC_CLINICAL"))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing required files:\n", paste(missing_files, collapse = "\n"))

source(MAIN_SCRIPT)
source(RUNDE_FILE)

if (!exists("dge")) stop("Object 'dge' was not created by main script.")
if (!exists("algos")) stop("Object 'algos' was not created by runDE.R.")
if (!exists("B_FPR")) B_FPR <- as.integer(Sys.getenv("SCLC_FPR_B", unset = "100"))
if (!exists("B_CAT")) B_CAT <- as.integer(Sys.getenv("SCLC_CAT_B", unset = "10"))
if (!exists("FPR_MODE")) FPR_MODE <- Sys.getenv("SCLC_FPR_MODE", unset = "normal_only")

canonical_method <- function(x) ifelse(x == "T-test", "T.test", x)
sanitize_method <- function(x) gsub("[^A-Za-z0-9_.-]", "_", canonical_method(x))

requested_raw <- Sys.getenv("SCLC_ALL_METHODS", unset = paste(names(algos), collapse = ","))
requested_methods <- trimws(strsplit(requested_raw, ",")[[1]])
requested_methods <- requested_methods[nzchar(requested_methods)]

resolve_method <- function(m, available) {
  if (m %in% available) return(m)
  if (m == "T.test" && "T-test" %in% available) return("T-test")
  if (m == "T-test" && "T.test" %in% available) return("T.test")
  NA_character_
}

resolved_methods <- vapply(requested_methods, resolve_method, character(1), available = names(algos))
missing_methods <- requested_methods[is.na(resolved_methods)]
if (length(missing_methods) > 0) {
  stop("Unknown methods: ", paste(missing_methods, collapse = ", "),
       "\nAvailable methods in runDE.R algos: ", paste(names(algos), collapse = ", "))
}
algos_run <- algos[resolved_methods]
names(algos_run) <- resolved_methods
expected_output_methods <- canonical_method(names(algos_run))

contains_try_error <- function(x) {
  if (inherits(x, "try-error")) return(TRUE)
  if (inherits(x, "error")) return(TRUE)
  if (is.list(x)) return(any(vapply(x, contains_try_error, logical(1))))
  FALSE
}

is_good_rds <- function(path, expected_len = NULL, reject_try_error = FALSE) {
  if (!file.exists(path)) return(FALSE)
  x <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(x, "error")) return(FALSE)
  if (!is.null(expected_len) && length(x) != expected_len) return(FALSE)
  if (reject_try_error && contains_try_error(x)) return(FALSE)
  TRUE
}

safe_saveRDS <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp.", Sys.getpid())
  saveRDS(object, tmp)
  if (file.exists(file)) file.remove(file)
  ok <- file.rename(tmp, file)
  if (!ok) {
    file.copy(tmp, file, overwrite = TRUE)
    file.remove(tmp)
  }
  invisible(TRUE)
}

run_one_method <- function(y, method_name, method_fun) {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
  if (requireNamespace("dplyr", quietly = TRUE)) suppressPackageStartupMessages(library(dplyr))

  keep <- rowSums(y$counts >= 5) >= 2
  y <- y[keep, , keep.lib.sizes = FALSE]
  genes <- rownames(y$counts)

  res <- tryCatch(method_fun(y), error = function(e) e)
  if (inherits(res, "error")) {
    return(structure(paste0("Method failed: ", method_name, " | ", conditionMessage(res)),
                     class = "try-error", condition = res))
  }

  if (!is.data.frame(res) || !"PValue" %in% colnames(res)) {
    return(structure(paste0("Method did not return data.frame with PValue: ", method_name),
                     class = "try-error"))
  }

  pv <- rep(1, length(genes))
  names(pv) <- genes

  common <- intersect(genes, rownames(res))
  if (length(common) > 0) {
    pv[common] <- suppressWarnings(as.numeric(res[common, "PValue"]))
  }

  bad <- is.na(pv) | is.nan(pv) | is.infinite(pv) | pv < 0 | pv > 1
  pv[bad] <- 1

  out <- data.frame(pv, row.names = genes)
  colnames(out) <- canonical_method(method_name)
  out
}

batch_file_path <- function(batch_dir, start_i, end_i) {
  file.path(batch_dir, sprintf("batch_%03d_%03d.rds", start_i, end_i))
}

method_final_file <- function(output_base, method_name, prefix) {
  file.path(output_base, sanitize_method(method_name), "datp", paste0(prefix, ".rds"))
}

method_batch_dir <- function(output_base, method_name, prefix) {
  file.path(output_base, sanitize_method(method_name), "batch", prefix)
}

generate_fpr_counts_if_needed <- function(dge, work_dir) {
  counts_dir <- file.path(work_dir, "counts")
  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

  if (FPR_MODE == "normal_only") {
    source_counts <- dge$counts[, dge$samples$group == "Normal", drop = FALSE]
  } else if (FPR_MODE == "all_random_labels") {
    source_counts <- dge$counts
  } else {
    stop("Unsupported FPR_MODE: ", FPR_MODE)
  }

  for (s in NSAMPLE) {
    count_file <- file.path(counts_dir, sprintf("%s-SCLC.rds", s))
    if (is_good_rds(count_file, B_FPR)) {
      next
    }

    if (ncol(source_counts) < 2 * s) {
      warning("Skipping FPR n=", s, ": insufficient samples.")
      next
    }

    set.seed(SEED_VALUE)
    myseed <- sample(1:2000000, B_FPR)
    counts_list <- vector("list", B_FPR)

    for (i in seq_along(myseed)) {
      set.seed(myseed[i])
      inds <- sample(colnames(source_counts), 2 * s, replace = FALSE)
      group <- rep(c(1, 2), each = s)
      dat <- source_counts[, inds, drop = FALSE]
      counts_list[[i]] <- edgeR::DGEList(counts = dat, group = group)
    }

    safe_saveRDS(counts_list, count_file)
  }
}

generate_cat_counts_if_needed <- function(dge, work_dir) {
  counts_dir <- file.path(work_dir, "counts")
  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

  cancer_idx <- which(dge$samples$group == "Cancer")
  normal_idx <- which(dge$samples$group == "Normal")

  set.seed(SEED_VALUE)
  myseed <- sample(1:9999999, B_CAT)

  for (n in NSAMPLE) {
    if (length(cancer_idx) < 2 * n || length(normal_idx) < 2 * n) {
      warning("Skipping CAT n=", n, ": insufficient samples.")
      next
    }

    for (subset_name in c("subset1", "subset2")) {
      count_file <- file.path(counts_dir, sprintf("%s-%s.rds", n, subset_name))
      if (is_good_rds(count_file, B_CAT)) {
        next
      }

      counts_list <- lapply(myseed, function(seed) {
        set.seed(seed)
        index.a <- sample(cancer_idx, 2 * n)
        index.b <- sample(normal_idx, 2 * n)

        chosen <- if (subset_name == "subset1") {
          c(head(index.a, n = n), head(index.b, n = n))
        } else {
          c(tail(index.a, n = n), tail(index.b, n = n))
        }

        counts <- dge$counts[, chosen, drop = FALSE]
        group <- c(rep("a", n), rep("b", n))
        edgeR::DGEList(counts = counts, group = group)
      })

      safe_saveRDS(counts_list, count_file)
    }
  }
}

combine_batches_to_method_final <- function(batch_dir, final_file, B) {
  batch_files <- sort(list.files(batch_dir, pattern = "^batch_.*\\.rds$", full.names = TRUE))
  if (length(batch_files) == 0) {
    return(FALSE)
  }

  out <- vector("list", B)
  filled <- rep(FALSE, B)

  for (bf in batch_files) {
    x <- readRDS(bf)
    idx <- attr(x, "replicate_index")
    if (is.null(idx) || length(x) != length(idx)) {
      warning("Bad replicate_index in batch file: ", bf)
      next
    }
    out[idx] <- x
    filled[idx] <- TRUE
  }

  if (!all(filled)) {
    return(FALSE)
  }

  safe_saveRDS(out, final_file)
  TRUE
}

run_method_batches_for_count_file <- function(count_file, method_name, method_fun, output_base, expected_len) {
  prefix <- sub("\\.rds$", "", basename(count_file))
  final_file <- method_final_file(output_base, method_name, prefix)

  if (is_good_rds(final_file, expected_len = expected_len, reject_try_error = TRUE)) {
    return(TRUE)
  }

  dat_list <- readRDS(count_file)
  if (length(dat_list) != expected_len) {
    warning("Count file length differs from expected_len: ", count_file,
            " length=", length(dat_list), " expected=", expected_len)
  }

  batch_dir <- method_batch_dir(output_base, method_name, prefix)
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

  B <- length(dat_list)
  starts <- seq(1, B, by = BATCH_SIZE)
  ends <- pmin(starts + BATCH_SIZE - 1, B)

  for (k in seq_along(starts)) {
    idx <- starts[k]:ends[k]
    batch_file <- batch_file_path(batch_dir, starts[k], ends[k])

    if (is_good_rds(batch_file, expected_len = length(idx), reject_try_error = TRUE)) {
      next
    }

    t0 <- Sys.time()

    ans <- parallel::mclapply(
      idx,
      function(i) run_one_method(dat_list[[i]], method_name, method_fun),
      mc.cores = min(WORKERS, length(idx))
    )

    if (contains_try_error(ans)) {
      err_dir <- file.path(batch_dir, "errors")
      dir.create(err_dir, recursive = TRUE, showWarnings = FALSE)
      err_file <- file.path(err_dir, sprintf("batch_%03d_%03d_error_%s.txt",
                                             starts[k], ends[k], format(Sys.time(), "%Y%m%d_%H%M%S")))
      writeLines(capture.output(str(ans, max.level = 3)), err_file)
      next
    }

    attr(ans, "replicate_index") <- idx
    safe_saveRDS(ans, batch_file)

    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  combine_batches_to_method_final(batch_dir, final_file, B)
}

combine_methods_for_count_file <- function(work_dir, prefix, expected_len) {
  by_method_dir <- file.path(work_dir, "datp_by_method")
  combined_dir <- file.path(work_dir, "datp_method_combined")
  original_dir <- file.path(work_dir, "datp")

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(original_dir, recursive = TRUE, showWarnings = FALSE)

  method_files <- vapply(names(algos_run), function(m) method_final_file(by_method_dir, m, prefix), character(1))
  names(method_files) <- names(algos_run)

  existing <- method_files[file.exists(method_files)]
  if (length(existing) == 0) {
    return(FALSE)
  }

  method_lists <- lapply(existing, readRDS)
  ok_len <- vapply(method_lists, function(x) length(x) == expected_len, logical(1))
  if (!all(ok_len)) {
    warning("Some method final files have unexpected length for ", prefix, ": ",
            paste(names(existing)[!ok_len], collapse = ", "))
  }

  method_lists <- method_lists[ok_len]
  if (length(method_lists) == 0) return(FALSE)

  combined <- vector("list", expected_len)

  for (i in seq_len(expected_len)) {
    dfs <- lapply(method_lists, `[[`, i)
    all_genes <- unique(unlist(lapply(dfs, rownames), use.names = FALSE))

    merged <- data.frame(row.names = all_genes)

    for (df in dfs) {
      m <- colnames(df)[1]
      pv <- rep(1, length(all_genes))
      names(pv) <- all_genes
      pv[rownames(df)] <- suppressWarnings(as.numeric(df[[1]]))
      bad <- is.na(pv) | is.nan(pv) | is.infinite(pv) | pv < 0 | pv > 1
      pv[bad] <- 1
      merged[[m]] <- pv
    }

    present <- intersect(expected_output_methods, colnames(merged))
    extra <- setdiff(colnames(merged), present)
    merged <- merged[, c(present, extra), drop = FALSE]
    combined[[i]] <- merged
  }

  combined_file <- file.path(combined_dir, paste0(prefix, ".rds"))
  safe_saveRDS(combined, combined_file)

  if (WRITE_ORIGINAL_DATP) {
    original_file <- file.path(original_dir, paste0(prefix, ".rds"))
    safe_saveRDS(combined, original_file)
  }

  TRUE
}

compute_cat_concordance <- function(work_dir) {

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)

  setwd(work_dir)
  source(CONC_FILE)

  dir.create("./conc", recursive = TRUE, showWarnings = FALSE)

  datp_files <- list.files("./datp", full.names = TRUE, pattern = "\\.rds$")
  if (length(datp_files) == 0) {
    return(invisible(FALSE))
  }

  for (file in datp_files) {
    target <- gsub("datp", "conc", file)

    if (is_good_rds(target, expected_len = B_CAT, reject_try_error = TRUE)) {
      next
    }

    get_conc(file)
  }

  for (n in NSAMPLE) {
    target <- sprintf("./conc/%s-1vs2.rds", n)

    if (is_good_rds(target, expected_len = B_CAT, reject_try_error = TRUE)) {
      next
    }

    get_conc_w(n)
  }

  invisible(TRUE)
}

run_fpr_module <- function() {
  work_dir <- file.path(BASE_DIR, "FPR_SCLC")
  counts_dir <- file.path(work_dir, "counts")
  by_method_dir <- file.path(work_dir, "datp_by_method")

  generate_fpr_counts_if_needed(dge, work_dir)

  count_files <- file.path(counts_dir, sprintf("%s-SCLC.rds", NSAMPLE))
  count_files <- count_files[file.exists(count_files)]

  for (method_name in names(algos_run)) {

    for (cf in count_files) {
      run_method_batches_for_count_file(
        count_file = cf,
        method_name = method_name,
        method_fun = algos_run[[method_name]],
        output_base = by_method_dir,
        expected_len = B_FPR
      )

      combine_methods_for_count_file(
        work_dir = work_dir,
        prefix = sub("\\.rds$", "", basename(cf)),
        expected_len = B_FPR
      )
    }
  }
}

run_cat_module <- function() {
  work_dir <- file.path(BASE_DIR, "CAT_SCLC")
  counts_dir <- file.path(work_dir, "counts")
  by_method_dir <- file.path(work_dir, "datp_by_method")

  generate_cat_counts_if_needed(dge, work_dir)

  count_files <- character(0)
  for (n in NSAMPLE) {
    count_files <- c(
      count_files,
      file.path(counts_dir, sprintf("%s-subset1.rds", n)),
      file.path(counts_dir, sprintf("%s-subset2.rds", n))
    )
  }
  count_files <- count_files[file.exists(count_files)]

  for (method_name in names(algos_run)) {

    for (cf in count_files) {
      run_method_batches_for_count_file(
        count_file = cf,
        method_name = method_name,
        method_fun = algos_run[[method_name]],
        output_base = by_method_dir,
        expected_len = B_CAT
      )

      combine_methods_for_count_file(
        work_dir = work_dir,
        prefix = sub("\\.rds$", "", basename(cf)),
        expected_len = B_CAT
      )
    }
  }

  if (COMPUTE_CAT_CONC && WRITE_ORIGINAL_DATP) {
    compute_cat_concordance(work_dir)
  } else if (COMPUTE_CAT_CONC && !WRITE_ORIGINAL_DATP) {
  }
}

if (grepl("FPR", RUN_MODULES)) run_fpr_module()
if (grepl("CAT", RUN_MODULES)) run_cat_module()

