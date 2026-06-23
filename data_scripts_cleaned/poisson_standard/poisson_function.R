library(zoo)
library(tidyverse)
library(data.table)
library(edgeR)
library(limma)

# Poisson simulation function for revision.
# This script mimics sim1_function.R as closely as possible.
# Difference from NB simulation:
#   NB:      counts ~ NB(mu, dispersion)
#   Poisson: counts ~ Poisson(lambda = mu)
# Therefore dispersion is not used to generate counts here.

get_param <- function(dge) {
  dge$counts <- dge$counts[rowMeans(dge$counts) > 7, ]
  cond1 <- unique(dge$sample$group)[1]
  cond2 <- unique(dge$sample$group)[2]

  dge.c1 <- DGEList(
    counts = dge$counts[, dge$sample$group == cond1],
    group = rep(1, sum(dge$sample$group == cond1))
  )
  dge.c2 <- DGEList(
    counts = dge$counts[, dge$sample$group == cond2],
    group = rep(2, sum(dge$sample$group == cond2))
  )

  mean.c1 <- apply(dge.c1$counts, 1, mean)
  mean.c2 <- apply(dge.c2$counts, 1, mean)
  mean.total <- apply(dge$counts, 1, mean)

  return(list(
    mean.c1 = mean.c1,
    mean.c2 = mean.c2,
    mean.total = mean.total
  ))
}

create_countmat <- function(param_list, gs, s = s, n.var = 10000,
                            fixedfold = FALSE) {
  mean.total <- param_list$mean.total

  n.diffexp <- 3000
  fraction.upregulated <- 0.5

  sample.mean1 <- mean.total[gs]
  sample.mean2 <- sample.mean1
  truth <- rep(FALSE, n.var)

  if (n.diffexp != 0) {
    upindex <- 1:round(n.diffexp * fraction.upregulated)
    dnindex <- round(n.diffexp * fraction.upregulated + 1):n.diffexp
    truth[c(upindex, dnindex)] <- TRUE

    # Keep the same fold-change generation form as sim1_function.R.
    factor1 <- 2 + rexp(n = length(upindex), rate = 1)
    factor2 <- 2 + rexp(n = length(dnindex), rate = 1)

    sample.mean2[upindex] <- sample.mean2[upindex] * factor1
    sample.mean2[dnindex] <- sample.mean2[dnindex] / factor2
  }

  sample.mean1[is.na(sample.mean1) | sample.mean1 < 0] <- 0
  sample.mean2[is.na(sample.mean2) | sample.mean2 < 0] <- 0

  counts <- matrix(nrow = n.var, ncol = 2 * s)

  for (i in 1:n.var) {
    counts[i, 1:s] <- rpois(s, lambda = sample.mean2[i])
    counts[i, (s + 1):(2 * s)] <- rpois(s, lambda = sample.mean1[i])
  }

  rownames(counts) <- paste("g", 1:n.var, sep = "")
  colnames(counts) <- seq(2 * s)
  dat_alt <- DGEList(counts = counts, group = c(rep(1, s), rep(2, s)))

  params <- list(
    sample.mean1 = sample.mean1,
    sample.mean2 = sample.mean2
  )

  res <- list("dge" = dat_alt, "truth" = truth, "params" = params)
  return(res)
}

simulateOneSplit <- function(X, rseed = 8848, param = param,
                             s = s, n.var = 10000) {
  n.total <- length(param[[1]])
  set.seed(as.numeric(X) * as.numeric(rseed))
  gs <- sample(
    1:n.total, n.var,
    replace = ifelse(n.total < n.var, TRUE, FALSE)
  )

  dat <- create_countmat(param, gs, s = s, n.var = n.var)
  return(dat)
}
