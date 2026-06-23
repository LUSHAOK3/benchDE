library(zoo)
library(tidyverse)
library(data.table)
library(edgeR)
library(limma)

# Poisson Strategy 2: Poisson + random outlier simulation.
# This script mimics sim2_function.R as closely as possible.
#
# Original sim2:
#   1) Generate NB counts.
#   2) Randomly select RO.prop% count entries.
#   3) Multiply selected entries by runif(0, 10), then round.
#
# This Poisson version:
#   1) Generate Poisson counts using lambda = mean.
#   2) Apply the same random outlier mechanism.
#
# Only the count-generating distribution changes from NB to Poisson.

get_param <- function(dge){

  dge$counts <- dge$counts[rowMeans(dge$counts) > 7, ]  
  cond1 <- unique(dge$sample$group)[1]
  cond2 <- unique(dge$sample$group)[2]  

  dge.c1 <- DGEList(counts = dge$counts[, dge$sample$group == cond1], group = rep(1, sum(dge$sample$group == cond1)))
  dge.c2 <- DGEList(counts = dge$counts[, dge$sample$group == cond2], group = rep(2, sum(dge$sample$group == cond2)))

  mean.c1 <- apply(dge.c1$counts, 1, mean)
  mean.c2 <- apply(dge.c2$counts, 1, mean)
  mean.total <- apply(dge$counts, 1, mean)

  return(list(mean.c1 = mean.c1, mean.c2 = mean.c2, mean.total = mean.total))
}

create_countmat <- function(param_list, gs, mode = mode, random_sampling = FALSE, 
                            s = s, n.var = 10000, RO.prop = RO.prop,
                            Large_sample = FALSE, fixedfold = FALSE){

  mean.c1 <- param_list$mean.c1
  mean.c2 <- param_list$mean.c2
  mean.total <- param_list$mean.total

  n.diffexp <- 3000
  fraction.upregulated <- 0.5

  sample.mean1 <- mean.total[gs]
  sample.mean2 <- sample.mean1
  truth <- rep(FALSE, n.var)

  # Keep the fold-change setting used by sim2_function.R.
  foldchage0 <- ifelse(s <= 10, 2, 1.25)

  if(n.diffexp != 0){
    upindex <- 1:round(n.diffexp * fraction.upregulated)
    dnindex <- round(n.diffexp * fraction.upregulated + 1):n.diffexp
    truth[c(upindex, dnindex)] <- TRUE

    factor1 <- foldchage0 + rexp(n = length(upindex), rate = 1) 
    factor2 <- foldchage0 + rexp(n = length(dnindex), rate = 1)

    sample.mean2[upindex] <- sample.mean2[upindex] * factor1
    sample.mean2[dnindex] <- sample.mean2[dnindex] / factor2
  } 

  sample.mean1[is.na(sample.mean1) | sample.mean1 < 0] <- 0
  sample.mean2[is.na(sample.mean2) | sample.mean2 < 0] <- 0

  counts <- matrix(nrow = n.var, ncol = 2 * s)

  if(random_sampling == TRUE){
    rand1 <- runif(s, min = 0.7, max = 1.3)
    rand2 <- runif(s, min = 0.7, max = 1.3)
  } else {
    rand1 <- rep(1, s)
    rand2 <- rep(1, s)
  }

  for(i in 1:n.var){
    counts[i, 1:s] <- sapply(rand1, FUN = function(x) rpois(1, lambda = sample.mean2[i] * x))
    counts[i, (s + 1):(2 * s)] <- sapply(rand2, FUN = function(x) rpois(1, lambda = sample.mean1[i] * x))
  }

  ### Random Outlier, same mechanism as sim2_function.R
  if(mode == "R"){
    RO <- matrix(runif(n.var * 2 * s , min = 0, max = 100), nrow = n.var, ncol = 2 * s)
    index.outlier <- which(RO < RO.prop)
    if(length(index.outlier) > 0){
      counts[index.outlier] <- counts[index.outlier] * runif(n = length(index.outlier), min = 0, max = 10)
      counts <- round(counts)
    }
  }

  rownames(counts) <- paste("g", 1:n.var, sep = "")
  colnames(counts) <- seq(2 * s)
  dat_alt <- DGEList(counts = counts, group = c(rep(1, s), rep(2, s)))

  params <- list(sample.mean1 = sample.mean1, sample.mean2 = sample.mean2,
                 RO.prop = RO.prop)

  res <- list("dge" = dat_alt, "truth" = truth, "params" = params)                              
  return(res)
}

simulateOneSplit <- function(X, rseed, param = param, RO.prop = RO.prop, s = s, n.var = 10000){

  n.total <- length(param[[1]])
  set.seed(as.numeric(X) * as.numeric(rseed))
  gs <- sample(1:n.total, n.var,
               replace = ifelse(n.total < n.var , TRUE , FALSE))

  dat <- create_countmat(param, gs, s = s, n.var = n.var, RO.prop = RO.prop, mode = "R")
  return(dat)
}
