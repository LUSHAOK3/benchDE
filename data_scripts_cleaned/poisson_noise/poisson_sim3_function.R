library(zoo)
library(tidyverse)
library(data.table)
library(edgeR)
library(limma)

# Poisson Strategy 3: Poisson + noise simulation.
# This script is the Poisson counterpart of the restored NB Strategy 3.
#
# Noise definition:
#   For each gene-by-sample observation, a multiplicative noise factor is sampled:
#       e_ij ~ Uniform(1 - Noise, 1 + Noise)
#   The Poisson mean is perturbed as:
#       lambda'_ij = lambda_g * e_ij
#   Counts are then generated as:
#       Y_ij ~ Poisson(lambda'_ij)
#
# Noise = 0 recovers the baseline Poisson simulation.

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

create_countmat <- function(param_list, gs, mode = "Noise", random_sampling = FALSE,
                            s = s, n.var = 10000, Noise = Noise,
                            Large_sample = FALSE, fixedfold = FALSE){

  if (is.null(Noise) || length(Noise) != 1 || is.na(Noise)) {
    stop("Noise must be a single numeric value.")
  }
  if (Noise < 0 || Noise >= 1) {
    stop("Noise must be in [0, 1). Suggested values: 0, 0.2, 0.4, 0.6, 0.8.")
  }

  mean.c1 <- param_list$mean.c1
  mean.c2 <- param_list$mean.c2
  mean.total <- param_list$mean.total

  n.diffexp <- 3000
  fraction.upregulated <- 0.5

  sample.mean1 <- mean.total[gs]
  sample.mean2 <- sample.mean1
  truth <- rep(FALSE, n.var)

  # Keep the fold-change setting used by sim2/sim3 style scripts.
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

  if(mode == "Noise"){
    noise_factor2 <- matrix(runif(n.var * s, min = 1 - Noise, max = 1 + Noise), nrow = n.var, ncol = s)
    noise_factor1 <- matrix(runif(n.var * s, min = 1 - Noise, max = 1 + Noise), nrow = n.var, ncol = s)

    for(i in 1:n.var){
      lambda2 <- pmax(sample.mean2[i] * noise_factor2[i, ], 0)
      lambda1 <- pmax(sample.mean1[i] * noise_factor1[i, ], 0)
      counts[i, 1:s] <- rpois(s, lambda = lambda2)
      counts[i, (s + 1):(2 * s)] <- rpois(s, lambda = lambda1)
    }
  } else {
    stop("This poisson_sim3_function.R is designed for mode = 'Noise'.")
  }

  rownames(counts) <- paste("g", 1:n.var, sep = "")
  colnames(counts) <- seq(2 * s)
  dat_alt <- DGEList(counts = counts, group = c(rep(1, s), rep(2, s)))

  params <- list(sample.mean1 = sample.mean1, sample.mean2 = sample.mean2,
                 Noise = Noise)

  res <- list("dge" = dat_alt, "truth" = truth, "params" = params)                              
  return(res)
}

simulateOneSplit <- function(X, rseed, param = param, Noise = Noise, s = s, n.var = 10000){

  n.total <- length(param[[1]])
  set.seed(as.numeric(X) * as.numeric(rseed))
  gs <- sample(1:n.total, n.var,
               replace = ifelse(n.total < n.var , TRUE , FALSE))

  dat <- create_countmat(param, gs, s = s, n.var = n.var, Noise = Noise, mode = "Noise")
  return(dat)
}
