library(zoo)
library(tidyverse)
library(data.table)
library(edgeR)
library(limma)

# This function is adapted from the simulation framework proposed by Baik et al. (2020) for benchmarking RNA-seq differential expression methods using spike-in and simulation data.
# Reference: Baik B, Yoon S, Nam D. Benchmarking RNA-seq differential expression analysis methods using spike-in and simulation data. PLOS ONE. 2020;15:e0232271.
get_param <- function(dge) {
  dge$counts <- dge$counts[rowMeans(dge$counts) > 7, ]  
  cond1 <- unique(dge$sample$group)[1]
  cond2 <- unique(dge$sample$group)[2]  
  
  dge.c1 <- DGEList(counts = dge$counts[, dge$sample$group == cond1], group = rep(1, sum(dge$sample$group == cond1)))
  dge.c2 <- DGEList(counts = dge$counts[, dge$sample$group == cond2], group = rep(2, sum(dge$sample$group == cond2)))
  
  dge.c1 <- calcNormFactors(dge.c1)
  dge.c1 <- estimateCommonDisp(dge.c1)
  dge.c1 <- estimateTagwiseDisp(dge.c1)
  disp.c1 <- dge.c1$tagwise.dispersion # Dispersion
  mean.c1 <- apply(dge.c1$counts, 1, mean)
  
  dge.c2 <- calcNormFactors(dge.c2)
  dge.c2 <- estimateCommonDisp(dge.c2)
  dge.c2 <- estimateTagwiseDisp(dge.c2)
  disp.c2 <- dge.c2$tagwise.dispersion
  mean.c2 <- apply(dge.c2$counts, 1, mean)
  
  mean.total <- apply(dge$counts, 1, mean)
  dge <- calcNormFactors(dge)
  dge <- estimateCommonDisp(dge)
  dge <- estimateTagwiseDisp(dge)
  disp.total <- dge$tagwise.dispersion
  
  return(list(disp.c1 = disp.c1, disp.c2 = disp.c1, disp.total = disp.total, 
              mean.c1 = mean.c1, mean.c2 = mean.c2, mean.total = mean.total))
}

# This function is adapted from the simulation framework proposed by Baik et al. (2020) for benchmarking RNA-seq differential expression methods using spike-in and simulation data.
# Reference: Baik B, Yoon S, Nam D. Benchmarking RNA-seq differential expression analysis methods using spike-in and simulation data. PLOS ONE. 2020;15:e0232271.
getDisp <- function(mean, mean.condition, disp.condition) {
  pool <- disp.condition[which(mean.condition > (mean - 20) & mean.condition < (mean + 20))]
  if (length(pool) == 0) {
    value <- disp.condition[which.min(abs(mean.condition - mean))]
  } else {
    value <- sample(pool, 1)
  }
  return(value)
}

# This function is adapted from the simulation framework proposed by Baik et al. (2020) for benchmarking RNA-seq differential expression methods using spike-in and simulation data.
# Reference: Baik B, Yoon S, Nam D. Benchmarking RNA-seq differential expression analysis methods using spike-in and simulation data. PLOS ONE. 2020;15:e0232271.
create_countmat <- function(param_list, gs, dispType = dispType, mode = mode, random_sampling = FALSE, 
                            s = s, n.var = 10000, RO.prop = 1, Large_sample = FALSE, fixedfold = FALSE) {
  disp.c1 <- param_list$disp.c1
  disp.c2 <- param_list$disp.c1
  disp.total <- param_list$disp.total
  
  mean.c1 <- param_list$mean.c1
  mean.c2 <- param_list$mean.c2
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
    
    factor1 <- 2 + rexp(n = length(upindex), rate = 1) 
    factor2 <- 2 + rexp(n = length(dnindex), rate = 1)
    
    sample.mean2[upindex] <- sample.mean2[upindex] * factor1
    sample.mean2[dnindex] <- sample.mean2[dnindex] / factor2
  }  
  
  if (dispType == 'different') {
    sample.disp1 <- sapply(sample.mean1, FUN = getDisp, mean.condition = mean.c1, disp.condition = disp.c1, simplify = TRUE, USE.NAMES = FALSE)
    sample.disp2 <- sapply(sample.mean2, FUN = getDisp, mean.condition = mean.c2, disp.condition = disp.c2, simplify = TRUE, USE.NAMES = FALSE)
  } else if (dispType == 'same') {
    sample.disp1 <- disp.total[gs]
    sample.disp2 <- sample.disp1
  }
  
  counts <- matrix(nrow = n.var, ncol = 2 * s)
  
  if (mode == "OS") {
    for (i in 1:n.var) {
      counts[i, 1:floor(s / 10)] <- rnbinom(floor(s / 10), 1 / (3 * sample.disp2[i]), mu = sample.mean2[i])
      counts[i, (floor(s / 10) + 1):s] <- rnbinom((s - floor(s / 10)), 1 / sample.disp2[i], mu = sample.mean2[i])
      counts[i, (s + 1):(2 * s)] <- rnbinom(s, 1 / sample.disp1[i], mu = sample.mean1[i])  
    }
  } else {
    if (random_sampling == TRUE) {
      rand1 <- runif(s, min = 0.7, max = 1.3)
      rand2 <- runif(s, min = 0.7, max = 1.3)
    } else {
      rand1 <- rep(1, s)
      rand2 <- rep(1, s)
    }
    for (i in 1:n.var) {
      counts[i, 1:s] <- sapply(rand1, FUN = function(x) rnbinom(1, 1 / sample.disp2[i], mu = sample.mean2[i] * x))
      counts[i, (s + 1):(2 * s)] <- sapply(rand2, FUN = function(x) rnbinom(1, 1 / sample.disp1[i], mu = sample.mean1[i] * x))
    }
  }
  
  ### Random Outlier
  if (mode == "R") {
    RO <- matrix(runif(n.var * 2 * s, min = 0, max = 100), nrow = n.var, ncol = 2 * s)
    index.outlier <- which(RO < RO.prop)
    counts[index.outlier] <- counts[index.outlier] * runif(n = length(index.outlier), min = 5, max = 10)
    counts <- round(counts)
  }
  
  rownames(counts) <- paste("g", 1:n.var, sep = "")
  colnames(counts) <- seq(2 * s)
  dat_alt <- DGEList(counts = counts, group = c(rep(1, s), rep(2, s)))
  
  params <- list(sample.disp1 = sample.disp1, sample.disp2 = sample.disp2,
                 sample.mean1 = sample.mean1, sample.mean2 = sample.mean2)
  
  res <- list("dge" = dat_alt, "truth" = truth, "params" = params)                              
  return(res)
}

simulateOneSplit <- function(X, rseed, param = param, s = s, n.var = 10000) {
  n.total <- length(param[[1]])
  set.seed(as.numeric(X) * as.numeric(rseed))
  gs <- sample(1:n.total, n.var,
               replace = ifelse(n.total < n.var, TRUE, FALSE))
  
  dat <- create_countmat(param, gs, s = s, dispType = 'different', n.var = n.var, RO.prop = 1, mode = "NB")
}
