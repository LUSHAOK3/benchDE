# sim3 function.

library(zoo)
library(tidyverse)
library(data.table)
library(edgeR)
library(limma)

get_param <- function(dge){

  dge$counts <- dge$counts[rowMeans(dge$counts) > 7, ]  
  cond1 <- unique(dge$sample$group)[1]
  cond2 <- unique(dge$sample$group)[2]  

  dge.c1 <- DGEList(counts = dge$counts[, dge$sample$group == cond1], group = rep(1, sum(dge$sample$group == cond1)))
  dge.c2 <- DGEList(counts = dge$counts[, dge$sample$group == cond2], group = rep(2, sum(dge$sample$group == cond2)))

  dge.c1 <- calcNormFactors(dge.c1)
  dge.c1 <- estimateCommonDisp(dge.c1)
  dge.c1 <- estimateTagwiseDisp(dge.c1)
  disp.c1 <- dge.c1$tagwise.dispersion
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

  return(list(disp.c1 = disp.c1, disp.c2 = disp.c2, disp.total = disp.total, 
              mean.c1 = mean.c1, mean.c2 = mean.c2, mean.total = mean.total))
}

getDisp <- function(mean, mean.condition, disp.condition){
  pool <- disp.condition[which(mean.condition > (mean - 20) & mean.condition < (mean + 20))]
  if(length(pool) == 0){
    value <- disp.condition[which.min(abs(mean.condition - mean))]
  } else {
    value <- sample(pool, 1)
  }
  return(value)
}

create_countmat <- function(param_list, gs, dispType = dispType, mode = mode,
                            random_sampling = FALSE, s = s, n.var = 10000,
                            RO.prop = RO.prop, Noise = 0, Large_sample = FALSE,
                            fixedfold = FALSE){

  if (is.null(Noise) || length(Noise) != 1 || is.na(Noise)) {
    stop("Noise must be a single numeric value.")
  }
  if (Noise < 0 || Noise >= 1) {
    stop("Noise must be in [0, 1). Suggested values: 0, 0.2, 0.4, 0.6, 0.8.")
  }

  disp.c1 <- param_list$disp.c1
  disp.c2 <- param_list$disp.c2
  disp.total <- param_list$disp.total

  mean.c1 <- param_list$mean.c1
  mean.c2 <- param_list$mean.c2
  mean.total <- param_list$mean.total

  n.diffexp <- 3000
  fraction.upregulated <- 0.5

  sample.mean1 <- mean.total[gs]
  sample.mean2 <- sample.mean1
  truth <- rep(FALSE, n.var)

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

  if(dispType == 'different'){
    sample.disp1 <- sapply(sample.mean1, FUN = getDisp, mean.condition = mean.c1, disp.condition = disp.c1, simplify = TRUE, USE.NAMES = FALSE)
    sample.disp2 <- sapply(sample.mean2, FUN = getDisp, mean.condition = mean.c2, disp.condition = disp.c2, simplify = TRUE, USE.NAMES = FALSE)
  } else if(dispType == 'same'){
    sample.disp1 <- disp.total[gs]
    sample.disp2 <- sample.disp1
  } else {
    stop("dispType must be 'different' or 'same'.")
  }

  counts <- matrix(nrow = n.var, ncol = 2 * s)

  if(mode == "Noise"){
    noise_factor2 <- matrix(runif(n.var * s, min = 1 - Noise, max = 1 + Noise), nrow = n.var, ncol = s)
    noise_factor1 <- matrix(runif(n.var * s, min = 1 - Noise, max = 1 + Noise), nrow = n.var, ncol = s)

    for(i in 1:n.var){
      mu2 <- pmax(sample.mean2[i] * noise_factor2[i, ], 0)
      mu1 <- pmax(sample.mean1[i] * noise_factor1[i, ], 0)
      counts[i, 1:s] <- rnbinom(s, 1 / sample.disp2[i], mu = mu2)
      counts[i, (s + 1):(2 * s)] <- rnbinom(s, 1 / sample.disp1[i], mu = mu1)
    }
  } else {
    stop("This sim3_function.R is designed for mode = 'Noise'.")
  }

  rownames(counts) <- paste("g", 1:n.var, sep = "")
  colnames(counts) <- seq(2 * s)
  dat_alt <- DGEList(counts = counts, group = c(rep(1, s), rep(2, s)))

  params <- list(sample.disp1 = sample.disp1, sample.disp2 = sample.disp2,
                 sample.mean1 = sample.mean1, sample.mean2 = sample.mean2,
                 Noise = Noise)

  res <- list("dge" = dat_alt, "truth" = truth, "params" = params)                              
  return(res)
}

simulateOneSplit <- function(X, rseed, param = param, Noise = Noise, s = s, n.var = 10000){

  n.total <- length(param[[1]])
  set.seed(as.numeric(X) * as.numeric(rseed))
  gs <- sample(1:n.total, n.var,
               replace = ifelse(n.total < n.var , TRUE , FALSE))

  dat <- create_countmat(param, gs, s = s, dispType = 'different',
                         n.var = n.var, RO.prop = 1, Noise = Noise,
                         mode = "Noise")
  return(dat)
}
