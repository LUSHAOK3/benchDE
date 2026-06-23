library(data.table)
library(tidyverse)
library(parallel)

get_conc <- function(file) {
  data <- readRDS(file)
  temp <- mclapply(data, function(dat) {
    nfeatures <- nrow(dat)
    k_max <- round(0.1 * nrow(dat))
    
    max_iters <- ncol(dat) * ncol(dat) * k_max
    df_res <- data.table(
      method1 = character(max_iters),
      method2 = character(max_iters),
      rank = integer(max_iters),
      conc = numeric(max_iters),
      comparision = character(max_iters),
      nfeatures = numeric(max_iters)
    )
    
    counter <- 1
    
    for (i in 1:ncol(dat)) {
      for (j in 1:ncol(dat)) {
        if (i != j) {      
          
          i.sorted <- names(sort(setNames(dat[, i], rownames(dat)), decreasing = FALSE))
          j.sorted <- names(sort(setNames(dat[, j], rownames(dat)), decreasing = FALSE))
          nfeatures <- nrow(dat)
          
          for (k in 1:k_max) {
            
            i.top <- i.sorted[1:k]
            j.top <- j.sorted[1:k]
            
            inters <- length(intersect(i.top, j.top))
            conc <- inters / k
            
            tmp <- data.table(
              method1 = colnames(dat)[i],
              method2 = colnames(dat)[j],
              rank = k,
              conc = conc,
              comparision = basename(file),
              nfeatures = nfeatures
            )
            df_res[counter, ] <- tmp
            
            counter <- counter + 1
          }
        }
      }
    }
    df_res <- df_res[1:(counter - 1), ]
    
  }, mc.cores = 10)
  
  saveRDS(temp, file = gsub("datp", "conc", file))
}

get_conc_w <- function(nsample) {
  
  data1 <- readRDS(sprintf("./datp/%s-subset1.rds", nsample))
  data2 <- readRDS(sprintf("./datp/%s-subset2.rds", nsample))
  temp <- mclapply(1:10, function(n) {
    dat1 <- data1[[n]]
    dat2 <- data2[[n]]
    
    dat1[is.na(dat1)] <- 1
    dat2[is.na(dat2)] <- 1
    dat1 <- data.frame(dat1)
    dat2 <- data.frame(dat2)
    
    nfeatures <- max(nrow(dat1), nrow(dat2))
    k_max <- max(nrow(dat1), nrow(dat2)) * 0.1
    
    sorted_a <- lapply(dat1, function(col) names(sort(setNames(col, rownames(dat1)), decreasing = FALSE)))
    sorted_b <- lapply(dat2, function(col) names(sort(setNames(col, rownames(dat2)), decreasing = FALSE)))
    
    list_res <- vector("list", ncol(dat1) * 1000)
    
    counter <- 1
    
    for (i in 1:ncol(dat1)) {
      for (k in 1:k_max) {
        a.top <- sorted_a[[i]][1:k]
        b.top <- sorted_b[[i]][1:k]
        
        inters <- length(intersect(a.top, b.top))
        conc <- inters / k
        
        list_res[[counter]] <- data.frame(
          method1 = colnames(dat1)[i],
          method2 = colnames(dat2)[i],
          rank = k,
          conc = conc,
          comparision = paste0(nsample, "n-1vs2"),
          nfeatures = nfeatures
        )
        
        counter <- counter + 1
      }
    }
    df_res <- do.call(rbind, list_res)
  }, mc.cores = 10)
  
  saveRDS(temp, file = sprintf("./conc/%s-1vs2.rds", nsample))
}
