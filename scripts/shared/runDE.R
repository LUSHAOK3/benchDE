library(DESeq2)
library(DESeq)
library(edgeR)
library(limma)
library(DSS)
library(NOISeq)
library(ABSSeq)
library(NBPSeq)
library(ROTS)
library(samr)
library(dplyr)
# library(baySeq)
library(parallel)

runDESeq2 <- function(d) {
  counts <- as.data.frame(d$counts)
  dds <- DESeqDataSetFromMatrix(counts, d$samples, ~ group)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(DESeq2::results(dds))
  res$pvalue[is.na(res$pvalue)] <- 1
  res$padj[is.na(res$padj)] <- 1
  res <- dplyr::rename(res,PValue="pvalue")
  return(res)
}

runDESeq <- function(d){
  cds <- newCountDataSet(d$counts,d$samples$group)
  cds <- DESeq::estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, fitType="local")
  cons <- as.character(unique(d$samples$group)) 
  res <- nbinomTest(cds,cons[1],cons[2])
  res <- dplyr::rename(res,PValue="pval")
  rownames(res) <- res$id
  res <- res[-1]
  return(res)
}

runEdgeR.lrt <- function(d) {
  group=d$samples$group
  design <- model.matrix(~0+ group)
  dgel <- DGEList(d$counts)
  dgel <- calcNormFactors(dgel, method="TMM")
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit,contrast=c(-1,1))
  res <- as.data.frame(topTags(edger.lrt,n=nrow(dgel)))
  return(res)
}

runEdgeR.qlf <- function(d) {
  group=d$samples$group
  design <- model.matrix(~0+ group)
  dgel <- DGEList(d$counts)
  dgel <- calcNormFactors(dgel, method="TMM")
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmQLFit(dgel, design)
  edger.lrt <- glmQLFTest(edger.fit,contrast=c(-1,1))
  res <- as.data.frame(topTags(edger.lrt,n=nrow(dgel)))
  return(res)
}

runDSS <- function(d) {
  design <- d$samples[1]
  colnames(design) <- "designs" 
  cons <- unique(d$samples$group)
  seqData <- newSeqCountSet(d$counts, design)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  suppressWarnings({
    res <- waldTest(seqData, cons[2], cons[1])
  })
  res <- res[match(rownames(seqData),rownames(res)),]
  res <- dplyr::rename(res,PValue="pval")
  return(res)
}

runVoom <- function(d) {
  group=factor(d$samples$group)
  design <- model.matrix( ~ group)
  colnames(design) <- c(levels(group)[2],paste(rev(levels(group)),collapse = "-"))
  d <- calcNormFactors(d, method="TMM")
  v <- voom(d,design,plot=FALSE,normalize.method="none")
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  res <- topTable(fit,coef=ncol(design),n=nrow(d),sort.by="none")
  res <- dplyr::rename(res,PValue="P.Value")
  return(res)
}

runWilcoxon <- function(d){
  lsize <- colSums(d$counts)
  datas <- as.data.frame(log2(t(t(d$counts)/lsize*1e6)+1))
  group=d$samples[1]
  p<- sapply(1:nrow(datas),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(datas[i,])),group=group$group)
    p <- wilcox.test(gene~group, data)$p.value
    return(p)
    })
  res <- data.frame(row.names=rownames(datas), PValue=p)
  res$PValue <- as.numeric(res$PValue)
  return(res)
}

runTtest <- function(d){
  lsize <- colSums(d$counts)
  datas <- as.data.frame(log2(t(t(d$counts)/lsize*1e6)+1))
  group <- d$samples[1]
  p_value <- sapply(1:nrow(datas),function(i){
    data <- data.frame(gene=as.numeric(t(datas[i,])),group)
    b1 <- length(unique(split(data$gene,data$group)[[1]]))==1
    b2 <- length(unique(split(data$gene,data$group)[[2]]))==1
    p <- ifelse(all(b1,b2),1,t.test(data$gene~data$group)$p.value)
    return(p)})
  res <- data.frame(PValue=p_value,row.names=rownames(datas))
  return(res)  
}

runNBPSeq.ex <- function(d){
  cons <- unique(d$samples$group)
  norm.factors <- estimate.norm.factors(d$counts)
  obj <- prepare.nbp(d$counts, 
                    d$samples$group, 
                    lib.size=d$samples$lib.size, 
                    norm.factors=norm.factors)
  obj <- estimate.disp(obj)
  res <- exact.nb.test(obj, cons[1], cons[2], print.level=0)
  res <- data.frame(lfc=res$log.fc,
                    PValue=ifelse(is.na(res$p.values),1,res$p.values),
                    QValue=ifelse(is.na(res$q.values),1,res$q.values))

  return(res) 
}

runNOISeq <- function(d){
  cds <- newCountDataSet(d$counts,d$samples$group)
  res <- NOISeq::noiseqbio(cds,norm="tmm",factor = "condition")
  res <- res@results[[1]]
  res <- dplyr::mutate(res,PValue=1-prob) #not true pvalue
  return(res)
}

runABSSeq <- function(d){
  abs <- ABSDataSet(counts =d$counts ,groups = d$samples$group)
  abs <- ABSSeq(abs)
  res <- ABSSeq::results(abs,c("Amean","Bmean","foldChange","pvalue","adj.pvalue"))
  res <- data.frame(res)
  res <- dplyr::rename(res,PValue="pvalue")
  return(res)  
}

runROTS <- function(d){
  group <- as.numeric(d$sample$group)
  results <- ROTS(data = d$counts, groups = group , B = 1000 , seed = 1234)
  res <- data.frame(results$data,PValue=results$pvalue)
  return(res)
}

runSAMSeq <- function(d){
  set.seed(1)
  group <- as.numeric(d$sample$group)
  samfit <- SAMseq(d$counts, group, resp.type = "Two class unpaired")
  pval <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  res=data.frame(d$counts,PValue=pval)  
  return(res)
}

algos <- list("DESeq2"=runDESeq2,
              "DESeq"=runDESeq,
              "edgeR.lrt"=runEdgeR.lrt,
              "edgeR.qlf"=runEdgeR.qlf,
              "DSS"=runDSS,
              "voom"=runVoom,
              "NBPSeq"=runNBPSeq.ex,
              "Wilcoxon"=runWilcoxon,
              "T.test"=runTtest,
              "NOISeq"=runNOISeq,
              "ROTS"=runROTS,
              # "SAMSeq"=runSAMSeq,
              "ABSSeq"=runABSSeq)
