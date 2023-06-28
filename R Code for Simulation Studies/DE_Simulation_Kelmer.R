######################################################
#SIEVE: One-stop differential expression, variability,
# and skewness using RNA-Seq data
#Authors: Hongxiang Li and Tsung Fei Khang
#Email: chelsea.divo@hotmail.com
#Last update: 22 May 2023
#R Codes for comparison of methods
#Part 3: DE simulation study using the Kelmer dataset
#        (subgroup: at age of 10 weeks)
######################################################

## Required R packages
# install.packages("compositions")
# install.packages("devtools")
# install.packages("BiocManager")
# BiocManager::install("edgeR", force = T)
# BiocManager::install("DESeq2")
# BiocManager::install("limma", force = T)
# BiocManager::install("ALDEx2")
# BiocManager::install("DSS")
# BiocManager::install("tweeDEseq")
# BiocManager::install("dearseq")
# BiocManager::install("NOISeq")
# BiocManager::install("polyester")
library(polyester); library(compositions)
library(edgeR); library(DESeq2); library(limma)
library(ALDEx2); library(tweeDEseq); library(DSS)
library(NOISeq); library(dearseq); library(SIEVE)

################################
### GSE150318_counts.csv.gz
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150318
### Kelmer Dataset
################################
Data3 <- read.csv('...GSE150318_counts.csv.gz', header = T, check.names = TRUE, row.names = 1)
dim(Data3)

data_10weeks <- Data3[, seq(1,228,2)]
data_20weeks <- Data3[, seq(2,228,2)]
Data3 <- cbind(data_10weeks, data_20weeks)

# filter
CPM <- cpm(Data3)
keep <- (rowMeans(CPM[, 1:114]) > 0.5 & rowMeans(CPM[, 115:228]) > 0.5 & apply(Data3[, 1:114], 1, function(k) length(k[k==0])/length(k)) < 0.85 & apply(Data3[, 115:228], 1, function(k) length(k[k==0])/length(k)) < 0.85)
Data3 <- Data3[keep, ]
data_10weeks <- Data3[, 1:114]
dim(Data3); dim(data_10weeks)

## clr-transformation
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}

#############################
## Comparisons of Methods ###
#############################
DE_Comparison_simu <- function(data = NULL,
                               N.genes = 2000,
                               N.samples = NULL,
                               p.de = 0.1,
                               N.simulations = 30,
                               seed = NULL){

  set.seed(seed)
  seeds <- sample(1:10000, size = 2*N.simulations)
  seeds1 <- seeds[1:N.simulations]
  seeds2 <- seeds[(1 + N.simulations):(2*N.simulations)]

  group = rep(c(0,1), each=N.samples)
  mod = model.matrix(~-1 + group)
  # log_FC range between 2 DE genes
  log_FC <- c(seq(log(1/4), log(1/2), 0.000001), seq(log(2), log(4), 0.000001))
  N.de <- round(p.de*N.genes)  # number of DE genes
  params <- get_params(data)
  n = N.samples  # Number of samples for each group

  names <- c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq",
             "ALDEx2", "NOISeq", "DSS", "dearseq", "Wilcoxon")

  FDP <- matrix(NA, nrow = length(names), ncol = N.simulations) # false discovery proportion
  Type_II_error <- matrix(NA, nrow = length(names), ncol = N.simulations)
  Time <- matrix(NA, nrow = length(names), ncol = N.simulations)
  row.names(FDP) <- row.names(Type_II_error)<- row.names(Time)<- names
  colnames(FDP) <- colnames(Type_II_error) <- colnames(Time) <- paste("simu", 1:N.simulations, sep = "")

  i <- 1
  while (i <= N.simulations) {
    set.seed(seeds1[i])
    coeffs = cbind(c(sample(log_FC, size = N.de), rep(0, (N.genes - N.de))))

    data0 = create_read_numbers(params$mu, params$fit, params$p0,
                                m=dim(data)[1], n=dim(data)[2],
                                beta=coeffs, mod=mod, seed=seeds2[i])
    row.names(data0) <- paste('gene', 1:N.genes, sep='')
    de_gene <- row.names(data0)[1:N.de]

    CPM <- cpm(data0)
    keep <- (rowMeans(CPM[,1:n]) > 0.5 & rowMeans(CPM[,(n+1):(2*n)]) > 0.5 & apply(data0[,1:n], 1, function(x) length(x[x == 0])/length(x)) < 0.85 & apply(data0[,(n+1):(2*n)], 1, function(x) length(x[x == 0])/length(x)) < 0.85 )

    N.genes_de_filter <- sum(keep)   # number of genes after filtering
    N.de_filter <- sum(keep[1:N.de]) # number of DV genes after filtering
    de_gene <- de_gene[keep[1:N.de]] # names of DV genes after filtering
    data1 <- data0[keep, ]


    ### Comparisons of Methods
    ########
    ## SIEVE
    t1 <- proc.time()
    clr.counts1 <- clr.transform(data = data1)

    mu_test <- clrDE(data = clr.counts1, group = group)
    mu_test <- na.omit(mu_test)
    de_table <- mu_test[mu_test$adj_pval < 0.05, ]
    de_gene_clrDE <- row.names(de_table) # DE genes called

    FDP_SIEVE <- (length(de_gene_clrDE) -  length(intersect(de_gene_clrDE, de_gene)))/length(de_gene_clrDE)
    type_II_error_SIEVE <- (N.de_filter - length(intersect(de_gene_clrDE, de_gene)))/N.de_filter

    FDP[1, i] <- FDP_SIEVE
    Type_II_error[1, i] <- type_II_error_SIEVE
    Time[1, i] <- as.numeric(proc.time() - t1)[3]


    ########
    ## edgeR
    design = model.matrix(~group)

    t2 <- proc.time()
    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM") # normalization factors 

    dat.edgeR <- DGEList(counts=data1, norm.factors=nf, group=group)
    dat.edgeR <- estimateDisp(dat.edgeR, design)
    fit.edgeR <- glmQLFit(dat.edgeR, design)
    test.edgeR <- glmQLFTest(fit.edgeR, contrast = c(0,1))                                                                                                                                                       
    tp <- topTags(test.edgeR, n = Inf)
    top_table <- na.omit(tp$table)

    de_table <- top_table[top_table$FDR < 0.05, ]
    de_gene_edgeR <- row.names(de_table) # DE genes called
    # length(intersect(de_gene, de_gene_edgeR)) # True Positive

    FDP_edgeR = (length(de_gene_edgeR)  - length(intersect(de_gene, de_gene_edgeR)))/length(de_gene_edgeR)
    type_II_error_edgeR <- (N.de_filter - length(intersect(de_gene, de_gene_edgeR)))/N.de_filter

    FDP[2, i] <- FDP_edgeR
    Type_II_error[2, i] <- type_II_error_edgeR
    Time[2, i] <- as.numeric(proc.time() - t2)[3]


    #########
    ## DESeq2
    t3 <- proc.time()

    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM") # normalisation factors for edgeR, limma
    els <- nf*libsizes # effective library sizes for baySeq
    sf <- els/exp(mean(log(libsizes))) # size factors for DESeq2

    dat.DESeq2 <- DESeqDataSetFromMatrix(countData = data1, colData = data.frame(group),
                                         design = ~group)
    sizeFactors(dat.DESeq2) <- sf
    fit.DESeq2 <- DESeq(dat.DESeq2, minReplicatesForReplace=Inf)
    res.DESeq2 <- results(fit.DESeq2, cooksCutoff=F, alpha=0.05)
    res.DESeq2 <- na.omit(res.DESeq2)

    de_table <- res.DESeq2[res.DESeq2$padj < 0.05, ]
    de_gene_DESeq2 <- row.names(de_table)

    FDP_DESeq2 <- (length(de_gene_DESeq2) - length(intersect(de_gene_DESeq2, de_gene)))/length(de_gene_DESeq2)
    type_II_error_DESeq2 <- (N.de_filter -  length(intersect(de_gene_DESeq2, de_gene)))/N.de_filter

    FDP[3, i] <- FDP_DESeq2
    Type_II_error[3, i] <- type_II_error_DESeq2
    Time[3, i] <- as.numeric(proc.time() - t3)[3]


    #######
    ## voom
    t4 <- proc.time()

    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM")

    dat.edgeR <- DGEList(counts=data1, norm.factors=nf, group=group)
    dat.edgeR <- estimateDisp(dat.edgeR, design)
    dat.voom <- voom(dat.edgeR) # uses same data object as edgeR
    fit.voom <- lmFit(dat.voom, design)
    res.voom <- eBayes(fit.voom)

    top.table <- topTable(res.voom, sort.by = "P", n = Inf)
    top.table <- na.omit(top.table)
    de_gene_voom <- row.names(top.table[which(top.table$adj.P.Val < 0.05),])

    FDP_voom <- (length(de_gene_voom) - length(intersect(de_gene, de_gene_voom)))/length(de_gene_voom)
    type_II_error_voom <- (N.de_filter - length(intersect(de_gene, de_gene_voom)))/N.de_filter

    FDP[4, i] <- FDP_voom
    Type_II_error[4, i] <- type_II_error_voom
    Time[4, i] <- as.numeric(proc.time() - t4)[3]


    ### tweeDEseq
    t5 <- proc.time()

    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM") # normalisation factors
    els <- nf * libsizes # effective library sizes 
    sf <- els / exp(mean(log(libsizes))) # size factors 
    norm.data1 <- t(t(data1) / sf)

    counts <- as.matrix(norm.data1)
    resPT <- tweeDE(counts, group = group)
    resPT <- na.omit(resPT)
    # sum(resPT$pval.adjust < 0.05) # number of DE called
    de_table <- resPT[resPT$pval.adjust < 0.05, ]
    de_gene_twee<- row.names(de_table)
    # length(intersect(de_gene, de_gene_twee)) # True Positive

    FDP_twee = (dim(de_table)[1]  - length(intersect(de_gene, de_gene_twee)))/dim(de_table)[1]
    type_II_error_twee <- (length(de_gene) - length(intersect(de_gene, de_gene_twee)))/length(de_gene)

    FDP[5, i] <- FDP_twee
    Type_II_error[5, i] <- type_II_error_twee
    Time[5,i] <- as.numeric(proc.time() - t5)[3]


    ########
    ## ALDEx2
    t6 <- proc.time()

    conds <- as.character(group)
    we.T_test <- aldex(data1, conds, mc.samples = 128, denom = "all",
                       test = "t", effect = TRUE, paired.test = TRUE)
    we.T_test <- na.omit(we.T_test)
    de_table_ALDEx2 <- we.T_test[we.T_test$we.eBH < 0.05, ]
    de_genes_ALDEx2 <- row.names(de_table_ALDEx2)

    FDP_ALDEx2 <- (length(de_genes_ALDEx2) - length(intersect(de_genes_ALDEx2, de_gene)))/length(de_genes_ALDEx2)
    type_II_error_ALDEx2 <- (N.de_filter - length(intersect(de_genes_ALDEx2, de_gene)))/N.de_filter

    FDP[6, i] <- FDP_ALDEx2
    Type_II_error[6, i] <- type_II_error_ALDEx2
    Time[6, i] <- as.numeric(proc.time() - t6)[3]


    ##########
    ### NOISeq
    t7 <- proc.time()
    noiseq_out <- noiseqbio(readData(data = data1,
                                     factors = as.data.frame(group)),
                            k = 0.5, norm = "tmm", factor = "group",
                            random.seed = 12345, filter = 0)
    # data1 was filtered already,
    # all methods use the same filtering criteria
    noiseq_result <- noiseq_out@results[[1]]
    noiseq_result <- na.omit(noiseq_result)
    de_table_noiseq <- noiseq_result[noiseq_result$prob > 0.95, ]
    # prob = 1 - fdr
    de_genes_noiseq <- row.names(de_table_noiseq)

    FDP_noiseq <- (length(de_genes_noiseq) - length(intersect(de_genes_noiseq, de_gene)))/length(de_genes_noiseq)
    type_II_error_noiseq <- (N.de_filter - length(intersect(de_genes_noiseq, de_gene)))/N.de_filter

    FDP[7, i] <- FDP_noiseq
    Type_II_error[7, i] <- type_II_error_noiseq
    Time[7, i] <- as.numeric(proc.time() - t7)[3]


    #######
    ### DSS
    t8 <- proc.time()
    #    colnames(data1) <- NULL
    #    seqData <- newSeqCountSet(data1, group)
    #    seqData <- estNormFactors(seqData, method = "lr") # similar to TMM. NO TMM choice in DSS package
    #    seqData <- estDispersion(seqData)

    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM") # normalization factors
    els <- nf*libsizes # effective library sizes
    sf <- els/exp(mean(log(libsizes))) # size factors

    seqData <- newSeqCountSet(counts = data1, designs = group,
                              normalizationFactor = sf)
    seqData <- estDispersion(seqData)

    result_DSS <- waldTest(seqData, 0, 1)
    result_DSS <- na.omit(result_DSS)

    de_table_DSS <- result_DSS[result_DSS$fdr < 0.05, ]
    de_genes_DSS <- row.names(de_table_DSS)

    FDP_DSS <- (length(de_genes_DSS) - length(intersect(de_genes_DSS, de_gene)))/length(de_genes_DSS)
    type_II_error_DSS <- (N.de_filter - length(intersect(de_genes_DSS, de_gene)))/N.de_filter

    FDP[8, i] <- FDP_DSS
    Type_II_error[8, i] <- type_II_error_DSS
    Time[8, i] <- as.numeric(proc.time() - t8)[3]


    ## dearseq
    t9 <- proc.time()
    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM") # normalisation factors 
    els <- nf * libsizes # effective library sizes 
    sf <- els / exp(mean(log(libsizes))) # size factors 
    norm.data1 <- t(t(data1) / sf) # normalised count matrix 

    conditions <- matrix(group, ncol=1)
    res <- dear_seq(exprmat = norm.data1,
                    variables2test = conditions,
                    which_test = "asymptotic",
                    parallel_comp = F,
                    preprocessed = T)

    res_pvals <- res$pvals
    # sum(res_pvals$adjPval < 0.05)

    de_table_dearseq <- res_pvals[res_pvals$adjPval < 0.05, ]
    de_gene_dearseq <- row.names(de_table_dearseq) # DE genes called

    FDP_dearseq = (length(de_gene_dearseq)  - length(intersect(de_gene, de_gene_dearseq)))/length(de_gene_dearseq)
    type_II_error_dearseq <- (N.de_filter - length(intersect(de_gene, de_gene_dearseq)))/N.de_filter

    FDP[9, i] <- FDP_dearseq
    Type_II_error[9, i] <- type_II_error_dearseq
    Time[9, i] <- as.numeric(proc.time() - t9)[3]



    ## Wilcoxon rank-sum test
    t10 <- proc.time()
    libsizes <- colSums(data1)
    nf <- calcNormFactors(data1, method="TMM") # normalisation factors 
    els <- nf * libsizes # effective library sizes 
    sf <- els / exp(mean(log(libsizes))) # size factors 
    norm.data1 <- t(t(data1) / sf) # normalised count matrix 

    pvalues <- sapply(1:nrow(norm.data1), function(i){
      data<-cbind.data.frame(gene=as.numeric(t(norm.data1[i,])), group)
      p=wilcox.test(gene~group, data)$p.value
      return(p)
    })
    adj_pval <- p.adjust(pvalues, method = "BH")

    result_Wilcoxon <- data.frame(pval=pvalues, adj_pval=adj_pval)
    rownames(result_Wilcoxon)=rownames(data1)

    de_table_Wilcoxon <- result_Wilcoxon[result_Wilcoxon$adj_pval < 0.05, ]
    de_gene_Wilcoxon <- row.names(de_table_Wilcoxon)

    FDP_Wilcoxon <- (length(de_gene_Wilcoxon) -  length(intersect(de_gene_Wilcoxon, de_gene)))/length(de_gene_Wilcoxon)
    type_II_error_Wilcoxon <- (N.de_filter - length(intersect(de_gene_Wilcoxon, de_gene)))/N.de_filter

    FDP[10, i] <- FDP_Wilcoxon
    Type_II_error[10, i] <- type_II_error_Wilcoxon
    Time[10, i] <- as.numeric(proc.time() - t10)[3]
    #

    i <- i+1
  }
  return(list("FDP" =  FDP, "Type_II_error" = Type_II_error, "Time" = Time))
}


### Methods Comparison
## Comparison One, 2*30 samples
DE_test1 <- DE_Comparison_simu(data = data_10weeks, N.genes = 2000,
                               N.samples = 30, p.de = 0.1,
                               N.simulations = 30, seed = 3003)
DE_test1$FDP; DE_test1$Type_II_error; DE_test1$Time

write.csv(round(DE_test1$FDP, 4), "DE_FDP_10weeks_30samples.csv")
write.csv(round(DE_test1$Type_II_error, 4), "DE_Type_II_error_10weeks_30samples.csv")
write.csv(round(DE_test1$Time, 2), "DE_Time_10weeks_30samples.csv")


## Comparison Two, 2*50 samples
DE_test2 <- DE_Comparison_simu(data = data_10weeks, N.genes = 2000,
                               N.samples = 50, p.de = 0.1,
                               N.simulations = 30, seed = 4003)
DE_test2$FDP; DE_test2$Type_II_error; DE_test2$Time

write.csv(round(DE_test2$FDP, 4), "DE_FDP_10weeks_50samples.csv")
write.csv(round(DE_test2$Type_II_error, 4), "DE_Type_II_error_10weeks_50samples.csv")
write.csv(round(DE_test2$Time, 2), "DE_Time_10weeks_50samples.csv")


## Comparison Three, 2*100 samples
DE_test3 <- DE_Comparison_simu(data = data_10weeks, N.genes = 2000,
                               N.samples = 100, p.de = 0.1,
                               N.simulations = 30, seed = 6003)
DE_test3$FDP; DE_test3$Type_II_error; DE_test3$Time

write.csv(round(DE_test3$FDP, 4), "DE_FDP_10weeks_100samples.csv")
write.csv(round(DE_test3$Type_II_error, 4), "DE_Type_II_error_10weeks_100samples.csv")
write.csv(round(DE_test3$Time, 2), "DE_Time_10weeks_100samples.csv")


################
### Analysis ###
################
### Scater plot, Fig. 3
setwd("~/.../Kelmer_10_weeks")
FDR.Kelmer.30samples <- read.csv('DE_FDP_10weeks_30samples.csv', row.names = 1)
FDR.Kelmer.50samples <- read.csv('DE_FDP_10weeks_50samples.csv', row.names = 1)
FDR.Kelmer.100samples <- read.csv('DE_FDP_10weeks_100samples.csv', row.names = 1)

type2.Kelmer.30samples <- read.csv('DE_Type_II_error_10weeks_30samples.csv', row.names = 1)
type2.Kelmer.50samples <- read.csv('DE_Type_II_error_10weeks_50samples.csv', row.names = 1)
type2.Kelmer.100samples <- read.csv('DE_Type_II_error_10weeks_100samples.csv', row.names = 1)

### Jitter Plots
### 2*30 samples, Fig.3 (a)
plot(jitter(as.numeric(FDR.Kelmer.30samples[4, ]), amount = 0.002),
     jitter(as.numeric(type2.Kelmer.30samples[4, ]), amount = 0.002),
     lwd = 1.5, cex = 1, col = "#1b9e77", pch = 43,
     xlim = c(0, 0.125), ylim = c(0, 0.1),
     xlab = "FDR", ylab = "Type II Error", cex.lab = 1.15)

points(jitter(as.numeric(FDR.Kelmer.30samples[5, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[5, ]), amount = 0.002),
       lwd =1.5, cex=1, col = "#7570b3", pch =15)

points(jitter(as.numeric(FDR.Kelmer.30samples[2, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[2, ]), amount = 0.002),
       col = "#1b9e77", lwd =1.5, cex=1, pch = 2 )

points(jitter(as.numeric(FDR.Kelmer.30samples[3, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[3, ]), amount = 0.002),
       lwd =1.5, cex=0.75, col = "black", pch =8)

points(jitter(as.numeric(FDR.Kelmer.30samples[1, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[1, ]), amount = 0.002),
       lwd =1.5, cex=1.25, pch = 1, col = "#d95f02")

points(jitter(as.numeric(FDR.Kelmer.30samples[6, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[6, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 92)

points(jitter(as.numeric(FDR.Kelmer.30samples[7, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[7, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 95)

points(jitter(as.numeric(FDR.Kelmer.30samples[8, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[8, ]), amount = 0.002),
       lwd = 1, cex = 1, col = "#d95f02", pch = 6)

points(jitter(as.numeric(FDR.Kelmer.30samples[9, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[9, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "#1b9e77", pch = 1)

points(jitter(as.numeric(FDR.Kelmer.30samples[9, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.30samples[9, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 47)

abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("topright", 
       c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq", 
         "ALDEx2", "NOISeq", "DSS", "dearseq", "Wilcoxon"),
       pch = c(1, 2, 8, 43, 15, 
               92, 95, 6, 1, 47), 
       pt.cex = c(1.25,1,1,1,1,1,1,1,1,1),
       cex= 1.1, pt.lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("#d95f02", "#1b9e77", "black", "#1b9e77", "#7570b3", 
               "black", "black", "#d95f02", "#1b9e77", "black"))
title(adj = 0, "(a)")


### Jitter Plots
### 2*50 samples, Fig.3 (b)
plot(jitter(as.numeric(FDR.Kelmer.50samples[4, ]), amount = 0.002),
     jitter(as.numeric(type2.Kelmer.50samples[4, ]), amount = 0.002),
     lwd = 1.5, cex = 1, col = "#1b9e77", pch = 43,
     xlim = c(0, 0.1), ylim = c(-0.001, 0.035),
     xlab = "FDR", ylab = "Type II Error", cex.lab = 1.15)

points(jitter(as.numeric(FDR.Kelmer.50samples[5, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[5, ]), amount = 0.002),
       lwd =1.5, cex=1, col = "#7570b3", pch =15)

points(jitter(as.numeric(FDR.Kelmer.50samples[2, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[2, ]), amount = 0.002),
       col = "#1b9e77", lwd =1.5, cex=1, pch = 2 )

points(jitter(as.numeric(FDR.Kelmer.50samples[3, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[3, ]), amount = 0.002),
       lwd =1.5, cex=0.75, col = "black", pch =8)

points(jitter(as.numeric(FDR.Kelmer.50samples[1, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[1, ]), amount = 0.002),
       lwd =1.5, cex=1.25, pch = 1, col = "#d95f02")

points(jitter(as.numeric(FDR.Kelmer.50samples[6, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[6, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 92)

points(jitter(as.numeric(FDR.Kelmer.50samples[7, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[7, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 95)

points(jitter(as.numeric(FDR.Kelmer.50samples[8, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[8, ]), amount = 0.002),
       lwd = 1, cex = 1, col = "#d95f02", pch = 6)

points(jitter(as.numeric(FDR.Kelmer.50samples[9, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[9, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "#1b9e77", pch = 1)

points(jitter(as.numeric(FDR.Kelmer.50samples[9, ]), amount = 0.002),
       jitter(as.numeric(type2.Kelmer.50samples[9, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 47)

abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("topright", 
       c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq", 
         "ALDEx2", "NOISeq", "DSS", "dearseq", "Wilcoxon"),
       pch = c(1, 2, 8, 43, 15, 
               92, 95, 6, 1, 47), 
       pt.cex = c(1.25,1,1,1,1,1,1,1,1,1),
       cex= 1.1, pt.lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("#d95f02", "#1b9e77", "black", "#1b9e77", "#7570b3", 
               "black", "black", "#d95f02", "#1b9e77", "black"))
title(adj = 0, "(b)")


### Jitter Plots
### 2*100 samples, Fig.3 (c)
plot(jitter(as.numeric(FDR.Kelmer.100samples[4, ]), amount = 0.0005),
     jitter(as.numeric(type2.Kelmer.100samples[4, ]), amount = 0.0005),
     lwd = 1.5, cex = 1, col = "#1b9e77", pch = 43,
     xlim = c(-0.0002, 0.105), ylim = c(-0.0002, 0.0105),
     xlab = "FDR", ylab = "Type II Error", cex.lab = 1.15)

points(jitter(as.numeric(FDR.Kelmer.100samples[5, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[5, ]), amount = 0.0005),
       lwd =1.5, cex=1, col = "#7570b3", pch =15)

points(jitter(as.numeric(FDR.Kelmer.100samples[2, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[2, ]), amount = 0.0005),
       col = "#1b9e77", lwd =1.5, cex=1, pch = 2 )

points(jitter(as.numeric(FDR.Kelmer.100samples[3, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[3, ]), amount = 0.0005),
       lwd =1.5, cex=0.75, col = "black", pch =8)

points(jitter(as.numeric(FDR.Kelmer.100samples[1, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[1, ]), amount = 0.0005),
       lwd =1.5, cex=1.25, pch = 1, col = "#d95f02")

points(jitter(as.numeric(FDR.Kelmer.100samples[6, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[6, ]), amount = 0.0005),
       lwd = 1.5, cex = 1, col = "black", pch = 92)

points(jitter(as.numeric(FDR.Kelmer.100samples[7, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[7, ]), amount = 0.0005),
       lwd = 1.5, cex = 1, col = "black", pch = 95)

points(jitter(as.numeric(FDR.Kelmer.100samples[8, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[8, ]), amount = 0.0005),
       lwd = 1, cex = 1, col = "#d95f02", pch = 6)

points(jitter(as.numeric(FDR.Kelmer.100samples[9, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[9, ]), amount = 0.0005),
       lwd = 1.5, cex = 1, col = "#1b9e77", pch = 1)

points(jitter(as.numeric(FDR.Kelmer.100samples[9, ]), amount = 0.0005),
       jitter(as.numeric(type2.Kelmer.100samples[9, ]), amount = 0.0005),
       lwd = 1.5, cex = 1, col = "black", pch = 47)

abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("topright", 
       c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq", 
         "ALDEx2", "NOISeq", "DSS", "dearseq", "Wilcoxon"),
       pch = c(1, 2, 8, 43, 15, 
               92, 95, 6, 1, 47), 
       pt.cex = c(1.25,1,1,1,1,1,1,1,1,1),
       cex= 1.1, pt.lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("#d95f02", "#1b9e77", "black", "#1b9e77", "#7570b3", 
               "black", "black", "#d95f02", "#1b9e77", "black"))
title(adj = 0, "(c)")


######################
### Heatmap, Fig.5 (b)
stat_fdr <- vector("list",3)
stat_t2e <- vector("list",3)
stat_time <- vector("list",3)

names(stat_fdr) <- c(30,50,100)
setwd("~/.../Kelmer_10_weeks")
dir() # be careful with the order of csv files in the folder

id <- c(2,3,1) # FDP, (30, 50, 100)
for(i in 1:3){
  dat <- read.csv(file=(dir())[id[i]], row.names=1)
  dat <- t(dat)
  stat_fdr[[i]] <- round(t(rbind(apply(dat, 2, mean), apply(dat, 2, sd))),4)
  colnames(stat_fdr[[i]]) <- c("Mean","SD")
}

id <- c(9,10,8) # Type II error, (30,m 50, 100)
for(i in 1:3){
  dat <- read.csv(file=(dir())[id[i]], row.names=1)
  dat <- t(dat)
  stat_t2e[[i]] <- round(t(rbind(apply(dat, 2, mean), apply(dat, 2, sd))),4)
  colnames(stat_t2e[[i]]) <- c("Mean","SD")
}

id <- c(6,7,5) # computing time, (30, 50, 100)
for(i in 1:3){
  dat <- read.csv(file=(dir())[id[i]], row.names=1)
  dat <- t(dat)
  stat_time[[i]] <- round(t(rbind(apply(dat, 2, mean), apply(dat, 2, sd))),3)
  colnames(stat_time[[i]]) <- c("Mean","SD")
}

fdr_table <- sd_fdr_table <- mean_t2e_table <- sd_t2e_table <- mean_time_table <- sd_time_table <- matrix(0,10,3)
rownames(fdr_table) <- rownames(sd_fdr_table) <- rownames(mean_t2e_table) <- rownames(sd_t2e_table) <- rownames(mean_time_table) <- rownames(sd_time_table) <- row.names(stat_fdr[[1]])
colnames(fdr_table) <- colnames(sd_fdr_table) <- colnames(mean_t2e_table) <- colnames(sd_t2e_table) <- colnames(mean_time_table) <- colnames(sd_time_table) <- c("n=30", "n=50", "n=100")

for(i in 1:10){
  fdr_table[i,] <- unlist(lapply(stat_fdr, function(k) k[i,1]))
  sd_fdr_table[i,] <- unlist(lapply(stat_fdr, function(k) k[i,2]))
  mean_t2e_table[i,] <- unlist(lapply(stat_t2e, function(k) k[i,1]))
  sd_t2e_table[i,] <- unlist(lapply(stat_t2e, function(k) k[i,2]))
  mean_time_table[i,] <- unlist(lapply(stat_time, function(k) k[i,1]))
  sd_time_table[i,] <- unlist(lapply(stat_time, function(k) k[i,2]))
}

mean_time_table[c(1,5,6), ] <- round(mean_time_table[c(1,5,6), ],0)
mean_time_table[-c(1,5,6), ] <- round(mean_time_table[-c(1,5,6), ], 1)
sd_time_table[c(1,5), ] <- round(sd_time_table[c(1,5), ], 0)
sd_time_table[-c(1,5), ] <- round(sd_time_table[-c(1,5), ], 2)

# Table S1 (superscript *)
# write.csv(round(fdr_table, 3), file="fdr_table_Kelmer.csv", quote=FALSE, row.names=T)
# write.csv(round(sd_fdr_table, 3), file="sd_fdr_table_Kelmer.csv", quote=FALSE, row.names=T)

# Table S1 (superscript \dag)
# write.csv(round(mean_t2e_table, 3), file="mean_t2e_table_Kelmer.csv", quote=FALSE, row.names=T)
# write.csv(round(sd_t2e_table, 3), file="sd_t2e_table_Kelmer.csv", quote=FALSE, row.names=T)

# Table S3 (superscript \dag)
# write.csv(mean_time_table, file="mean_time_table_Kelmer.csv", quote=FALSE, row.names=T)
# write.csv(sd_time_table, file="sd_time_table_Kelmer.csv", quote=FALSE, row.names=T)

###
#Sometimes, heatmap clustering can be used to study similarity of methods
X <- data.frame(fdr_table, sd_fdr_table, mean_t2e_table, sd_t2e_table)
colnames(X) <- c("fdr_m30","fdr_m50","fdr_m100",
                 "fdr_sd30","fdr_sd50","fdr_sd100",
                 "t2e_m30","t2e_m50","t2e_m100",
                 "t2e_sd30","t2e_sd50","t2e_sd100")

#install.packages(gplots)
library(gplots)

#set color tones
greensp <- c("#ffffcc", "#d9f0a3", "#addd8e","#78c679", "#31a354")
greensp2 <- colorRampPalette(greensp, space = "Lab")(5)

X_norm <- apply(X,2,function(k) (k-min(k))/(max(k) - min(k)))
pca <- prcomp(X_norm)

#heatmap, Fig.5 (b)
heatmap.2(t(X_norm),
          cellnote = round(t(X_norm), 2),
          notecex = 0.7, notecol = "black",
          reorderfun = function(d, w) reorder(d, w, agglo.FUN = min),
          hclustfun = function(k) hclust(dist(k), method="ward.D"),
          dist=function(k)dist(k,method="euclidean"),Rowv=TRUE,
          symkey=FALSE, density.info = "none",trace="none",Colv=TRUE,
          col=greensp2, labCol=rownames(pca$x),dendrogram="both",
          xlab="Methods",srtRow=30, cexRow=1, srtCol=30, cexCol=1)


################################
### Percentage of FDR < 0.05 and probability of Type II error < 0.05
# marginal prob.
#data.frame(cbind(
#  apply(FDR.Kelmer.30samples, 1, function(k)  (sum(k<0.05)/30)*100),
#  apply(FDR.Kelmer.50samples, 1, function(k)  (sum(k<0.05)/30)*100),
#  apply(FDR.Kelmer.100samples, 1, function(k)  (sum(k<0.05)/30)*100)))

#data.frame(cbind(
#  apply(type2.Kelmer.30samples, 1, function(k)  (sum(k<0.05)/30)*100),
#  apply(type2.Kelmer.50samples, 1, function(k)  (sum(k<0.05)/30)*100),
#  apply(type2.Kelmer.100samples, 1, function(k)  (sum(k<0.05)/30)*100)))

# joint prob.  Table S1 (superscript, \ddag)
tab1 <- (FDR.Kelmer.30samples < 0.05) + (type2.Kelmer.30samples < 0.05)
tab2 <- (FDR.Kelmer.50samples < 0.05) + (type2.Kelmer.50samples < 0.05)
tab3 <- (FDR.Kelmer.100samples < 0.05) + (type2.Kelmer.100samples < 0.05)

apply(tab1, 1, function(k)  (sum(k == 2)/30)*100)
apply(tab2, 1, function(k)  (sum(k == 2)/30)*100)
apply(tab3, 1, function(k)  (sum(k == 2)/30)*100)
###END###
