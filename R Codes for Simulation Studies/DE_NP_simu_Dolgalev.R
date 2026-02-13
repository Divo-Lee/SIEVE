######################################################
#SIEVE: One-stop differential expression, variability,
# and skewness using RNA-Seq data
#Authors: Hongxiang Li and Tsung Fei Khang
#Email: chelsea.divo@hotmail.com
#Last update: 23 Nov. 2025
#R Codes for comparison of methods
#Part 4: DE nonPara. simulation study using the Koebbe dataset
#        (GSE242339)
######################################################

## Required R packages
# install.packages("compositions")
# install.packages("devtools")
# install.packages("BiocManager")
# install.packages("SimSeq")
# install.packages("gplots")
# BiocManager::install("edgeR", force = T)
# BiocManager::install("DESeq2")
# BiocManager::install("limma", force = T)
# BiocManager::install("ALDEx2")
# BiocManager::install("DSS")
# BiocManager::install("tweeDEseq")
# BiocManager::install("dearseq")
# BiocManager::install("NOISeq")
# BiocManager::install("polyester")


library(edgeR); library(DESeq2); library(limma)
library(ALDEx2); library(tweeDEseq); library(DSS)
library(NOISeq); library(dearseq); library(SIEVE)
library(compositions); library(SimSeq); library(gplots)


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
DE_Comparison_simseq_simu <- function(data = NULL,
                                      groups.0 = NULL,
                                      N.genes = 2000,
                                      N.samples = NULL,
                                      N.de = 200,
                                      lfc.de = 1,
                                      N.simulations = 30,
                                      seed = NULL){

  set.seed(seed)
  seeds1 <- sample(1:10000, size = N.simulations)


  ### simseq pre-process
  ####
  groups.0 <- as.factor(groups.0)
  nf0 <- calcNormFactors(data, method="TMM")
  sort.list <- SortData(counts = data,
                        treatment = groups.0,
                        sort.method = "unpaired",
                        norm.factors = nf0)

  counts.0 <- sort.list$counts
  treatment <- sort.list$treatment
  nf1 <- sort.list$norm.factors

  probs <- CalcPvalWilcox(counts= data,
                          treatment = groups.0,
                          sort.method = "unpaired",
                          sorted = T,
                          norm.factors = nf1,
                          exact = FALSE)
  weights <- 1 - fdrtool(probs,
                         statistic = "pvalue",
                         plot = FALSE,
                         verbose = FALSE)$lfdr

  mean_group1 <- rowMeans(log2( (counts.0[, treatment==0] + 1) %*% diag(1/nf1[treatment==0])))
  mean_group2 <- rowMeans(log2( (counts.0[, treatment==1] + 1) %*% diag(1/nf1[treatment==1])))

  lfc <- mean_group1 - mean_group2
  # sum(abs(lfc) > 1) #
  N.de <- min(N.de, sum(abs(lfc) > lfc.de))


  weights.zero <- abs(lfc) <= lfc.de
  weights[weights.zero] <- 0
  ##


  n = N.samples  # Number of samples for each group

  names <- c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq",
             "ALDEx2", "DSS", "dearseq", "NOISeq")

  FDP <- matrix(NA, nrow = length(names), ncol = N.simulations) # false discovery proportion
  Type_II_error <- matrix(NA, nrow = length(names), ncol = N.simulations)
  Time <- matrix(NA, nrow = length(names), ncol = N.simulations)
  row.names(FDP) <- row.names(Type_II_error)<- row.names(Time)<- names
  colnames(FDP) <- colnames(Type_II_error) <- colnames(Time) <- paste("simu", 1:N.simulations, sep = "")


  i <- 1
  while (i <= N.simulations) {

    set.seed(seeds1[i])
    data.sim <- SimData(counts = counts.0,
                        treatment = treatment,
                        sort.method = "unpaired",
                        k.ind = N.samples,
                        n.genes = N.genes,
                        n.diff = N.de,
                        weights = weights,
                        norm.factors = nf1)

    sim_count1 <- data.sim$counts
    group <- data.sim$treatment

    # filtering
    keep <- (rowMeans(cpm(sim_count1[ , 1:n])) > 0.5  &
               rowMeans(cpm(sim_count1[ , (1+n):(2*n)])) > 0.5 &
               apply(sim_count1[ , 1:n], 1, function(x) length(x[x==0])/length(x)) < 0.85 &
               apply(sim_count1[ , (1+n):(2*n)], 1, function(x) length(x[x==0])/length(x)) < 0.85)
    data1 <- sim_count1 <- sim_count1[keep, ] # sim_count1

    spike.DE.genes <- rownames(data)[data.sim$DE.genes]
    de_gene <- intersect(spike.DE.genes, rownames(data1)) # names of DE genes after filtering
    N.de_filter <- length(de_gene) # number of DE genes after filtering



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

    dat.DESeq2 <- DESeqDataSetFromMatrix(countData = data1,
                                         colData = data.frame(group),
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



    #######
    ### DSS
    t7 <- proc.time()
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

    FDP[7, i] <- FDP_DSS
    Type_II_error[7, i] <- type_II_error_DSS
    Time[7, i] <- as.numeric(proc.time() - t7)[3]


    ## dearseq
    t8 <- proc.time()
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

    FDP[8, i] <- FDP_dearseq
    Type_II_error[8, i] <- type_II_error_dearseq
    Time[8, i] <- as.numeric(proc.time() - t8)[3]


    ##########
    ### NOISeq
    t9 <- proc.time()

    groups <- c(rep("A", n), rep("B", n))
    groups_df <- as.data.frame(groups)
    # 明确强制 rownames 等于样本名
    rownames(groups_df) <-  paste(groups, 1:(2*n),  sep='')
    colnames(data1) <-  rownames(groups_df)


    noiseq_out <- noiseqbio(readData(data = data1,
                                     factors = as.data.frame(groups_df)),
                            k = 0.5, norm = "tmm", factor = "groups",
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

    FDP[9, i] <- FDP_noiseq
    Type_II_error[9, i] <- type_II_error_noiseq
    Time[9, i] <- as.numeric(proc.time() - t9)[3]

    #

    i <- i+1
  }
  return(list("FDP" =  FDP, "Type_II_error" = Type_II_error, "Time" = Time))
}


############################
### Data 4
### GSE229705_counts-raw.csv.gz
### normal = 123, tumor = 123
############################
Data1 = read.csv(gzfile('GSE229705_counts-raw.csv.gz'), header = T, check.names = TRUE, row.names = 1)
dim(Data1) # 60591 246

# Filtering
CPM <- cpm(Data1)
keep <- (rowMeans(CPM[,1:123]) > 0.5 &
         rowMeans(CPM[,124:246]) > 0.5 &
         apply(Data1[,1:123], 1, function(x) length(x[x==0])/length(x)) < 0.85 &
         apply(Data1[,124:246], 1, function(x) length(x[x==0])/length(x)) < 0.85)
data1 <- Data1[keep, ]; dim(data1) # 20151

row.names(data1) <- paste('gene', 1:nrow(data1), sep='')
groups.0 <- c(rep(0, 123), rep(1, 123))


## Comparison One, 2*30 samples
DE_test1 <- DE_Comparison_simseq_simu(data = data1,
                                      groups.0 = groups.0,
                                      N.genes = 2000,
                                      N.samples = 30,
                                      N.de = 200,
                                      lfc.de = 1,
                                      N.simulations = 3,
                                      seed = 104)
DE_test1$FDP; DE_test1$Type_II_error; DE_test1$Time

rowMeans(DE_test1$FDP); rowSds(DE_test1$FDP)
rowMeans(DE_test1$Type_II_error); rowSds(DE_test1$Type_II_error)

write.csv(round(DE_test1$FDP, 5), "DE_FDP_Dolgalev_30samples.csv", quote=FALSE, row.names=T)
write.csv(round(DE_test1$Type_II_error, 4), "DE_Type_II_error_Dolgalev_30samples.csv", quote=FALSE, row.names=T)
write.csv(round(DE_test1$Time, 2), "DE_Time_Dolgalev_30samples.csv", quote=FALSE, row.names=T)


## Comparison Two, 2*50 samples
DE_test2 <- DE_Comparison_simseq_simu(data = data1,
                                      groups.0 = groups.0,
                                      N.genes = 2000,
                                      N.samples = 50,
                                      N.de = 200,
                                      lfc.de = 1,
                                      N.simulations = 3,
                                      seed = 105)
DE_test2$FDP; DE_test2$Type_II_error; DE_test2$Time

rowMeans(DE_test2$FDP); rowSds(DE_test2$FDP)
rowMeans(DE_test2$Type_II_error); rowSds(DE_test2$Type_II_error)

write.csv(round(DE_test2$FDP, 5), "DE_FDP_Dolgalev_50samples.csv", quote=FALSE, row.names=T)
write.csv(round(DE_test2$Type_II_error, 4), "DE_Type_II_error_Dolgalev_50samples.csv", quote=FALSE, row.names=T)
write.csv(round(DE_test2$Time, 2), "DE_Time_Dolgalev_50samples.csv", quote=FALSE, row.names=T)

################
### Analysis ###
################
stat_fdr <- vector("list",2)
stat_t2e <- vector("list",2)
stat_time <- vector("list",2)

names(stat_fdr) <- c(30,50)
# dir()

id <- c(1,2)
for(i in 1:2){
  dat <- read.csv(file=(dir())[id[i]], row.names=1)
  dat <- t(dat)
  stat_fdr[[i]] <- round(t(rbind(apply(dat, 2, mean), apply(dat, 2, sd))),4)
  colnames(stat_fdr[[i]]) <- c("Mean","SD")
}

id <- c(6,7)
for(i in 1:2){
  dat <- read.csv(file=(dir())[id[i]], row.names=1)
  dat <- t(dat)
  stat_t2e[[i]] <- round(t(rbind(apply(dat, 2, mean), apply(dat, 2, sd))),4)
  colnames(stat_t2e[[i]]) <- c("Mean","SD")
}

id <- c(4,5)
for(i in 1:2){
  dat <- read.csv(file=(dir())[id[i]], row.names=1)
  dat <- t(dat)
  stat_time[[i]] <- round(t(rbind(apply(dat, 2, mean), apply(dat, 2, sd))),3)
  colnames(stat_time[[i]]) <- c("Mean","SD")
}


fdr_table <- matrix(0,9,2)
sd_fdr_table <-matrix(0,9,2)
mean_t2e_table <- matrix(0,9,2)
sd_t2e_table <- matrix(0,9,2)
mean_time_table <- matrix(0,9,2)
sd_time_table <- matrix(0,9,2)

rownames(fdr_table) <- row.names(stat_fdr[[1]])
rownames(sd_fdr_table) <- row.names(stat_fdr[[1]])
rownames(mean_t2e_table) <- row.names(stat_fdr[[1]])
rownames(sd_t2e_table) <- row.names(stat_fdr[[1]])
rownames(mean_time_table) <- row.names(stat_fdr[[1]])
rownames(sd_time_table) <- row.names(stat_fdr[[1]])
colnames(fdr_table) <- c("n=30", "n=50")
colnames(sd_fdr_table) <- c("n=30", "n=50")
colnames(mean_t2e_table) <-c("n=30", "n=50")
colnames(sd_t2e_table) <- c("n=30", "n=50")
colnames(mean_time_table) <- c("n=30", "n=50")
colnames(sd_time_table) <- c("n=30", "n=50")


for(i in 1:9){
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


X <- data.frame(fdr_table, sd_fdr_table, mean_t2e_table, sd_t2e_table)
colnames(X) <- c("fdr_m30","fdr_m50",
                 "fdr_sd30","fdr_sd50",
                 "t2e_m30","t2e_m50",
                 "t2e_sd30","t2e_sd50")


#set color tones
greensp <- c("#ffffcc", "#d9f0a3", "#addd8e","#78c679", "#31a354")
greensp2 <- colorRampPalette(greensp, space = "Lab")(5)

#X_norm <- apply(X,2,function(k) (k-mean(k))/sd(k) )
X_norm <- apply(X,2,function(k) (k-min(k))/(max(k) - min(k)))
pca <- prcomp(X_norm)

#heatmap
heatmap.2(t(X_norm),
          cellnote = round(t(X_norm), 2),
          notecex = 0.7, notecol = "black",
          reorderfun = function(d, w) reorder(d, w, agglo.FUN = min),
          hclustfun = function(k) hclust(dist(k), method="ward.D"),
          dist=function(k)dist(k,method="euclidean"),Rowv=TRUE,
          symkey=FALSE, density.info = "none",trace="none",Colv=TRUE,
          col=greensp2, labCol=rownames(pca$x),dendrogram="both",
          xlab="Methods",srtRow=30, cexRow=1, srtCol=30, cexCol=1)

###
### Jitter Plots
FDR.Dolgalev.30samples <- read.csv('DE_FDP_Dolgalev_30samples.csv', row.names = 1)
FDR.Dolgalev.50samples <- read.csv('DE_FDP_Dolgalev_50samples.csv', row.names = 1)
type2.Dolgalev.30samples <- read.csv('DE_Type_II_error_Dolgalev_30samples.csv', row.names = 1)
type2.Dolgalev.50samples <- read.csv('DE_Type_II_error_Dolgalev_50samples.csv', row.names = 1)


### 2*30 samples
plot(jitter(as.numeric(FDR.Dolgalev.30samples[5, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[5, ]), amount = 0.002),
       lwd =1.5, cex=1, col = "#7570b3", pch =15,
       xlim = c(0, 0.5), ylim = c(0, 0.9),
       xlab = "FDR", ylab = "Type II Error", cex.lab = 1.15)

points(jitter(as.numeric(FDR.Dolgalev.30samples[4, ]), amount = 0.002),
     jitter(as.numeric(type2.Dolgalev.30samples[4, ]), amount = 0.002),
     lwd = 1.5, cex = 1, col = "#1b9e77", pch = 43)

points(jitter(as.numeric(FDR.Dolgalev.30samples[2, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[2, ]), amount = 0.002),
       col = "#1b9e77", lwd =1.5, cex=1, pch = 2 )

points(jitter(as.numeric(FDR.Dolgalev.30samples[3, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[3, ]), amount = 0.002),
       lwd =1.5, cex=0.75, col = "black", pch =8)


points(jitter(as.numeric(FDR.Dolgalev.30samples[6, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[6, ]), amount = 0.002),
       lwd = 1.5, cex = 1.25, col = "black", pch = 92)


points(jitter(as.numeric(FDR.Dolgalev.30samples[7, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[7, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 95)


points(jitter(as.numeric(FDR.Dolgalev.30samples[8, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[8, ]), amount = 0.002),
       lwd = 1, cex = 1, col = "black", pch = 6)


points(jitter(as.numeric(FDR.Dolgalev.30samples[9, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[9, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "#1b9e77", pch = 1)


points(jitter(as.numeric(FDR.Dolgalev.30samples[1, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.30samples[1, ]), amount = 0.002),
       lwd =1.5, cex=1.25, pch = 1, col = "#d95f02")


abline(h = 0.2, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("topright",
       c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq",
         "ALDEx2", "NOISeq", "DSS", "dearseq"),
       pch = c(1, 2, 8, 43, 15,
               92, 95, 6, 1),
       pt.cex = c(1.25,1,1,1,1,1,1,1,1),
       cex= 0.95, pt.lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("#d95f02", "#1b9e77", "black", "#1b9e77", "#7570b3",
               "black", "black", "black", "#1b9e77"))
title(adj = 0, "(a)")


### 2*50 samples
plot(jitter(as.numeric(FDR.Dolgalev.50samples[4, ]), amount = 0.002),
     jitter(as.numeric(type2.Dolgalev.50samples[4, ]), amount = 0.002),
     lwd = 1.5, cex = 1, col = "#1b9e77", pch = 43,
     xlim = c(0, 0.7), ylim = c(-0.001, 0.8),
     xlab = "FDR", ylab = "Type II Error", cex.lab = 1.15)


points(jitter(as.numeric(FDR.Dolgalev.50samples[5, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[5, ]), amount = 0.002),
       lwd =1.5, cex=1, col = "#7570b3", pch =15)

points(jitter(as.numeric(FDR.Dolgalev.50samples[2, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[2, ]), amount = 0.002),
       col = "#1b9e77", lwd =1.5, cex=1, pch = 2 )

points(jitter(as.numeric(FDR.Dolgalev.50samples[3, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[3, ]), amount = 0.002),
       lwd =1.5, cex=0.75, col = "black", pch =8)


points(jitter(as.numeric(FDR.Dolgalev.50samples[6, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[6, ]), amount = 0.002),
       lwd = 1.5, cex = 1.25, col = "black", pch = 92)


points(jitter(as.numeric(FDR.Dolgalev.50samples[7, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[7, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "black", pch = 95)


points(jitter(as.numeric(FDR.Dolgalev.50samples[8, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[8, ]), amount = 0.002),
       lwd = 1, cex = 1, col = "black", pch = 6)


points(jitter(as.numeric(FDR.Dolgalev.50samples[9, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[9, ]), amount = 0.002),
       lwd = 1.5, cex = 1, col = "#1b9e77", pch = 1)

points(jitter(as.numeric(FDR.Dolgalev.50samples[1, ]), amount = 0.002),
       jitter(as.numeric(type2.Dolgalev.50samples[1, ]), amount = 0.002),
       lwd =1.5, cex=1.25, pch = 1, col = "#d95f02")



abline(h = 0.2, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("topright",
       c("SIEVE", "edgeR", "DESeq2", "voom", "tweeDEseq",
         "ALDEx2", "NOISeq", "DSS", "dearseq"),
       pch = c(1, 2, 8, 43, 15,
               92, 95, 6, 1),
       pt.cex = c(1.25,1,1,1,1,1,1,1,1),
       cex= 0.95, pt.lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("#d95f02", "#1b9e77", "black", "#1b9e77", "#7570b3",
               "black", "black", "black", "#1b9e77"))
title(adj = 0, "(b)")


######
# marginal prob.
data.frame(cbind(
  apply(FDR.Dolgalev.30samples, 1, function(k)  (sum(k<0.05)/30)*100),
  apply(FDR.Dolgalev.50samples, 1, function(k)  (sum(k<0.05)/30)*100)))


data.frame(cbind(
  apply(type2.Dolgalev.30samples, 1, function(k)  (sum(k<0.2)/30)*100),
  apply(type2.Dolgalev.50samples, 1, function(k)  (sum(k<0.2)/30)*100)))

# joint prob.
tab1 <- (FDR.Dolgalev.30samples < 0.05) + (type2.Dolgalev.30samples < 0.2)
tab2 <- (FDR.Dolgalev.50samples < 0.05) + (type2.Dolgalev.50samples < 0.2)

apply(tab1, 1, function(k)  (sum(k == 2)/30)*100)
apply(tab2, 1, function(k)  (sum(k == 2)/30)*100)

###END###