#######################################################
#SIEVE: One-stop differential expression, variability, 
#       and skewness analyses using RNA-Seq data
#Authors: Hongxiang Li and Dr. Tsung Fei Khang
#Email: chelsea.divo@hotmail.com (H. Li)
#Date: 08 April 2024
#Part 6: R Codes for DE, DV, and DS tests on Mayo 
#RNA-Seq dataset control vs. AD comparison
#######################################################

## R packages downloaded
# install.packages("readr")
# install.packages("devtools")
# install.packages("compositions")
# BiocManager::install("polyester")
# BiocManager::install("edgeR", force = T)
# devtools::install_github("Divo-Lee/SIEVE")
# install.packages("gamlss")
# BiocManager::install("cqn")
# devtools::install_github("zjdaye/MDSeq")
# install.packages("vioplot")
# install.packages("VennDiagram")
# install.packages("gridExtra")
# install.packages("httr")
# install.packages("jsonlite")
# install.packages("MASS")

library(readr); library(compositions); library(sn)
library(edgeR); library(DESeq2); library(limma)
library(tweeDEseq); library(ALDEx2); library(DSS)
library(vioplot); library(VennDiagram); library(gridExtra)
library(dplyr); library(httr); library(jsonlite)
library(SIEVE); library(MDSeq); library(gamlss)
library(MASS)

### Read data
# The use of the data files MayoRNAseq_individual_metadata_031422.csv
# and MayoRNAseq_RNAseq_TCX_geneCounts.csv requires permission from the data owners.
# Request permission from https://adknowledgeportal.synapse.org/ (look for Mayo RNAseq study)

# raw count table
ad_counts <- read.csv('...MayoRNAseq_RNAseq_TCX_geneCounts.csv', row.names = 1)
dim(ad_counts)

# meta-data
ad_meta <- read.csv('...MayoRNAseq_individual_metadata_031422.csv')
sum(is.na(ad_meta$diagnosis)) # check NA in disease name column
# remove samples which give NA in disease name column, in meta-data
ad_meta <- ad_meta[-which(is.na(ad_meta$diagnosis)), ]

# samples id in meta-data
control_id <- (ad_meta$individualID)[ad_meta$diagnosis == "control"]
ad_id <- (ad_meta$individualID)[ad_meta$diagnosis == "Alzheimer Disease"]
# samples id in raw count table
ad_counts_id <-  colnames(ad_counts)

control_vector <- sapply(ad_counts_id, function(k) k %in% control_id)
control_counts <-  ad_counts[, control_vector]
ad_vector <- sapply(ad_counts_id, function(k) k %in% ad_id)
ad_counts1 <- ad_counts[, ad_vector]
dim(control_counts); dim(ad_counts1)
N_ad <- length(control_counts) + length(ad_counts1)
mayo_counts1 <- as.matrix(cbind(control_counts, ad_counts1),
                          nrows=dim(ad_counts)[1], ncol = N_ad)

## Filter
CPM2 <- cpm(mayo_counts1)

keep <- rowMeans(CPM2[,1:length(control_counts)]) > 0.5 &
  rowMeans(CPM2[,(length(control_counts)+1):N_ad]) > 0.5 &
  apply(mayo_counts1[,1:length(control_counts)], 1, function(k) length(k[k == 0])/length(k)) < 0.85 &
  apply(mayo_counts1[,(length(control_counts)+1):N_ad], 1, function(k) length(k[k == 0])/length(k)) < 0.85

mayo_counts_filter <- mayo_counts1[keep, ]
dim(mayo_counts_filter) # 18664 genes left
                          



#######################
### DE Test Methods ###
#######################
group2 = c(rep(0, length(control_counts)), rep(1, length(ad_counts1))) # 78 control, 82 AD


###########################
### DE Test using SIEVE ###
###########################
# clr-transformation
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}

#########
### SIEVE
t1 <- proc.time()
clr_counts <- clr.transform(data = mayo_counts_filter)
clrSeq_result <- clrSeq(clr_counts, group = group2)
clrSIEVE_result <- clrSIEVE(clrSeq_result = clrSeq_result,
                            alpha_level = 0.05,
                            order_DE = F,
                            order_LFC = F,
                            order_DS = F,
                            order_sieve = F)
 # clrSIEVE() simultaneously detects DE, DV, and DS genes
 # clrSIEVE() results in 4 lists
 # clrSIEVE_result includes the results of DE, DV, and DS tests

 # clrDE is the DE test in SIEVE
 # function clrDE() only identifies DE genes
clrDE_result <- clrSIEVE_result$clrDE_test # DE test result
 # clrDE_result <- na.omit(clrDE_result) # here, no NA value

# apply Xiao et al (2014) method
criteria.1 <- (-log10(0.05))*quantile(clrDE_result$DE, 0.025)
criteria.2 <- (-log10(0.05))*quantile(clrDE_result$DE, 0.975)
criteria_clrDE <- ((-log10(clrDE_result$adj_pval_DE))*clrDE_result$DE < criteria.1 |
                     (-log10(clrDE_result$adj_pval_DE))*clrDE_result$DE > criteria.2)
 # sum(criteria_clrDE) # number of DE genes called

DE_table_clrDE_AD <- clrDE_result[criteria_clrDE, ]
DE_genes_clrDE_AD <- row.names(DE_table_clrDE_AD)
as.numeric(proc.time() - t1)[3] # run time, in seconds

dim(clrDE_result) # it shows that no NA value in clrDE test
length(DE_genes_clrDE_AD) # number of DE genes called



## volcano plot, DE genes, SIEVE, AD vs. control
plot(clrDE_result$DE, -log10(clrDE_result$adj_pval_DE),
     xlim = c(-1.5, 1.5),
     ylim = c(0, 15),
     xlab = expression(paste(Delta, hat(mu))),
     ylab = expression(-log[10](p)),
     col = c(rgb(0,0,0,0.3)), pch=1, lwd = 1.5, cex = 0.6)
points(clrDE_result$DE[which(-log10(clrDE_result$adj_pval_DE) > 15)],
       rep(15, dim(clrDE_result[-log10(clrDE_result$adj_pval_DE) > 15, ])[1]),
       col = c(rgb(0,0,0,0.3)),
       pch = 2, lwd = 1.5, cex = 0.6)
points(rep(1.5, dim(clrDE_result[clrDE_result$DE > 1.5, ])[1]),
       -log10(clrDE_result$adj_pval_DE[which(clrDE_result$DE > 1.5)]),
       col = c(rgb(0,0,0,0.3)),
       pch = 6, lwd = 1.5, cex = 0.6)
abline(h = -log10(0.05), lty = 2, lwd = 1.5, col = "blue")
abline(v = quantile(clrDE_result$DE , c(0.025, 0.975)),
       col = "blue", lty = 2, lwd = 1.5)

x1 <- seq(-2, 0, 0.001); x2 <- seq(0, 2, 0.001)
lines(x1, criteria.1/x1, col = "red", lwd = 1.25)
lines(x2, criteria.2/x2, col = "red", lwd = 1.25)




############################
### competing DE methods ###
############################

########
## edgeR
data1 <- mayo_counts_filter 
 # hereafter, we use data1 instead of using mayo_counts_filter
design = model.matrix(~group2)

t2 <- proc.time()
libsizes <- colSums(data1)
nf <- calcNormFactors(data1, method="TMM") # normalization factors for edgeR, limma

dat.edgeR <- DGEList(counts=data1, norm.factors=nf, group=group2)
dat.edgeR <- estimateDisp(dat.edgeR, design)
fit.edgeR <- glmQLFit(dat.edgeR, design)
test.edgeR <- glmQLFTest(fit.edgeR, coef=2)
tp <- topTags(test.edgeR, n = Inf)

criteria.1_edgeR <- quantile(tp$table$logFC, 0.025)*(-log10(0.05))
criteria.2_edgeR <- quantile(tp$table$logFC, 0.975)*(-log10(0.05))
criteria_edgeR <- (tp$table$logFC*(-log10(tp$table$FDR)) < criteria.1_edgeR | tp$table$logFC*(-log10(tp$table$FDR)) > criteria.2_edgeR)
sum(criteria_edgeR)

DE_table_edgeR <- tp$table[criteria_edgeR, ]
DE_genes_edgeR <- row.names(DE_table_edgeR)
as.numeric(proc.time() - t2)[3] # run time, in seconds

length(DE_genes_edgeR) # 3942 DE genes called by edgeR
 # length(intersect(DE_genes_edgeR, DE_genes_clrDE_AD))



#########
## DESeq2
t3 <- proc.time()
libsizes <- colSums(data1)
nf <- calcNormFactors(data1, method="TMM") # normalization factors for edgeR, limma
els <- nf*libsizes # effective library sizes
sf <- els/exp(mean(log(libsizes))) # size factors for DESeq2, factors should multiple to 1

dat.DESeq2 <- DESeqDataSetFromMatrix(countData = data1, colData = data.frame(group2),
                                     design = design)
sizeFactors(dat.DESeq2) <- sf
fit.DESeq2 <- DESeq(dat.DESeq2, minReplicatesForReplace=Inf)
res.DESeq2 <- results(fit.DESeq2, cooksCutoff=F, alpha=0.05)

criteria.1_DESeq2 <- quantile(res.DESeq2$log2FoldChange,0.025)*(-log10(0.05))
criteria.2_DESeq2 <- quantile(res.DESeq2$log2FoldChange,0.975)*(-log10(0.05))
criteria_DESeq2 <- (res.DESeq2$log2FoldChange*(-log10(res.DESeq2$padj)) < criteria.1_DESeq2 | res.DESeq2$log2FoldChange*(-log10(res.DESeq2$padj)) > criteria.2_DESeq2)

DE_table_DESeq2 <- res.DESeq2[criteria_DESeq2, ]
DE_genes_DESeq2 <- row.names(DE_table_DESeq2)
as.numeric(proc.time() - t3)[3]

length(DE_genes_DESeq2) # 4155 DE genes called by DESeq2
 # length(intersect(DE_genes_DESeq2, DE_genes_clrDE_AD))



##############
### limma-voom
t4 <- proc.time()
design = model.matrix(~group2)
# Normalisation & Transformation
libsizes <- colSums(data1)
nf <- calcNormFactors(data1, method="TMM")

dat.edgeR <- DGEList(counts=data1, norm.factors=nf, group=group2)
dat.edgeR <- estimateDisp(dat.edgeR, design)

dat.voom <- voom(dat.edgeR) # uses same data object as edgeR
fit.voom <- lmFit(dat.voom, design)
res.voom <- eBayes(fit.voom)
 # sum(is.na(res.voom)) 
top.table <- topTable(res.voom, sort.by = "P", n = Inf)

criteria.1_voom <- quantile(top.table$logFC, 0.025)*(-log10(0.05))
criteria.2_voom <- quantile(top.table$logFC, 0.975)*(-log10(0.05))
criteria_voom <- (top.table$logFC*(-log10(top.table$adj.P.Val)) < criteria.1_voom | top.table$logFC*(-log10(top.table$adj.P.Val)) > criteria.2_voom)

DE_table_voom <- top.table[criteria_voom, ]
DE_genes_voom <- row.names(DE_table_voom)
as.numeric(proc.time() - t4)[3]

length(DE_genes_voom) # 3636 DE genes called by voom
 # length(intersect(DE_genes_voom, DE_genes_clrDE_AD))



#############
### tweeDEseq
t5 <- proc.time()
libsizes <- colSums(data1)
nf <- calcNormFactors(data1, method="TMM")
els <- nf * libsizes
sf <- els / exp(mean(log(libsizes)))
norm.data1 <- t(t(data1) / sf)

counts_tweeDEseq <- as.matrix(norm.data1)
resPT <- tweeDE(counts_tweeDEseq, group = group2)

criteria.1_tweeDEseq <- quantile(resPT$log2fc, 0.025)*(-log10(0.05))
criteria.2_tweeDEseq <- quantile(resPT$log2fc, 0.975)*(-log10(0.05))
criteria_tweeDEseq <- (resPT$log2fc*(-log10(resPT$pval.adjust)) < criteria.1_tweeDEseq | resPT$log2fc*(-log10(resPT$pval.adjust)) > criteria.2_tweeDEseq)
DE_table_twee <- resPT[criteria_tweeDEseq, ]
DE_genes_twee <-  row.names(de_table_twee)

as.numeric(proc.time() - t5)[3] # about 2.5 hours

length(DE_genes_twee) # 3553 DE genes called by tweeDEseq
 # length(intersect(DE_genes_clrDE_AD, DE_genes_twee))




########
## ALDEx2
t6 <- proc.time()
conds <- as.character(group2)
y <- aldex(data1, conds, mc.samples=128, denom="all",test="t", effect=TRUE,  paired.test=FALSE)

criteria.1_ALDEx2 <- quantile(y$diff.btw, 0.025)*(-log10(0.05))
criteria.2_ALDEx2 <- quantile(y$diff.btw, 0.975)*(-log10(0.05))
criteria_ALDEx2 <- (y$diff.btw*(-log10(y$we.eBH)) < criteria.1_ALDEx2 | y$diff.btw*(-log10(y$we.eBH)) > criteria.2_ALDEx2)
DE_table_ALDEx2 <- y[criteria_ALDEx2,]
DE_genes_ALDEx2 <- row.names(de_table_ALDEx2)
as.numeric(proc.time() - t6)[3] # about 37 minutes

length(DE_genes_ALDEx2) # 3230 DE genes called by ALDEx2
 # length(intersect(DE_genes_ALDEx2, DE_genes_clrDE_AD))




#######
### DSS
colnames(data1) <- NULL
t7 <- proc.time()
seqData_AD <- newSeqCountSet(data1, group2)
seqData_AD <- estNormFactors(seqData_AD, method="lr") # similar to TMM. No TMM choice in DSS package.
seqData_AD <- estDispersion(seqData_AD)
result_DSS_AD <- waldTest(seqData_AD, 0, 1)

criteria.1_DSS <- quantile(result_DSS_AD$lfc, 0.025)*(-log10(0.05))
criteria.2_DSS <- quantile(result_DSS_AD$lfc, 0.975)*(-log10(0.05))
criteria_DSS <- (result_DSS_AD$lfc*(-log10(result_DSS_AD$fdr)) < criteria.1_DSS |  result_DSS_AD$lfc*(-log10(result_DSS_AD$fdr)) > criteria.2_DSS)
sum(criteria_DSS)

DE_table_DSS <- result_DSS_AD[criteria_DSS, ]
DE_genes_DSS <- row.names(DE_table_DSS)
as.numeric(proc.time() - t7)[3]

length(DE_genes_DSS) # 4135 DE genes called by DSS
 # length(intersect(DE_genes_DSS, DE_genes_clrDE_AD))



###################
### DE Analysis ###
###################

### DE genes list, detected by SIEVE
DE_genes_list_AD <- DE_table_clrDE_AD[order(DE_table_clrDE_AD$DE), ]
DE_genes_list_AD <- select(DE_genes_list_AD, -"de_indicator")
#DE_genes_id_AD <- row.names(DE_genes_list_AD)
## Gene id convert to gene symbol
#url = "https://biotools.fr/human/ensembl_symbol_converter/"

#ids_AD = DE_genes_id_AD
#ids_json_AD <- toJSON(ids_AD)
#body_AD <- list(api=1, ids=ids_json_AD)
#r_AD <- POST(url, body = body_AD)
#output_AD <- fromJSON(content(r_AD, "text"), flatten=TRUE)
#output_AD[sapply(output_AD, is.null)] <- NA # some genes do not have gene symbol
#gene_symbol_AD <- unlist(output_AD, use.names = FALSE)
#for (i in 1:length(ids_AD)) {
#  if(is.na(gene_symbol_AD[i])){
#    gene_symbol_AD[i] <- ids_AD[i]
#  }
#}
#DE_genes_list_AD <- cbind.data.frame(gene_symbol = gene_symbol_AD, DE_genes_list_AD)
# head(DE_genes_list_AD); tail(DE_genes_list_AD)

# write.csv(DE_genes_list_AD, "AD_DE_genes_table_SIEVE.csv", quote=FALSE, row.names=T)



##################################
### Union and unique DE genes list
## Union DE genes
## DE_genes_union_AD (Table S9)
union_DE_genes_AD <- unique(c(DE_genes_clrDE_AD, # SIEVE
                              DE_genes_edgeR,
                              DE_genes_DESeq2,
                              DE_genes_voom,
                              DE_genes_twee,
                              DE_genes_ALDEx2,
                              DE_genes_DSS))
length(union_DE_genes_AD)

AD_union_DE_table <- clrDE_result[union_DE_genes_AD, ]
AD_union_DE_table$SIEVE_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_clrDE_AD)
AD_union_DE_table$edgeR_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_edgeR)
AD_union_DE_table$DESeq2_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_DESeq2)
AD_union_DE_table$voom_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_voom)
AD_union_DE_table$tweeDEseq_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_twee)
AD_union_DE_table$ALDEx2_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_ALDEx2)
AD_union_DE_table$DSS_indicator <- as.integer(union_DE_genes_AD %in% DE_genes_DSS)
AD_union_DE_table <- select(AD_union_DE_table, -"de_indicator")

# convert gene id to gene symbol
url = "https://www.biotools.fr/human/ensembl_symbol_converter/"
ids_AD_union = union_DE_genes_AD
ids_json_AD_union <- toJSON(ids_AD_union)
body_AD_union <- list(api=1, ids=ids_json_AD_union)
r_AD_union <- POST(url, body = body_AD_union)
output_AD_union <- fromJSON(content(r_AD_union, "text"), flatten=TRUE)
output_AD_union[sapply(output_AD_union, is.null)] <- NA # some genes do not have gene symbol
gene_symbol_AD_union <- unlist(output_AD_union, use.names = FALSE)
for (i in 1:length(ids_AD_union)) {
  if(is.na(gene_symbol_AD_union[i])){
    gene_symbol_AD_union[i] <- ids_AD_union[i]
  }
}
AD_union_DE_table <- cbind.data.frame(gene_symbol = gene_symbol_AD_union, AD_union_DE_table)
AD_union_DE_table <- AD_union_DE_table[order(AD_union_DE_table$DE),]

#head(AD_union_DE_table) 
#tail(AD_union_DE_table)
# online supplementary table S9
# write.csv(AD_union_DE_table, "DE_genes_union_AD (Table S9).csv", quote=FALSE, row.names=T)
 


## Unique DE genes, detected by SIEVE
Unique_SIEVE_AD <- AD_union_DE_table[(AD_union_DE_table$SIEVE_indicator == 1 &
                                      AD_union_DE_table$edgeR_indicator == 0 &
                                      AD_union_DE_table$DESeq2_indicator == 0 &
                                      AD_union_DE_table$voom_indicator == 0 &
                                      AD_union_DE_table$tweeDEseq_indicator == 0 &
                                      AD_union_DE_table$ALDEx2_indicator == 0 &
                                      AD_union_DE_table$DSS_indicator == 0), ]
Unique_SIEVE_AD <- select(Unique_SIEVE_AD, gene_symbol:mu2)
dim(Unique_SIEVE_AD) # 24 DE genes uniquely detected by SIEVE
Unique_SIEVE_AD 



## High confidence DE genes
overlap_DE_genes_AD <- AD_union_DE_table[(AD_union_DE_table$SIEVE_indicator == 1 &
                                          AD_union_DE_table$edgeR_indicator == 1 &
                                          AD_union_DE_table$DESeq2_indicator == 1 &
                                          AD_union_DE_table$voom_indicator == 1 &
                                          AD_union_DE_table$tweeDEseq_indicator == 1 &
                                          AD_union_DE_table$ALDEx2_indicator == 1 &
                                          AD_union_DE_table$DSS_indicator == 1), ]


dim(overlap_DE_genes_AD)
 # 2439 high-confidence DE genes


###################
### DE Venn Diagram
vd1 <- venn.diagram(
  x = list(DE_genes_clrDE_AD,
           DE_genes_edgeR,
           DE_genes_DESeq2,
           DE_genes_voom),
  category.names = c("SIEVE" , "edgeR",  "DESeq2", "voom"),
  fill=c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  cat.col="black",
  filename = NULL
)

#
vd2 <- draw.pairwise.venn(length(DE_genes_clrDE_AD),
                          length(DE_genes_ALDEx2),
                          length(intersect(DE_genes_clrDE_AD,DE_genes_ALDEx2)),
                          category = c("SIEVE", "ALDEx2"),
                          fill = c("#0073C2FF", "green"),
                          cat.pos = c(0, 0),
                          alpha = rep(0.5, 2),
                          scaled = FALSE)

#
vd3 <- draw.pairwise.venn(length(DE_genes_clrDE_AD),
                          length(DE_genes_twee),
                          length(intersect(DE_genes_clrDE_AD,DE_genes_twee)),
                          category = c("SIEVE", "tweeDEseq"),
                          fill = c("#0073C2FF", "orange"),
                          cat.pos = c(0, 0),
                          alpha = rep(0.5, 2),
                          scaled = FALSE)


vd4 <- draw.pairwise.venn(length(DE_genes_clrDE_AD),
                          length(DE_genes_DSS),
                          length(intersect(DE_genes_clrDE_AD,DE_genes_DSS)),
                          category = c("SIEVE", "DSS"),
                          fill = c("#0073C2FF", "pink1"),
                          cat.pos = c(0, 0),
                          alpha = rep(0.5, 2),
                          scaled = FALSE)

grid.arrange(gTree(children=vd1), top=textGrob(expression(bold("(a)")), x = 0.1, hjust = 0))
grid.arrange(gTree(children=vd2), top=textGrob(expression(bold("(b)")), x = 0.1, hjust = 0))
grid.arrange(gTree(children=vd3), top=textGrob(expression(bold("(c)")), x = 0.1, hjust = 0))
grid.arrange(gTree(children=vd4), top=textGrob(expression(bold("(d)")), x = 0.1, hjust = 0))
  ### DE Analysis End ###




###########################
### DV Methods Analysis ###
###########################

### DV Test using SIEVE
 # call the DV test result from "clrSIEVE_result"
SIEVE_DV_test_AD <- clrSIEVE_result$clrDV_test
 # head(SIEVE_DV_test_AD) # all observed genes
 # clrDV/clrDV() is the DV test in SIEVE and clrDV package

dv_criteria.1 <- -log10(0.05)*(-1.25) # LFC_DV = -1.25
dv_criteria.2 <- -log10(0.05)*1.25 # LFC_DV = 1.25
SIEVE_dv_criteria <- ((-log10(SIEVE_DV_test_AD$adj_pval_DV))*SIEVE_DV_test_AD$LFC < dv_criteria.1 | (-log10(SIEVE_DV_test_AD$adj_pval_DV))*SIEVE_DV_test_AD$LFC > dv_criteria.2)
 sum(SIEVE_dv_criteria)

# call DV genes  
DV_table_SIEVE <- SIEVE_DV_test_AD[SIEVE_dv_criteria, ]
  # head(DV_table_SIEVE); tail(DV_table_SIEVE)
DV_genes_SIEVE <- row.names(DV_table_SIEVE)
length(DV_genes_SIEVE)

 # sum(DV_table_SIEVE$LFC > 0) 
  # 15 genes
 # sum(DV_table_SIEVE$LFC < 0) 
  # 2261 genes
 

## Volcano plot, DV genes, SIEVE
plot(SIEVE_DV_test_AD$LFC,
     -log10(SIEVE_DV_test_AD$adj_pval_DV),
     xlim = c(-2, 2),
     xlab = expression(hat(psi)),
     ylab = expression(-log[10](p)),
     col = c(rgb(0,0,0,0.3)), pch=1, lwd = 1.5, cex = 0.6)
abline(h = -log10(0.05), lty = 2, lwd = 1.5, col = "blue")
abline(v = c(-1.25,1.25), col = "blue", lty = 2, lwd = 1.5)

x1 <- seq(-2, 0, 0.001); x2 <- seq(0, 2, 0.001)
lines(x1, dv_criteria.1/x1, col = "red", lwd = 1.25)
lines(x2, dv_criteria.2/x2, col = "red", lwd = 1.25)
 # title(adj = 0, "(b)")






###########################
### DS Test using SIEVE ###
###########################
DS_test_AD <- clrSIEVE_result$clrDS_test
# sum(DS_test_AD$adj_pval < 0.05) # only consider statistical significance

ds_criteria.1 <- (-log10(0.05))*quantile(DS_test_AD$DS, 0.025)
ds_criteria.2 <- (-log10(0.05))*quantile(DS_test_AD$DS, 0.975)
ds_criteria <- ((-log10(DS_test_AD$adj_pval_DS))*DS_test_AD$DS < ds_criteria.1 | (-log10(DS_test_AD$adj_pval_DS))*DS_test_AD$DS > ds_criteria.2)
sum(ds_criteria) # 1204 DS genes called by SIEVE


### Online supplementary table
### SIEVE_DE_DV_DS_genes_control_AD (Table S10)
full_table <- clrSIEVE_result$clrSIEVE_tests
 # the list, clrSIEVE_tests, gives DE, DV and DS tests in one table
 # all the 18664 genes

# bs: biological and statistical significance
# de_bs, dv_bs, ds_bs are indicators by using Xiao et al (2014) method
full_table$de_bs <- as.integer(criteria_clrDE)
full_table$dv_bs <- as.integer(SIEVE_dv_criteria)
full_table$ds_bs <- as.integer(ds_criteria)

full_table <- full_table[order(full_table$DS), ] # ordered by DS
genes_id_AD <- row.names(full_table)

## Gene id convert to gene symbol
url = "https://biotools.fr/human/ensembl_symbol_converter/"
ids_AD = genes_id_AD
ids_json_AD <- toJSON(ids_AD)
body_AD <- list(api=1, ids=ids_json_AD)
r_AD <- POST(url, body = body_AD)
output_AD <- fromJSON(content(r_AD, "text"), flatten=TRUE)
output_AD[sapply(output_AD, is.null)] <- NA # some genes do not have gene symbol
gene_symbol_AD <- unlist(output_AD, use.names = FALSE)
for (i in 1:length(ids_AD)) {
  if(is.na(gene_symbol_AD[i])){
    gene_symbol_AD[i] <- ids_AD[i]
  }
}
full_table <- cbind.data.frame(gene_symbol = gene_symbol_AD, full_table)
 # head(full_table); tail(full_table)
# Supplementary Table S10
# write.csv(full_table, "SIEVE_DE_DV_DS_genes_control_AD (Table S10).csv", quote=FALSE, row.names=T)

# ordered by SD_ratio (or, LFC for variability)
#full_table2 <- full_table1[order(full_table1$SD_ratio), ]
#head(full_table2, 10)

# ordered by DE
#full_table3 <- full_table1[order(full_table1$DE), ]
#head(full_table3, 10)



###########################
## Volcano plot for DS test
plot(DS_test_AD$DS[which(-log10(DS_test_AD$adj_pval_DS) <= 15)],
     -log10(DS_test_AD$adj_pval_DS[which(-log10(DS_test_AD$adj_pval_DS) <= 15)]),
     xlim = c(-2, 2),
     ylim = c(0, 15),
     xlab = expression(paste(Delta, hat(gamma))),
     ylab = expression(-log[10](p)),
     col = c(rgb(0,0,0,0.3)), pch=1, lwd = 1.5, cex = 0.6)
points(DS_test_AD$DS[which(-log10(DS_test_AD$adj_pval_DS) > 15)],
       rep(15, dim(DS_test_AD[-log10(DS_test_AD$adj_pval_DS) > 15, ])[1]),
       col = c(rgb(0,0,0,0.3)),
       pch = 17, lwd = 1.5, cex = 0.6)
abline(h = -log10(0.05), lty = 2, lwd = 1.5, col = "blue")
abline(v = quantile(DS_test_AD$DS , c(0.025, 0.975)),
       col = "blue", lty = 2, lwd = 1.5)
x1 <- seq(-2, 0, 0.001); x2 <- seq(0, 2, 0.001)
lines(x1, ds_criteria.1/x1, col = "red", lwd = 1.25)
lines(x2, ds_criteria.2/x2, col = "red", lwd = 1.25)
#title(adj = 0, "(c)")



######
# Fig. Frequency of estimated skewness papameters

#par(mfrow=c(1,2))
#esti_skew_control <- clrSeq_result$gamma1
#esti_skew_AD <- clrSeq_result$gamma2

#hist(esti_skew_control, breaks = 20,
#     main = NULL, probability = T,
#     ylim = c(0,1),
#     xlab = expression(paste(hat(gamma))))
#title(adj=0, "(a)")

#hist(esti_skew_AD, breaks = 20,
#     main = NULL,  probability = T,
#     ylim = c(0,1),
#     xlab = expression(paste(hat(gamma))))
#title(adj=0, "(b)")



######
# Fig. Distribution of 8 classes of gene (converted this Fig. to Table 3 already)
genes_000 <- row.names(full_table[full_table$de_bs == 0 & full_table$dv_bs==0 & full_table$ds_bs==0 ,])
genes_100 <- row.names(full_table[full_table$de_bs == 1 & full_table$dv_bs==0 & full_table$ds_bs==0 ,])
genes_010 <- row.names(full_table[full_table$de_bs == 0 & full_table$dv_bs==1 & full_table$ds_bs==0 ,])
genes_001 <- row.names(full_table[full_table$de_bs == 0 & full_table$dv_bs==0 & full_table$ds_bs==1 ,])
genes_110 <- row.names(full_table[full_table$de_bs == 1 & full_table$dv_bs==1 & full_table$ds_bs==0 ,])
genes_101 <- row.names(full_table[full_table$de_bs == 1 & full_table$dv_bs==0 & full_table$ds_bs==1 ,])
genes_011 <- row.names(full_table[full_table$de_bs == 0 & full_table$dv_bs==1 & full_table$ds_bs==1 ,])
genes_111 <- row.names(full_table[full_table$de_bs == 1 & full_table$dv_bs==1 & full_table$ds_bs==1 ,])
# length(genes_000)+length(genes_100)+length(genes_010)+length(genes_001)+length(genes_110)+length(genes_101)+length(genes_011)+length(genes_111)

x <- c("000", "100", "010", "001",
       "110", "101", "011", "111")
freq <- c(length(genes_000), length(genes_100),
          length(genes_010), length(genes_001),
          length(genes_110), length(genes_101),
          length(genes_011), length(genes_111))

b <- barplot(freq, names = x, xlab= "Gene classes", 
             ylab = "Frequency",
             ylim = c(0, 14000), 
             border ="grey83",
             col = "grey83")
text(b, freq + 500, freq, font=6, col="black")
# (converted this Fig. to Table 3 already)



#################
### Contour Plots (Fig. 7)
x1 <- esti_mean_control <- clrSeq_result$mu1
x2 <- esti_mean_AD <- clrSeq_result$mu2
y1 <- esti_sd_control <- clrSeq_result$sigma1
y2 <- esti_sd_AD <- clrSeq_result$sigma2
z1 <- esti_skew_control <- clrSeq_result$gamma1
z2 <- esti_skew_AD <- clrSeq_result$gamma2
# range(x1); range(log2(y1)); range(z1)
# range(x2); range(log2(y2)); range(z2)


xx_2 <- kde2d(x1, x2)  # kde2d: Two-Dimensional Kernel Density Estimation (MASS package)
yy_log2_2 <- kde2d(log2(y1), log2(y2))
zz_2 <- kde2d(z1, z2)

### mu1 vs. mu2
par(mar= c(5, 4.6, 4 ,1))

ffilled.contour(xx_2,
               xlab = expression(paste(hat(mu)[1])),
               ylab = expression(paste(hat(mu)[2])),
               color.palette = function(n) hcl.colors(n=10, "YlOrRd", rev = T),
               nlevels = 10,
               plot.axes={
                 points(x1[which(full_table$de_bs==1)], 
                        x2[which(full_table$de_bs==1)], 
                        col = c(rgb(0,0,0,0.1)), 
                        pch=".", lwd = 1, cex = 2.75)
                 axis(side=1)
                 axis(side=2)
               }
)
title(adj=0, "(a)") 

### log2(sigma1) vs. log2(sigma2)
filled.contour(yy_log2_2,
               xlab = expression(paste(log[2](hat(sigma)[1]))),
               ylab = expression(paste(log[2](hat(sigma)[2]))),
              # color.palette = function(n) hcl.colors(n=14, "YlOrRd", rev = T),
               xlim = c(log2(min(y1)), log2(max(y1))),
               ylim = c(log2(min(y2)), log2(max(y2))),
               #nlevels = 14,
               plot.axes={
                 points(y1[which(full_table$dv_bs==1)], y2[which(full_table$dv_bs==1)], 
                        col = c(rgb(0,0,0,0.075)), 
                        pch=".", lwd = 2, cex = 3.5)
                 axis(side=1)
                 axis(side=2)
               }
)
title(adj=0, "(b)")


### gamma1 vs. gamma2
filled.contour(zz_2, 
               xlab = expression(paste(hat(gamma)[1])),
               ylab = expression(paste(hat(gamma)[2])),
               # color.palette = function(n) hcl.colors(n=20, "YlOrRd", rev = T),
               xlim = c(-1, 1),
               ylim = c(-1, 1),
               plot.axes={
                 points(z1[which(full_table$ds_bs==1)], z2[which(full_table$ds_bs==1)], 
                        col = c(rgb(0,0,0,0.1)), 
                        pch=".", lwd = 1, cex = 3.5)
                 axis(side=1)
                 axis(side=2)
               }
)
title(adj=0, "(c)")


###################################
### Violin plots for selected genes
### Supplementary Figures
###################################
### Pure DE genes
DE_100_class_genes <- c("ENSG00000147571", "ENSG00000157005",
                        "ENSG00000163221", "ENSG00000188488")
# ENSG00000147571	CRH
# ENSG00000157005	SST
# ENSG00000163221	S100A12
# ENSG00000188488	SERPINA5
labels <- c("(a)","(b)","(c)","(d)")

for (i in 1:length(DE_100_class_genes)) {
  violin.plot.SIEVE(data = clr_counts, DE_100_class_genes[i],
                    group = group2,
                    group.names = c("control", "AD"))
  title(adj=0, labels[i])
}


### DE genes uniquely detected by SIEVE
DE_unique_genes <- c("ENSG00000176165", # FOXG1
                     "ENSG00000168830", # HTR1E
                     "ENSG00000112855", # HARS2
                     "ENSG00000204371", # EHMT2
                     "ENSG00000103197", # TSC2
                     "ENSG00000130684") # ZNF337

# uniquely detected
# ENSG00000176165	FOXG1 # 110 class
# ENSG00000168830	HTR1E
# ENSG00000130684	ZNF337
# ENSG00000204371	EHMT2
# ENSG00000112855	HARS2 # 110 class
# ENSG00000103197	TSC2
labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")


for (i in 1:length(DE_unique_genes)) {
  violin.plot.SIEVE(data = clr_counts, DE_unique_genes[i],
                    group = group2,
                    group.names = c("control", "AD"))
  title(adj=0, labels[i])
}



### Pure DV genes
DV_010_class_genes <- c("ENSG00000172349", 
                        "ENSG00000104969",
                        "ENSG00000203618")
# ENSG00000172349	IL16
# ENSG00000104969	SGTA
# ENSG00000203618	GP1BB
lables <- c("(a)","(b)","(c)")

for (i in 1:length(DV_010_class_genes)) {
  violin.plot.SIEVE(data = clr_counts, 
                    DV_010_class_genes[i],
                    group = group2,
                    group.names = c("control", "AD"))
  title(adj=0, labels[i])
}



### Pure DS genes
DS_001_class_gene <- c("ENSG00000198668", "ENSG00000121297",
                       "ENSG00000051523", "ENSG00000164050",
                       "ENSG00000188603", "ENSG00000185736",
                       "ENSG00000104774", "ENSG00000278768")

# ENSG00000198668	CALM1
# ENSG00000121297	TSHZ3
# ENSG00000051523	CYBA
# ENSG00000164050	PLXNB1
# ENSG00000188603	CLN3
# ENSG00000185736	ADARB2
# ENSG00000104774	MAN2B1
# ENSG00000278768	BACE1-AS

labels <- c("(a)","(b)","(c)","(d)",
            "(e)","(f)","(g)","(h)")

for (i in 1:length(DS_001_class_gene)) {
  violin.plot.SIEVE(data = clr_counts, DS_001_class_gene[i],
                    group = group2,
                    group.names = c("control", "AD"))
  title(adj=0, labels[i])
}



###
### simultaneous DE, DV and DS analyses
simultaneous_genes <- c("ENSG00000162747", 
                        "ENSG00000128564", 
                        "ENSG00000139117", 
                        "ENSG00000017427", 
                        "ENSG00000116678", 
                        "ENSG00000153266",
                        "ENSG00000134532",
                        "ENSG00000133110") 


# ENSG00000162747	FCGR3B  # 110 class
# ENSG00000128564	VGF     # 101 class
# ENSG00000139117	CPNE8   # 101 class 
# ENSG00000017427	IGF1    # 101 class
# ENSG00000116678	LEPR    # 101 class
# ENSG00000153266	FEZF2   # 111 class 
# ENSG00000134532	SOX5    # 111 class
# ENSG00000133110	POSTN   # 111 class
labels <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)")


for (i in 1:length(simultaneous_genes)) {
  violin.plot.SIEVE(data = clr_counts, simultaneous_genes[i],
                    group = group2,
                    group.names = c("control", "AD"))
  title(adj=0, labels[i])
}


###END##
