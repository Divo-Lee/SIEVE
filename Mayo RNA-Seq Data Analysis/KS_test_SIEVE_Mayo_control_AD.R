#########################################################
#SIEVE: One-stop differential expression, variability, 
#       and skewness analyses using RNA-Seq data
#Author: Hongxiang Li and by Tsung Fei Khang
#Email: chelsea.divo@hotmail.com (H. Li)
#Latest update: 23 March 2024
#Part 5: R Codes for DE, DV, and DS tests on Mayo RNA-Seq
#dataset Goodness-of-fit Assessment of the skew-normal model  
#on the distribution of the CLR-transformed Mayo RNA-Seq
#data (control vs. AD) 
#########################################################

## R packages downloaded
# install.packages("BiocManager")
# install.packages("devtools")
# install.packages("compositions")
# BiocManager::install("edgeR", force = T)
# devtools::install_github("Divo-Lee/SIEVE")
# install.packages("sn")
library(compositions); library(edgeR)
library(SIEVE); library(sn)


## Read data
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
dim(mayo_counts_filter)


################
group2 = c(rep(0, length(control_counts)), rep(1, length(ad_counts1))) # 78 control, 82 AD

# CLR-transformation
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}

##############
### MLE/MPLE
t1 <- proc.time()
clr.counts2 <- clr.transform(data = mayo_counts_filter)
clrseq_AD <- clrSeq(clr.counts2, group = group2)
as.numeric(proc.time() - t1)[3] # computing time, in seconds

control_sn_CP <- clrseq_AD[, c("mu1", "sigma1", "gamma1")]
AD_sn_CP <- clrseq_AD[, c("mu2", "sigma2", "gamma2")]


########################################################
## map centered parameters (CP) to direct parameters (DP)
cp_to_dp <- function(mean=NULL, sd=NULL, skewness=NULL){
  b <- sqrt(2/pi)
  
  if(skewness >= 0){
    r <- (2*skewness/(4-pi))^(1/3)
  } else {
    r <- -(2*(- skewness)/(4-pi))^(1/3)
  }
  
  alpha <- r/sqrt(2/pi - (1-2/pi)*r^2)
  delta <- alpha/sqrt(1 + alpha^2)
  omega <- sd/sqrt(1 - (b^2)*delta^2)
  xi <- mean - b*omega*delta
  return(c(xi, omega, alpha))
}


###
control_sn_DP <- matrix(NA, nrow = dim(control_sn_CP)[1], ncol = 3)
AD_sn_DP <- matrix(NA, nrow = dim(AD_sn_CP)[1], ncol = 3)

for (i in 1:dim(control_sn_CP)[1]) {
  control_sn_DP[i,] <- c(cp_to_dp(control_sn_CP[i,1],
                                  control_sn_CP[i,2],
                                  control_sn_CP[i,3]))
}

colnames(control_sn_DP) <- c("xi1", "omega1", "alpha1")
control_sn_DP <- as.data.frame(control_sn_DP)


for (i in 1:dim(AD_sn_CP)[1]) {
  AD_sn_DP[i,] <- c(cp_to_dp(AD_sn_CP[i,1],
                             AD_sn_CP[i,2],
                             AD_sn_CP[i,3]))
}

colnames(AD_sn_DP) <- c("xi2", "omega2", "alpha2")
AD_sn_DP <- as.data.frame(AD_sn_DP)


# KS test for the control group (control vs. AD)
control_KS_test_pvalue <- vector()
for (i in 1:dim(control_sn_CP)[1]) {
  ks <- ks.test(clr.counts2[i, 1:78],
                "psn",
                xi = control_sn_DP$xi1[i],
                omega = control_sn_DP$omega1[i],
                alpha = control_sn_DP$alpha1[i])
  control_KS_test_pvalue[i] <- ks$p.value
}
sum(control_KS_test_pvalue < 0.05)
# 279/18664 = 0.015
# the skew-normal distribution fits 98.5% of genes in the control group 

# Fig. (a) control group
hist(control_KS_test_pvalue,
     freq = F, breaks = 20,
     xlab =expression(p), ylim = c(0, 5),
     main = NULL, border = "grey83")
abline(v = 0.05, lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(a)")


## KS test for the AD group (control vs. AD)
AD_KS_test_pvalue <- vector()
for (i in 1:dim(control_sn_CP)[1]) {
  ks <- ks.test(clr.counts2[i, 79:160],
                "psn",
                xi = AD_sn_DP$xi2[i],
                omega = AD_sn_DP$omega2[i],
                alpha = AD_sn_DP$alpha2[i])
  AD_KS_test_pvalue[i] <- ks$p.value
}
sum(AD_KS_test_pvalue < 0.05)
# 93/18664 = 0.005
# the skew-normal distribution fits 99.5% of genes well in the AD group


# Fig. (b) AD group
hist(AD_KS_test_pvalue, freq = F, breaks = 20,
     main = NULL, border = "grey83",
     xlab = expression(p))
abline(v = 0.05, lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(b)")



# union of genes which did not fit skew-normal well
ks_pvalue_mat <- cbind.data.frame(control_KS_test_pvalue,
                                  AD_KS_test_pvalue)
colnames(ks_pvalue_mat) <- c("control", "AD")
dim(ks_pvalue_mat[(ks_pvalue_mat$control < 0.05 | ks_pvalue_mat$AD < 0.05), ])
# 323/18664 = 0.017
# overall, 98.3% of the genes for control vs. AD fitted skew-normal well


###END###