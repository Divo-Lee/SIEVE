#' Differential Expression (DE) test in RNA-Seq data based on skew-normal distribution
#'
#' @description Model CLR-transformed RNA-Seq data using the skew-normal distribution,
#'              and then conduct a statistical test for finding genes/transcripts with
#'              differential expression using the Wald test.
#'
#' @param data CLR-transformed counts matrix for genes in rows and
#'             samples in columns
#' @param group 2 groups (control vs. treatment)
#'
#' @return An object of `\code{clrDE}' class that contains the results of the DV test and
#'         associated information:
#'  \item{DE}{The difference of \code{mu} between
#'         gourp 2 and group 1 (\code{mu2} \code{-} \code{mu1}.}
#'  \item{se}{The standard error of \code{DE}.}
#'  \item{z}{The observed Wald statistic.}
#'  \item{pval}{The unadjusted p-value of the Wald test.}
#'  \item{adj_pval}{The p-value of the Wald test adjusted using the Benjamini-Yekutieli procedure.}
#'  \item{mu1}{The maximum likelihood estimate of the standard deviation parameter for group 1.}
#'  \item{se.mu1}{The standard error of the maximum likelihood estimate of \code{mu1}.}
#'  \item{z.mu1}{The Wald statistic for \code{mu1}.}
#'  \item{p.mu1}{The p-value of the Wald test for \code{mu1}.}
#'  \item{mu2}{The maximum likelihood estimate of the standard deviation parameter for group 2.}
#'  \item{se.mu2}{The standard error of the maximum likelihood estimate of \code{mu2}.}
#'  \item{z.mu2}{The Wald statistic for \code{mu2}.}
#'  \item{p.mu2}{The p-value of the Wald test for \code{mu2}.}
#'
#'
#' @import sn
#'
#' @examples
#'  library(SIEVE)
#'  data(clrCounts3) # The first 50 genes (gene1 to gene50) are DE genes
#'  groups <- c(rep(0, 200), rep(1, 200))
#'  clrDE_test <- clrDE(clrCounts3, group = groups)
#'  sum(is.na(clrDE_test))  # check NA values
#'  head(clrDE_test, 5) # adj_pval < 0.05, DE genes
#'  tail(clrDE_test, 5) # adj_pval > 0.05, non-DE genes
#'
#' @export
clrDE<- function(data = NULL,
                 group = NULL){
  # data = CLR-transformed counts
  # DE test for two groups only

  if (is.factor(group) == F){group = as.factor(group)}
  colnames(data) <- as.factor(group)

  clr_count_group1 <- data[, group == levels(group)[1]]
  clr_count_group2 <- data[, group == levels(group)[2]]
  d2 <- dim(data)[1]

  mean_Matrix_group1 <- matrix(nrow = d2, ncol = 4)
  mean_Matrix_group2 <- matrix(nrow = d2, ncol = 4)

  for (i in 1:d2) {
    fit_sn <- selm(clr_count_group1[i,] ~ 1,
                   family="SN")
    if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
      mean_Matrix_group1[i, ] <- summary(fit_sn)@param.table[1, ]
    } else {
      fit_sn <- selm(clr_count_group1[i,] ~ 1,
                     method = "MPLE",
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        mean_Matrix_group1[i, ] <- summary(fit_sn)@param.table[1, ]
      } else {
        fit_sn <- selm(clr_count_group1[i,] ~ 1,
                       method = "MPLE",
                       penalty = "MPpenalty",
                       family="SN")
        mean_Matrix_group1[i, ] <- summary(fit_sn)@param.table[1, ]
      }
    }
  }
  colnames(mean_Matrix_group1) <- c("mu1", "se.mu1", "z.mu1", "p.mu1")

  for (i in 1:d2) {
    fit_sn <- selm(clr_count_group2[i,] ~ 1,
                   family="SN")
    if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
      mean_Matrix_group2[i, ] <- summary(fit_sn)@param.table[1, ]
    } else {
      fit_sn <- selm(clr_count_group2[i,] ~ 1,
                     method = "MPLE",
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        mean_Matrix_group2[i, ] <- summary(fit_sn)@param.table[1, ]
      } else {
        fit_sn <- selm(clr_count_group2[i,] ~ 1,
                       method = "MPLE",
                       penalty = "MPpenalty",
                       family="SN")
        mean_Matrix_group2[i, ] <- summary(fit_sn)@param.table[1, ]
      }
    }
  }
  colnames(mean_Matrix_group2) <- c("mu2", "se.mu2", "z.mu2", "p.mu2")

  DE_test <- matrix(nrow = d2, ncol = 4)
  colnames(DE_test) <- c("DE", "se" ,"z", "pval")

  for (i in 1:d2) {
    diff <- mean_Matrix_group2[i, 1] - mean_Matrix_group1[i, 1]
    sd2 <- sqrt((mean_Matrix_group1[i, 2])^2 + (mean_Matrix_group2[i, 2])^2)
    DE_test[i, 1] <- diff
    DE_test[i, 2] <- sd2
    DE_test[i, 3] <- diff/sd2
    DE_test[i, 4] <- 2*pnorm(-abs(DE_test[i, 3]), 0, 1)
  }

  adjusted_p <- p.adjust(DE_test[,4], "BY")
  DE_test <- cbind(DE_test, adjusted_p)
  colnames(DE_test)[5] <- "adj_pval"
  result <- cbind(DE_test, mean_Matrix_group1, mean_Matrix_group2)
  row.names(result) <- row.names(data)
  result <- data.frame(result)
  return(result)
}
