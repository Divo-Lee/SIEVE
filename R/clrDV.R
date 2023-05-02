#' Differential variability (DV) test of CLR-transformed RNA-Seq data
#'
#' @description Model CLR-transformed RNA-Seq data using the skew-normal distribution,
#'              and then conduct a statistical test for finding genes/transcripts with
#'              differential variability using the Wald test.
#'
#' @param data A CLR-transformed count matrix with genes/transcripts on the rows and
#'             samples on the columns.
#' @param group A vector specifying the group labels of the data.
#'
#'
#' @return An object of `\code{clrDV}' class that contains the results of the DV test and
#'         associated information:
#'  \item{DV}{The difference of standard deviation (\code{sigma}) between
#'         group 2 and group 1 (\code{sigma2} \code{-} \code{sigma1}).}
#'  \item{se}{The standard error of \code{DV}.}
#'  \item{z}{The observed Wald statistic.}
#'  \item{pval}{The unadjusted p-value of Wald test.}
#'  \item{adj_pval}{The p-value of the Wald test adjusted using the Benjamini-Yekutieli procedure.}
#'  \item{sigma1}{The maximum likelihood estimate of the standard deviation parameter for group 1.}
#'  \item{se.sigma1}{The standard error of the maximum likelihood estimate of \code{sigma1}.}
#'  \item{z.sigma1}{The Wald statistic for \code{sigma1}.}
#'  \item{p.sigma1}{The p-value of the Wald test for \code{sigma1}.}
#'  \item{sigma2}{The maximum likelihood estimate of the standard deviation parameter for group 2.}
#'  \item{se.sigma2}{The standard error of the maximum likelihood estimate of \code{gamma2}.}
#'  \item{z.sigma2}{The Wald statistic for \code{sigma2}.}
#'  \item{p.sigma2}{The p-value of the Wald test for \code{sigma2}.}
#'
#'
#' @importFrom  stats p.adjust
#' @importFrom  stats pnorm
#' @importFrom  sn selm
#'
#'
#' @examples
#'    library(SIEVE)
#'    data(clrCounts2) # first 50 genes (gene1 to gene50) are DV genes
#'    groups <- c(rep(0, 200), rep(1, 200))
#'    clrDV_test <- clrDV(clrCounts2, group = groups)
#'    sum(is.na(clrDV_test))  # check NA values
#'    head(clrDV_test, 5)  # adj_pval < 0.05, DV genes
#'    tail(clrDV_test, 5)  # adj_pval > 0.05, non-DV genes
#'
#'
#' @export
clrDV <- function(data = NULL,
                  group = NULL){
  # data = clr-transformed count table
  # DV test for two groups only

  if (is.factor(group) == F){group = as.factor(group)}
  colnames(data) <- as.factor(group)

  clr_count_group1 <- data[, group == levels(group)[1]]
  clr_count_group2 <- data[, group == levels(group)[2]]
  d2 <- dim(data)[1]

  s.d._Matrix_group1 <- matrix(nrow = d2, ncol = 4)
  s.d._Matrix_group2 <- matrix(nrow = d2, ncol = 4)

  for (i in 1:d2) {
    fit_sn <- selm(clr_count_group1[i,] ~ 1,
                   family="SN")
    if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){ # check NA value
      s.d._Matrix_group1[i, ] <- summary(fit_sn)@param.table[2, ]
    } else {
      fit_sn <- selm(clr_count_group1[i,] ~ 1,
                     method = "MPLE", # Qpenalty (default)
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        s.d._Matrix_group1[i, ] <- summary(fit_sn)@param.table[2, ]
      } else {
        fit_sn <- selm(clr_count_group1[i,] ~ 1,
                       method = "MPLE",
                       penalty = "MPpenalty",
                       family="SN")
        s.d._Matrix_group1[i, ] <- summary(fit_sn)@param.table[2, ]
      }
    }
  }
  colnames(s.d._Matrix_group1) <- c("sigma1", "se.sigma1", "z.sigma1", "p.sigma1")

  for (i in 1:d2) {
    fit_sn <- selm(clr_count_group2[i,] ~ 1,
                   family="SN")
    if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
      s.d._Matrix_group2[i, ] <- summary(fit_sn)@param.table[2, ]
    } else {
      fit_sn <- selm(clr_count_group2[i,] ~ 1,
                     method = "MPLE",
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        s.d._Matrix_group2[i, ] <- summary(fit_sn)@param.table[2, ]
      } else {
        fit_sn <- selm(clr_count_group2[i,] ~ 1,
                       method = "MPLE",
                       penalty = "MPpenalty",
                       family="SN")
        s.d._Matrix_group2[i, ] <- summary(fit_sn)@param.table[2, ]
      }
    }
  }
  colnames(s.d._Matrix_group2) <- c("sigma2", "se.sigma2", "z.sigma2", "p.sigma2")


  s.d._test <- matrix(nrow = d2, ncol = 4)
  colnames(s.d._test) <- c("DV", "se" ,"z", "pval")

  for (i in 1:d2) {
    diff <- s.d._Matrix_group2[i, 1] - s.d._Matrix_group1[i, 1]
    sd2 <- sqrt((s.d._Matrix_group1[i, 2])^2 + (s.d._Matrix_group2[i, 2])^2)
    s.d._test[i, 1] <- diff
    s.d._test[i, 2] <- sd2
    s.d._test[i, 3] <- diff/sd2
    s.d._test[i, 4] <- 2*pnorm(-abs(s.d._test[i, 3]), 0, 1)
  }

  adjusted_p <- p.adjust(s.d._test[, 4], "BY")
  s.d._test <- cbind(s.d._test, adjusted_p)
  colnames(s.d._test)[5] <- "adj_pval"
  result <- cbind(s.d._test, s.d._Matrix_group1, s.d._Matrix_group2)
  row.names(result) <- row.names(data)
  result <- data.frame(result)
  return(result)
}
