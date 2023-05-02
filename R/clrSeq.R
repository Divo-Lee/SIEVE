#' Fit the skew-normal model to CLR-transformed RNA-Seq data for 2 groups
#'
#' @description Estimate the mean, standard deviation, and skewness parameters of
#'              the skew-normal distribution using centered log-ratio (CLR)
#'              transformed RNA-Seq data for 2 groups.
#'
#'
#' @param data A table of CLR-transformed count data, which genes/transcripts on the rows and
#'  samples on the columns for 2 groups.
#' @param group A vector specifying the group labels of the data.
#'
#' @return
#'  \item{mu1}{The maximum likelihood estimate of mean parameter for group 1.}
#'  \item{se.mu1}{The standard error of the maximum likelihood estimate of \code{mu1}.}
#'  \item{z.mu1}{The Wald statistic for \code{mu1}.}
#'  \item{p.mu1}{The p-value of the Wald test for \code{mu1}.}
#'  \item{sigma1}{The maximum likelihood estimate of the standard deviation parameter for group 1.}
#'  \item{se.sigma1}{The standard error of the maximum likelihood estimate of \code{sigma1}.}
#'  \item{z.sigma1}{The Wald statistic for \code{sigma1}.}
#'  \item{p.sigma1}{The p-value of the Wald test for \code{sigma1}.}
#'  \item{gamma1}{The maximum likelihood estimate of the skewness parameter for group 1.}
#'  \item{se.gamma1}{The standard error of the maximum likelihood estimate of \code{gamma1}.}
#'  \item{z.gamma1}{The Wald statistic for \code{gamma1}.}
#'  \item{p.gamma1}{The p-value of the Wald test for \code{gamma1}.}
#'  \item{mu2}{The maximum likelihood estimate of mean parameter for group 2.}
#'  \item{se.mu2}{The standard error of the maximum likelihood estimate of \code{mu2}.}
#'  \item{z.mu2}{The Wald statistic for \code{mu2}.}
#'  \item{p.mu2}{The p-value of the Wald test for \code{mu2}.}
#'  \item{sigma2}{The maximum likelihood estimate of the standard deviation parameter for group 2.}
#'  \item{se.sigma2}{The standard error of the maximum likelihood estimate of \code{sigma2}.}
#'  \item{z.sigma2}{The Wald statistic for \code{sigma2}.}
#'  \item{p.sigma2}{The p-value of the Wald test for \code{sigma2}.}
#'  \item{gamma2}{The maximum likelihood estimate of the skewness parameter for group 2.}
#'  \item{se.gamma2}{The standard error of the maximum likelihood estimate of \code{gamma2}.}
#'  \item{z.gamma2}{The Wald statistic for \code{gamma2}.}
#'  \item{p.gamma2}{The p-value of the Wald test for \code{gamma2}.}
#'
#'
#' @examples
#'    library(SIEVE)
#'    data("clrCounts2")
#'    groups <- c(rep(0, 200), rep(1, 200))
#'    clrSeq(clrCounts2[c(1:5),], group = groups)
#'
#'    clr_Seq <- clrSeq(clrCounts2, group = groups)
#'    tail(clr_Seq, 5)
#'
#' @export
 clrSeq <- function(data = NULL, group = NULL){
  # only for two group
  # data must be clr-transformed counts
  if (is.factor(group) == F){group = as.factor(group)}
  colnames(data) <- as.factor(group)

  clr_count_group1 <- data[, group == levels(group)[1]]
  clr_count_group2 <- data[, group == levels(group)[2]]

  group1_seq <- clr.SN.fit(clr_count_group1)
  group2_seq <- clr.SN.fit(clr_count_group2)

  result <- cbind(group1_seq, group2_seq)
  row.names(result) <- row.names(data)
  colnames(result) <-  c("mu1", "se.mu1", "z.mu1", "p.mu1",
                         "sigma1", "se.sigma1", "z.sigma1", "p.sigma1",
                         "gamma1", "se.gamma1", "z.gamma1", "p.gamma1",
                         "mu2", "se.mu2", "z.mu2", "p.mu2",
                         "sigma2", "se.sigma2", "z.sigma2", "p.sigma2",
                         "gamma2", "se.gamma2", "z.gamma2", "p.gamma2")
  result <- data.frame(result)
  return(result)
}


