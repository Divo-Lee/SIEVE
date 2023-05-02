#' Fit the skew-normal distribution to CLR-transformed RNA-Seq data
#'
#' @description Estimate the mean, standard deviation, and skewness parameters of
#'              the skew-normal distribution using centered log-ratio (CLR)
#'              transformed RNA-Seq data.
#'
#' @param data A table of CLR-transformed count data, where genes/transcripts on the rows and
#'             samples on columns.
#'
#' @return
#' \item{mu}{The maximum likelihood estimate of the mean parameter.}
#' \item{se.mu}{The standard error of the maximum likelihood estimate of \code{mu}.}
#' \item{z.mu}{The Wald statistic for \code{mu}.}
#' \item{p.mu}{The p-value of the Wald statistic for \code{mu}.}
#' \item{sigma}{The maximum likelihood estimate of the standard deviation parameter.}
#' \item{se.sigma}{The standard error of the maximum likelihood estimate of \code{sigma}.}
#' \item{z.sigma}{The Wald statistic for \code{sigma}.}
#' \item{p.sigma}{The p-value of the Wald statistic for \code{sigma}.}
#' \item{gamma}{The maximum likelihood estimate of the skewness parameter.}
#' \item{se.gamma}{The standard error of the maximum likelihood estimate of \code{gamma}.}
#' \item{z.gamma}{The Wald statistic for \code{gamma}.}
#' \item{p.gamma}{The p-value of the Wald statistic for \code{gamma}.}
#'
#'
#' @importFrom sn selm
#'
#' @examples
#'  library(SIEVE)
#'  data(clrCounts1)
#'  clr.SN.fit(clrCounts1[1:2, ])
#'  clr.SN.fit(clrCounts1[1, ])
#'
#'
#' @export
clr.SN.fit <- function(data){
  # data is only for one group
  # data must be clr-transformed counts

  if (is.matrix(data) == FALSE) {
    # only one particular gene
    esti._mat <- c(rep(NA, 12))
    names(esti._mat) <- c("mu", "se.mu", "z.mu", "p.mu",
                         "sigma", "se.sigma", "z.sigma", "p.sigma",
                         "gamma", "se.gamma.", "z.gamma", "p.gamma")
    fit_sn <- selm(data ~ 1,
                   family="SN") # regular MLE without penalty, centered parameters (CP) skew-normal (SN)
    if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
      esti._mat[1:4] <- summary(fit_sn)@param.table[1,]
      esti._mat[5:8] <- summary(fit_sn)@param.table[2,]
      esti._mat[9:12] <- summary(fit_sn)@param.table[3,]
    } else {
      fit_sn <- selm(data ~ 1,
                     method = "MPLE", # "Qpenalty"
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        esti._mat[1:4] <- summary(fit_sn)@param.table[1,]
        esti._mat[5:8] <- summary(fit_sn)@param.table[2,]
        esti._mat[9:12] <- summary(fit_sn)@param.table[3,]
      } else {
        fit_sn <- selm(data ~ 1,
                       method = "MPLE",
                       penalty = "MPpenalty",
                       family="SN")
        esti._mat[1:4] <- summary(fit_sn)@param.table[1,]
        esti._mat[5:8] <- summary(fit_sn)@param.table[2,]
        esti._mat[9:12] <- summary(fit_sn)@param.table[3,]
      }
    }
    return(esti._mat)
  } else {
    d1 <- dim(data)[1]
    # number of genes >= 2
    esti._mat <- matrix(nrow = d1, ncol = 12)
    colnames(esti._mat) <-c("mu", "se.mu", "z.mu", "p.mu",
                            "sigma", "se.sigma", "z.sigma", "p.sigma",
                            "gamma", "se.gamma.", "z.gamma", "p.gamma")
    for (i in 1:d1) {
      fit_sn <- selm(data[i,] ~ 1,
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        esti._mat[i, 1:4] <- summary(fit_sn)@param.table[1,]
        esti._mat[i, 5:8] <- summary(fit_sn)@param.table[2,]
        esti._mat[i, 9:12] <- summary(fit_sn)@param.table[3,]
      } else {
        fit_sn <- selm(data[i,] ~ 1,
                       method = "MPLE",
                       family="SN")
        if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
          esti._mat[i, 1:4] <- summary(fit_sn)@param.table[1,]
          esti._mat[i, 5:8] <- summary(fit_sn)@param.table[2,]
          esti._mat[i, 9:12] <- summary(fit_sn)@param.table[3,]
        } else {
          fit_sn <- selm(data[i,] ~ 1,
                         method = "MPLE",
                         penalty = "MPpenalty",
                         family="SN")
          esti._mat[i, 1:4] <- summary(fit_sn)@param.table[1,]
          esti._mat[i, 5:8] <- summary(fit_sn)@param.table[2,]
          esti._mat[i, 9:12] <- summary(fit_sn)@param.table[3,]
        }
      }
    }
    row.names(esti._mat) <- row.names(data)
    return(esti._mat)
  }
}

