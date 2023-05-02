#' Sieve Differential Expression, Variability and Skewness Genes
#' for RNA-Seq Data Based on Skew-normal Distribution
#'
#' @description Model CLR-transformed RNA-Seq data using the skew-normal
#'              distribution, and then conduct statistical tests for finding
#'              genes/transcripts with differential expression, variability
#'              and skewness using the Wald test.
#'
#' @param clrSeq_result The result of clrSeq() function.
#' @param alpha_level The adjusted p-value significance cutoff used for flagging
#'                    the differential expression, variability and skewness genes.
#' @param order_DE Logical string, "FALSE" for no ordering, "TRUE" for ordering by the value of \code{DE}.
#' @param order_LFC Logical string, "FALSE" for no ordering, "TRUE" for ordering by the value of \code{LFC}.
#' @param order_DS Logical string, "FALSE" for no ordering, "TRUE" for ordering by the value of \code{DS}.
#' @param order_sieve Character/logical string specifying the order method. Possibilities are "DE"
#'                    the value of \code{DE}, "LFC" for log2-fold-change of variability,
#'                    "DV" for the value of \code{DV} or "FALSE" for no ordering.
#'
#' @return Lists of \code{clrDE_test}, \code{clrDV_test}, \code{clrDS_test} and
#'         \code{clrSIEVE_tests} classes that contain the results of the DE, DV,
#'         DS tests and one-stop those three tests, and associated information:
#'  \item{clrDE_test}{A data.frame contating the results of differential expression test.}
#'  \item{clrDV_test}{A data.frame contating the results of differential variability test.}
#'  \item{clrDS_test}{A data.frame contating the results of differential skewness test.}
#'  \item{clrSIEVE_tests}{A data.frame contating the results of differential expression, variability and skewness tests.}
#'  \item{DE}{The difference of mean (\code{mu}) between group 2 and group 1 (\code{mu2} \code{-} \code{mu1}).}
#'  \item{se_DE}{The standard error of \code{DE}.}
#'  \item{z_DE}{The observed Wald statistic of \code{DE}.}
#'  \item{pval_DE}{The unadjusted p-value of Wald test of \code{DE}.}
#'  \item{adj_pval_DE}{The p-value of the Wald test of  adjusted using the Benjamini-Yekutieli procedure.}
#'  \item{mu1}{The maximum likelihood estimate of the mean parameter for group 1.}
#'  \item{mu2}{The maximum likelihood estimate of the mean parameter for group 2.}
#'  \item{de_indicator}{1: DE gene, 0: non-DE gene.}
#'  \item{SD_ratio}{The ratio of standard deviation (\code{sigma}) between group 2 and group 1 (\code{sigma2}\code{/}\code{sigma1}.). }
#'  \item{LFC}{\code{log2} fold change: \code{log2}(\code{sigma2}\code{/}\code{sigma1}).}
#'  \item{DV}{The difference of standard deviation (\code{sigma}) between group 2 and group 1 (\code{sigma2} \code{-} \code{sigma1}.).}
#'  \item{se_DV}{The standard error of \code{DV}.}
#'  \item{z_DV}{The observed Wald statistic of \code{DV}.}
#'  \item{pval_DV}{The unadjusted p-value of Wald test of \code{DV}.}
#'  \item{adj_pval_DV}{he p-value of the Wald test of \code{DV} adjusted using the Benjamini-Yekutieli procedure.}
#'  \item{sigma1}{The maximum likelihood estimate of the standard deviation parameter for group 1.}
#'  \item{sigma2}{The maximum likelihood estimate of the standard deviation parameter for group 2.}
#'  \item{dv_indicator}{1: DV gene, 0: non-DV gene.}
#'  \item{DS}{The difference of skewness parameter (\code{gamma}) between group 2 and group 1 (\code{gamma2} \code{-} \code{gamma1}.).}
#'  \item{se_DS}{The standard error of \code{DS}.}
#'  \item{z_DS}{The observed Wald statistic of \code{DS}.}
#'  \item{pval_DS}{The unadjusted p-value of Wald test of \code{DS}.}
#'  \item{adj_pval_DS}{The p-value of the Wald test of \code{DS} adjusted using the Benjamini-Yekutieli procedure.}
#'  \item{gamma1}{The maxiimum likelihood estimate of the skewness parameter for group 1.}
#'  \item{gamma2}{The maxiimum likelihood estimate of the skewness parameter for group 2.}
#'  \item{ds_indicator}{1: DS gene, 0: non-DS gene.}
#'
#'
#' @examples
#' library(SIEVE)
#' data(clrCounts3) # first 50 genes (gene1 to gene50) are DE genes
#' groups <- c(rep(0, 200), rep(1, 200))
#' clrSeq_result3 <- clrSeq(clrCounts3, group = groups) # DE dataset
#' clrSIEVE_result3 <- clrSIEVE(clrSeq_result = clrSeq_result3,
#'                              alpha_level = 0.05,
#'                              order_DE = FALSE,
#'                              order_LFC = FALSE,
#'                              order_DS = FALSE,
#'                              order_sieve = FALSE)
#' clrDE_test3 <- clrSIEVE_result3$clrDE_test # DE test
#' head(clrDE_test3, 5)
#' clrDS_test3 <- clrSIEVE_result3$clrDS_test # DS test
#' clrDS_test3[clrDS_test3$adj_pval_DS < 0.05, ]
#' clrSIEVE_tests3 <- clrSIEVE_result3$clrSIEVE_tests # Sieve DE, DV and DS genes
#' head(clrSIEVE_tests3, 5)
#' tail(clrSIEVE_tests3, 5)
#'
#' @export
#'
clrSIEVE <- function(clrSeq_result = NULL,
                     alpha_level = 0.05,
                     order_DE = FALSE,
                     order_LFC = FALSE,
                     order_DS = FALSE,
                     order_sieve = FALSE){
  # clrSeq_result = the result of clrSeq() function in clrDV package
  # DE, DV and DS tests for two groups only
  # order_sieve = "DE", "LFC", "DS" or FALSE

  if(alpha_level >=1 | alpha_level <= 0){
    print("Wrong input alpha_level!
           Default alpha_level (0.05) will be used.")
    alpha_level <- 0.05
  }

  ## DE test
  clrde_test <- data.frame(DE = clrSeq_result$mu2 - clrSeq_result$mu1,
                           se_DE = sqrt((clrSeq_result$se.mu2)^2 + (clrSeq_result$se.mu1)^2))
  clrde_test$z_DE <- clrde_test$DE/clrde_test$se_DE
  clrde_test$pval_DE <- 2*pnorm(-abs(clrde_test$z_DE), 0, 1)
  clrde_test$adj_pval_DE <- p.adjust(clrde_test$pval_DE, "BY")
  clrde_test$mu1 <- clrSeq_result$mu1
  clrde_test$mu2 <- clrSeq_result$mu2
  clrde_test$de_indicator <- as.integer(clrde_test$adj_pval_DE < alpha_level)
  row.names(clrde_test) <- row.names(clrSeq_result)

  ## DV test
  clrdv_test <- data.frame(SD_ratio = clrSeq_result$sigma2/clrSeq_result$sigma1,
                           LFC = log2(clrSeq_result$sigma2/clrSeq_result$sigma1),
                           DV = clrSeq_result$sigma2 - clrSeq_result$sigma1,
                           se_DV = sqrt((clrSeq_result$se.sigma2)^2 + (clrSeq_result$se.sigma1)^2))
  clrdv_test$z_DV <- clrdv_test$DV/clrdv_test$se_DV
  clrdv_test$pval_DV <- 2*pnorm(-abs(clrdv_test$z_DV), 0, 1)
  clrdv_test$adj_pval_DV <- p.adjust(clrdv_test$pval_DV, "BY")
  clrdv_test$sigma1 <- clrSeq_result$sigma1
  clrdv_test$sigma2 <- clrSeq_result$sigma2
  clrdv_test$dv_indicator <- as.integer(clrdv_test$adj_pval_DV < alpha_level)
  row.names(clrdv_test) <- row.names(clrSeq_result)

  ## DS
  clrds_test <- data.frame(DS = clrSeq_result$gamma2 - clrSeq_result$gamma1,
                           se_DS = sqrt((clrSeq_result$se.gamma2)^2 + (clrSeq_result$se.gamma1)^2))
  clrds_test$z_DS <- clrds_test$DS/clrds_test$se_DS
  clrds_test$pval_DS <- 2*pnorm(-abs(clrds_test$z_DS), 0, 1)
  clrds_test$adj_pval_DS <- p.adjust(clrds_test$pval_DS, "BY")
  clrds_test$gamma1 <- clrSeq_result$gamma1
  clrds_test$gamma2 <- clrSeq_result$gamma2
  clrds_test$ds_indicator <- as.integer(clrds_test$adj_pval_DS < alpha_level)
  row.names(clrds_test) <- row.names(clrSeq_result)


  ## diff_genes_SIEVE
  clrSIEVE_test <- data.frame(DE = clrde_test$DE,
                              adj_pval_DE= clrde_test$adj_pval_DE,
                              SD_ratio = clrdv_test$SD_ratio,
                              LFC = clrdv_test$LFC,
                              DV = clrdv_test$DV,
                              adj_pval_DV = clrdv_test$adj_pval_DV,
                              DS = clrds_test$DS,
                              adj_pval_DS = clrds_test$adj_pval_DS,
                              de_indicator = clrde_test$de_indicator,
                              dv_indicator = clrdv_test$dv_indicator,
                              ds_indicator = clrds_test$ds_indicator)
  row.names(clrSIEVE_test) <- row.names(clrSeq_result)


  if(order_sieve == "LFC" ){
    clrSIEVE_test <- clrSIEVE_test[order(clrSIEVE_test$LFC), ]
  }
  else if(order_sieve == "DS"){
    clrSIEVE_test <- clrSIEVE_test[order(clrSIEVE_test$DS), ]
  }
  else if(order_sieve == "DE"){
    clrSIEVE_test <- clrSIEVE_test[order(clrSIEVE_test$DE), ]
  }
  else if(order_sieve == FALSE){
    clrSIEVE_test <- clrSIEVE_test
  }


  if(order_DE == TRUE){
    clrde_test <- clrde_test[order(clrde_test$DE), ]
  }
  if(order_LFC == TRUE){
    clrdv_test <- clrdv_test[order(clrdv_test$LFC),]
  }
  if(order_DS == TRUE){
    clrds_test <- clrds_test[order(clrds_test$DS), ]
  }


  return(list("clrDE_test" = clrde_test,
              "clrDV_test" = clrdv_test,
              "clrDS_test" = clrds_test,
              "clrSIEVE_tests" = clrSIEVE_test))
}
