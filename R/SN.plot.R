#' Distribution of CLR-transformed counts
#'
#' @param data CLR-transformed counts for a particular gene/transcript.
#'
#' @description Produce a histogram of observed CLR-transformed counts, with the
#'          fitted skew-normal probability density function for a particular gene/transcript.
#'
#'
#' @import sn
#'
#' @examples
#'    library(SIEVE)
#'    data("clrCounts1")
#'    SN.plot(clrCounts1[1,])
#'    SN.plot(clrCounts1[2,])
#'    SN.plot(clrCounts1[3,])
#'
#' @export
SN.plot <- function(data){
  # data is the CLR-transformed counts of particular gene for 1 group
  # data is a vector
  fit_sn <-  selm(data ~ 1, family="SN")
  plot(fit_sn, which=2, caption = NULL)
}
