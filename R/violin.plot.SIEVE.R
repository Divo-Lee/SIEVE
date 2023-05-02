#' Violin plots
#'
#' @param data A CLR-transformed count table.
#' @param name.gene Gene/transcript name.
#' @param group  A vector specifying the group labels of the data.
#' @param group.names A vector specifying the group names.
#' @param xlab Name of the x-axis.
#' @param ylab Name of the y-axis.
#'
#' @description Produce violin plots of CLR-transformed count data for two or three groups.
#'
#' @import vioplot
#' @import grDevices
#' @import graphics
#'
#'
#' @examples
#'   library(SIEVE)
#'   data(clrCounts2)  # first 50 genes (gene1 to gene50) are DV genes
#'   data(clrCounts3)  # first 50 genes (gene1 to gene50) are DE genes
#'   group0 <- c(rep(0, 200), rep(1, 200))
#'   group1 <- c(rep(0, 200), rep(1, 100), rep(2, 100))
#'   violin.plot.SIEVE(data = clrCounts2, "gene1", group = group0,
#'                 group.names = c("control", "case")) # DV
#'   violin.plot.SIEVE(data = clrCounts2, "gene1", group = group1,
#'                 group.names = c("control", "case1", "case2")) # DV
#'   violin.plot.SIEVE(data = clrCounts3, "gene1", group = group0,
#'                 group.names = c("control", "case")) # DE
#'   violin.plot.SIEVE(data = clrCounts3, "gene2", group = group0,
#'                 group.names = c("control", "case")) # DE
#'   violin.plot.SIEVE(data = clrCounts3, "gene200", group = group0,
#'                 group.names = c("control", "case")) # non-DE
#' @export
violin.plot.SIEVE <- function(data = NULL, name.gene = NULL,
                          group = NULL, group.names = NULL,
                          xlab="CLR-transformed count",
                          ylab="Condition"){
  # data is the CLR-transformed count table
  # 2 or 3 groups only
  if (is.factor(group) == F){group = as.factor(group)}

  if (length(group.names) == 2){
    clr_count_group1 <- data[name.gene, ][group == levels(group)[1]]
    clr_count_group2 <- data[name.gene, ][group == levels(group)[2]]

    vioplot(clr_count_group1,
            clr_count_group2,
            names=group.names,  pchMed="",
            col = c(0,0), border ="black", horizontal = T,
            rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
            xlab=xlab, ylab=ylab)

    stripchart(list(clr_count_group1,
                    clr_count_group2),
               method="jitter", vertical=F,  add=TRUE,
               pch=1, cex=0.5, col="black")


  } else if (length(group.names) == 3){
    clr_count_group1 <- data[name.gene, ][group == levels(group)[1]]
    clr_count_group2 <- data[name.gene, ][group == levels(group)[2]]
    clr_count_group3 <- data[name.gene, ][group == levels(group)[3]]

    vioplot(clr_count_group1,
            clr_count_group2,
            clr_count_group3,
            names=group.names,  pchMed="",
            col = c(0,0), border ="black", horizontal = T,
            rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
            xlab=xlab, ylab=ylab)

    stripchart(list(clr_count_group1,
                    clr_count_group2,
                    clr_count_group3),
                    method="jitter", vertical=F,  add=TRUE,
                    pch=1, cex=0.5, col="black")

  }
}
