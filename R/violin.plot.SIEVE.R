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
#' @import ggplot2 
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
violin.plot.SIEVE <- function(data = NULL,
                               name.gene = NULL,
                               group = NULL,
                               group.names = NULL,
                               xlab = "CLR-transformed count",
                               ylab = "Condition"){
  # Ensure group is a factor with correct labels
  if (!is.factor(group)) {
    group <- as.factor(group)
  }

  # Extract expression values for the specified gene
  gene_values <- as.numeric(data[name.gene, ])
  df <- data.frame(
    gene_values,
    factor(group, levels = levels(group), labels = group.names)
  )

  # Define a custom color palette
  custom_colors <- c("#66c2a5", "#fc8d62", "#8da0cb")[seq_along(levels(df$Group))]

  # Create the ggplot
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group))
  p +  geom_violin(color = "black",
                   alpha = 0.6) +  # filled violins with black border
    geom_jitter(width = 0.2, size = 1, color = "black", shape = 1) +  # add points
    coord_flip() +  # horizontal violins
    scale_fill_manual(values = custom_colors) +
    xlab(ylab) +
    ylab(xlab) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # remove grid
      axis.line = element_line(color = "black"),  # show axes
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none"  # hide legend since x-axis already labels groups
    )
}
