###############################################################
#SIEVE: One-stop differential expression, variability,
# and skewness using RNA-Seq data
#Authors: Hongxiang Li and Tsung Fei Khang
#Email: hxli@ynnu.edu.cn
#Last update: 13 Feb. 2026
#R Codes for comparison of ORA and GSEA 
################################################################

## Required R packages
#devtools::install_github("alserglab/fgsea")
#BiocManager::install("org.Hs.eg.db")

library(fgsea)
library(org.Hs.eg.db)
library(dplyr)

###########################################################################
### 1. Calculate Composite Ranking Score (Absolute Max Selection)
###########################################################################

# Extract Wald statistics from clrSIEVE results
DE_AD_table <- clrSIEVE_result$clrDE_test
DV_AD_table <- clrSIEVE_result$clrDV_test
DS_AD_table <- clrSIEVE_result$clrDS_test

# Build a matrix of Z-scores (Wald statistics)
DE_DV_DS_Wald_Stat_Mat <- cbind(
  DE_AD_table$z_DE,
  DV_AD_table$z_DV,
  DS_AD_table$z_DS
)
row.names(DE_DV_DS_Wald_Stat_Mat) <- rownames(DE_AD_table)

# Identify the index of the column with the maximum absolute value for each row
idx <- apply(abs(DE_DV_DS_Wald_Stat_Mat), 1, which.max)

# Extract the original signed values corresponding to the absolute maximums
max_vec <- DE_DV_DS_Wald_Stat_Mat[
  cbind(seq_len(nrow(DE_DV_DS_Wald_Stat_Mat)), idx)
]
names(max_vec) <- rownames(DE_DV_DS_Wald_Stat_Mat)

# Rank the statistics in descending order
max_vec_ranked <- sort(max_vec, decreasing = TRUE)

# Filter out duplicate gene symbols to ensure unique gene-ranking
# Keeps only genes that appear exactly once in the dataset
max_vec_ranked <- max_vec_ranked[
  names(max_vec_ranked) %in% 
    names(which(table(names(max_vec_ranked)) == 1))
]

###########################################################################
### 2. Define Reference GO Terms (Top 20 ORA-defined)
###########################################################################

# Reference list of GO IDs identified via ORA
go_ids_top_20 <- c(
  "GO:1902531", "GO:0090150", "GO:0072657", "GO:0072599",
  "GO:0072359", "GO:0070972", "GO:0048514", "GO:0045047",
  "GO:0043062", "GO:0030198", "GO:0023056", "GO:0022610",
  "GO:0010647", "GO:0009967", "GO:0007155", "GO:0006614",
  "GO:0006613", "GO:0006612", "GO:0001568", "GO:0000184"
) 

# Descriptions for formatting Supplementary Table S14
go_descriptions <- c(
  "GO:1902531" = "regulation of intracellular signal transduction",
  "GO:0090150" = "establishment of protein localization to membrane",
  "GO:0072657" = "protein localization to membrane",
  "GO:0072599" = "establishment of protein localization to ER",
  "GO:0072359" = "circulatory system development",
  "GO:0070972" = "protein localization to ER",
  "GO:0048514" = "blood vessel morphogenesis",
  "GO:0045047" = "protein targeting to ER",
  "GO:0043062" = "extracellular structure organization",
  "GO:0030198" = "extracellular matrix organization",
  "GO:0023056" = "positive regulation of signaling",
  "GO:0022610" = "biological adhesion",
  "GO:0010647" = "positive regulation of cell communication",
  "GO:0009967" = "positive regulation of signal transduction",
  "GO:0007155" = "cell adhesion",
  "GO:0006614" = "SRP-dependent cotranslational protein targeting",
  "GO:0006613" = "cotranslational protein targeting to membrane",
  "GO:0006612" = "protein targeting to membrane",
  "GO:0001568" = "blood vessel development",
  "GO:0000184" = "nuclear-transcribed mRNA catabolic process"
) 

###########################################################################
### 3. Build Gene Sets (Mapping GO to ENSEMBL)
###########################################################################

# Primary mapping: GO to ENSEMBL
go_map_main <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_ids_top_20,
  keytype = "GOALL",
  columns = c("GOALL", "ENSEMBL")
) %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::distinct()

# Split into named list for fgsea
gene_sets_table <- split(go_map_main$ENSEMBL, go_map_main$GOALL)

# Remedial mapping for missing IDs using secondary environments
missing_ids <- setdiff(go_ids_top_20, names(gene_sets_table))
if (length(missing_ids) > 0) {
  for (go_id in missing_ids) {
    entrez_vec <- tryCatch({
      obj <- mget(go_id, envir = org.Hs.egGO2ALLEGS, ifnotfound = NA)[[1]]
      if (all(is.na(obj))) character(0) else as.character(obj)
    }, error = function(e) character(0))
    
    if (length(entrez_vec) == 0) next
    
    ens_map <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = entrez_vec,
      column = "ENSEMBL", keytype = "ENTREZID", multiVals = "list"
    )
    
    ens_vec <- unique(unlist(ens_map))
    ens_vec <- ens_vec[!is.na(ens_vec)]
    
    if (length(ens_vec) > 0) {
      gene_sets_table[[go_id]] <- ens_vec
    }
  }
}

###########################################################################
### 4. GSEA Calculation and Results Formatting
###########################################################################

set.seed(100)
GSEA_DE_DV_DS <- fgseaSimple(
  pathways = gene_sets_table,
  stats = max_vec_ranked,
  nperm = 5*10^5,
  minSize = 5,
  maxSize = 2000,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL
)

# Format the final data frame for Supplementary Table S14
final_output <- as.data.frame(GSEA_DE_DV_DS) %>%
  dplyr::filter(pathway %in% names(go_descriptions)) %>%
  dplyr::mutate(
    pathway = as.character(pathway),
    description = go_descriptions[pathway]
  ) %>%
  dplyr::select(
    GO_ID = pathway,
    description,
    NES,
    padj
  ) %>%
  dplyr::arrange(padj)

# Display final table
print(final_output)


### END ###

