# 04 - Build DESeq2 objects for both bulk models
#
# rRNA-biotype genes are removed, genes with fewer than 10 reads across all
# samples are discarded, and the paired structure of the data (independent
# experiment for J1.1, donor for CD4) is modelled as a replicate term.
#
# Input :  Merged_*_Counts_With_Length.csv, Counts.sample_info_*.csv
# Output:  dds_J1.1.rds, dds_CD4.rds, PCA_J1.1.jpg, PCA_CD4T.jpg

library(DESeq2)
library(biomaRt)
library(ggplot2)

out_dir <- "results"

# rRNA genes to exclude
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
rrna_genes <- getBM(
  filters    = "biotype", values = "rRNA",
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  mart       = mart
)

build_dds <- function(counts_file, coldata_file, label, plot_file) {
  counts  <- read.csv(file.path(out_dir, counts_file), header = TRUE, row.names = 1)
  counts  <- counts[, -1]                       # drop the Length column
  coldata <- read.csv(coldata_file, header = TRUE, row.names = 1)

  coldata$treatment <- as.factor(coldata$treatment)
  coldata$replicate <- as.factor(coldata$replicate)

  stopifnot(all(colnames(counts) == rownames(coldata)))

  counts <- counts[!(rownames(counts) %in% rrna_genes$ensembl_gene_id), ]

  dds <- DESeqDataSetFromMatrix(
    countData = counts, colData = coldata,
    design    = ~ 0 + treatment + replicate
  )
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)

  # PCA on variance-stabilised counts
  vsd     <- vst(dds, blind = TRUE)
  pca_dat <- plotPCA(vsd, intgroup = c("treatment", "replicate"), returnData = TRUE)
  pct     <- round(100 * attr(pca_dat, "percentVar"))

  p <- ggplot(pca_dat, aes(PC1, PC2, color = treatment, shape = replicate)) +
    geom_point(size = 3) +
    labs(title = paste("PCA of", label),
         x = paste0("PC1: ", pct[1], "%"),
         y = paste0("PC2: ", pct[2], "%")) +
    theme_minimal()

  ggsave(file.path(out_dir, plot_file), p, width = 8, height = 6, dpi = 300)
  dds
}

dds_j11 <- build_dds("Merged_J1.1_Counts_With_Length.csv",
                     "Counts.sample_info_J1.1.csv", "J1.1 cells", "PCA_J1.1.jpg")
saveRDS(dds_j11, file.path(out_dir, "dds_J1.1.rds"))

dds_cd4 <- build_dds("Merged_CD4_Counts_With_Length.csv",
                     "Counts.sample_info_CD4.csv", "CD4 T cells", "PCA_CD4T.jpg")
saveRDS(dds_cd4, file.path(out_dir, "dds_CD4.rds"))
