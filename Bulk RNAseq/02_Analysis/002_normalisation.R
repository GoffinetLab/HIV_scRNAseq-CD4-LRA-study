# 02 - TPM and RPKM normalisation
#
# Normalised matrices are descriptive output only. All differential expression
# testing is performed on raw counts (script 04).
#
# Input :  Merged_J1.1_Counts_With_Length.csv, Merged_CD4_Counts_With_Length.csv
# Output:  Counts.TPM.matrix.csv, Counts.RPKM.matrix.csv

library(tidyverse)

out_dir <- "results"

j11 <- read.csv(file.path(out_dir, "Merged_J1.1_Counts_With_Length.csv"))
cd4 <- read.csv(file.path(out_dir, "Merged_CD4_Counts_With_Length.csv"))

total        <- merge(j11, cd4, by = c("Geneid", "Length"))
gene_info    <- total[, c("Geneid", "Length")]
raw_counts   <- total[, setdiff(colnames(total), c("Geneid", "Length"))]
gene_lengths <- total$Length

# TPM: length-normalise first, then scale by the per-sample total
counts_per_kb <- raw_counts / (gene_lengths / 1000)
tpm <- t(t(counts_per_kb) / colSums(counts_per_kb)) * 1e6
write.csv(cbind(gene_info, tpm), file.path(out_dir, "Counts.TPM.matrix.csv"), row.names = FALSE)

# RPKM: scale by library size first, then length-normalise
rpkm <- sweep(raw_counts, 2, colSums(raw_counts) / 1e6, "/") / (gene_lengths / 1000)
write.csv(cbind(gene_info, rpkm), file.path(out_dir, "Counts.RPKM.matrix.csv"), row.names = FALSE)
