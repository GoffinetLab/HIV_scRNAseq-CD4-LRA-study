# 01 - Assemble gene count matrices from per-sample featureCounts output
#
# Input :  counts/*.clean.counts.txt   (Geneid, Length, count)
# Output:  Merged_J1.1_Counts_With_Length.csv
#          Merged_CD4_Counts_With_Length.csv
#          Merged_raw_counts.csv

library(tidyverse)

counts_dir <- "counts"
out_dir    <- "results"
dir.create(out_dir, showWarnings = FALSE)

# Read all per-sample count tables
files      <- list.files(counts_dir, pattern = "\\.clean\\.counts\\.txt$", full.names = TRUE)
count_list <- lapply(files, read.table, header = TRUE)
names(count_list) <- basename(files)

# Split by model: J1.1 samples are prefixed "Run", CD4 donors "D1"-"D3"
j11_names <- grep("^Run",   names(count_list), value = TRUE)
cd4_names <- grep("^D[123]", names(count_list), value = TRUE)

# Merge the count column of each sample on Geneid, then attach gene length
merge_counts <- function(sample_names) {
  cols <- lapply(sample_names, function(fname) {
    df <- count_list[[fname]][, c(1, 3)]
    colnames(df)[2] <- sub("\\.clean\\.counts\\.txt$", "", fname)
    df
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), cols)
  length_info <- count_list[[1]][, c("Geneid", "Length")]
  merge(length_info, merged, by = "Geneid")
}

j11_matrix <- merge_counts(j11_names)
cd4_matrix <- merge_counts(cd4_names)

write.csv(j11_matrix, file.path(out_dir, "Merged_J1.1_Counts_With_Length.csv"), row.names = FALSE)
write.csv(cd4_matrix, file.path(out_dir, "Merged_CD4_Counts_With_Length.csv"), row.names = FALSE)

# Combined matrix across both models
write.csv(
  merge(j11_matrix, cd4_matrix, by = c("Geneid", "Length"), all = TRUE),
  file.path(out_dir, "Merged_raw_counts.csv"), row.names = FALSE
)
