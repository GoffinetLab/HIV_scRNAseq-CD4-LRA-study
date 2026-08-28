# 03 - Annotate Ensembl gene identifiers
#
# Queries the Ensembl GRCh38.p14 October 2024 archive in batches, falling back to
# the current release and then GRCh37 for identifiers that remain unresolved.
#
# Input :  Counts.RPKM.matrix.csv
# Output:  ensembl_annotations.csv, Counts.RPKM.Annotated.csv

library(biomaRt)
library(dplyr)

out_dir     <- "results"
batch_size  <- 500

rpkm_matrix <- read.csv(file.path(out_dir, "Counts.RPKM.matrix.csv"))
gene_ids    <- unique(rpkm_matrix$Geneid)

connect_mart <- function(host) {
  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = host)
}

fetch <- function(ids, mart) {
  tryCatch(
    getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol",
                     "chromosome_name", "gene_biotype", "description"),
      filters = "ensembl_gene_id", values = ids, mart = mart
    ) %>% mutate(across(everything(), as.character)),
    error = function(e) data.frame()
  )
}

hosts <- c("https://oct2024.archive.ensembl.org",   # GRCh38.p14
           "https://www.ensembl.org",               # current release
           "https://grch37.ensembl.org")            # GRCh37

batches <- split(gene_ids, ceiling(seq_along(gene_ids) / batch_size))

annotations <- bind_rows(lapply(batches, function(batch) {
  result  <- data.frame()
  pending <- batch
  for (host in hosts) {
    if (length(pending) == 0) break
    found   <- fetch(pending, connect_mart(host))
    result  <- bind_rows(result, found)
    pending <- setdiff(pending, found$ensembl_gene_id)
  }
  result
})) %>% distinct(ensembl_gene_id, .keep_all = TRUE)

write.csv(annotations, file.path(out_dir, "ensembl_annotations.csv"), row.names = FALSE)

rpkm_annotated <- rpkm_matrix %>%
  left_join(annotations, by = c("Geneid" = "ensembl_gene_id")) %>%
  relocate(hgnc_symbol, chromosome_name, gene_biotype, description, .after = Geneid)

write.csv(rpkm_annotated, file.path(out_dir, "Counts.RPKM.Annotated.csv"), row.names = FALSE)
