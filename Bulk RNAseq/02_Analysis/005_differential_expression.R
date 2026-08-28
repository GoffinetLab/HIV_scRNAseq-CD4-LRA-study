# 05 - Differential expression for all contrasts and all three models
#
# Log2 fold changes are shrunken with ashr; p-values are adjusted by the
# Benjamini-Hochberg method. Single-cell results are pseudobulk profiles
# supplied as a list of data frames, reformatted here to share column names
# with the DESeq2 output so the three models can be compared.
#
# Input :  dds_J1.1.rds, dds_CD4.rds, deg.list_donor.treat.rds
# Output:  DEGs_vs_DMSO_all_models.xlsx, DEGs_extra.xlsx

library(DESeq2)
library(openxlsx)
library(org.Hs.eg.db)

out_dir <- "results"

dds_j11 <- readRDS(file.path(out_dir, "dds_J1.1.rds"))
dds_cd4 <- readRDS(file.path(out_dir, "dds_CD4.rds"))
sc_res  <- readRDS("deg.list_donor.treat.rds")

res_list <- list()

# Shrink, annotate and store one DESeq2 contrast
process_res <- function(dds, res, name) {
  df <- as.data.frame(lfcShrink(dds, res = res, type = "ashr"))
  df$ensID <- rownames(df)
  rownames(df) <- NULL

  df$hgnc_symbol <- mapIds(org.Hs.eg.db, keys = df$ensID,
                           column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  # fall back to the Ensembl identifier where no symbol exists
  df$hgnc_symbol <- ifelse(is.na(df$hgnc_symbol), df$ensID, df$hgnc_symbol)

  res_list[[name]] <<- df[, c("ensID", "hgnc_symbol",
                              setdiff(colnames(df), c("ensID", "hgnc_symbol")))]
}

# Pseudobulk single-cell results, keyed on gene symbol
process_sc <- function(df, name) {
  df <- as.data.frame(df)
  df$hgnc_symbol <- rownames(df)
  rownames(df) <- NULL
  res_list[[name]] <<- df[, c("hgnc_symbol", setdiff(colnames(df), "hgnc_symbol"))]
}

# Each treatment against the DMSO vehicle control
contrasts_vs_dmso <- list(
  c("treatmentPano",      "j1.1_pano_vs_dmso"),
  c("treatmentIFN",       "j1.1_ifn_vs_dmso"),
  c("treatmentPano_IFN",  "j1.1_panoifn_vs_dmso")
)
for (x in contrasts_vs_dmso) {
  process_res(dds_j11, results(dds_j11, contrast = list(x[1], "treatmentDMSO")), x[2])
  process_res(dds_cd4, results(dds_cd4, contrast = list(x[1], "treatmentDMSO")),
              sub("^j1\\.1", "cd4", x[2]))
}

process_sc(sc_res$Panobinostat_v_Mock, "scCD4_panobinostat_vs_mock")
process_sc(sc_res$Vorinostat_v_Mock,   "scCD4_vorinostat_vs_mock")
process_sc(sc_res$`IL2/PHA_v_Mock`,    "scCD4_IL2PHA_vs_mock")

write_sheets <- function(names, sheets, file) {
  wb <- createWorkbook()
  for (i in seq_along(names)) {
    addWorksheet(wb, sheets[i])
    writeData(wb, sheets[i], res_list[[names[i]]])
  }
  saveWorkbook(wb, file.path(out_dir, file), overwrite = TRUE)
}

write_sheets(
  c("j1.1_pano_vs_dmso", "j1.1_ifn_vs_dmso", "j1.1_panoifn_vs_dmso",
    "cd4_pano_vs_dmso",  "cd4_ifn_vs_dmso",  "cd4_panoifn_vs_dmso",
    "scCD4_panobinostat_vs_mock", "scCD4_vorinostat_vs_mock", "scCD4_IL2PHA_vs_mock"),
  c("J1.1_Pano_vs_DMSO", "J1.1_IFN_vs_DMSO", "J1.1_PanoIFN_vs_DMSO",
    "CD4_Pano_vs_DMSO",  "CD4_IFN_vs_DMSO",  "CD4_PanoIFN_vs_DMSO",
    "scCD4_Panobinostat_vs_Mock", "scCD4_Vorinostat_vs_Mock", "scCD4_IL2PHA_vs_Mock"),
  "DEGs_vs_DMSO_all_models.xlsx"
)

# Additional contrasts against the combined treatment
for (x in list(c("treatmentPano", "pano"), c("treatmentIFN", "ifn"))) {
  process_res(dds_j11, results(dds_j11, contrast = list(x[1], "treatmentPano_IFN")),
              paste0("j1.1_", x[2], "_vs_panoifn"))
  process_res(dds_cd4, results(dds_cd4, contrast = list(x[1], "treatmentPano_IFN")),
              paste0("cd4_",  x[2], "_vs_panoifn"))
}

write_sheets(
  c("j1.1_pano_vs_panoifn", "j1.1_ifn_vs_panoifn",
    "cd4_pano_vs_panoifn",  "cd4_ifn_vs_panoifn"),
  c("J1.1_Pano_vs_PanoIFN", "J1.1_IFN_vs_PanoIFN",
    "CD4_Pano_vs_PanoIFN",  "CD4_IFN_vs_PanoIFN"),
  "DEGs_extra.xlsx"
)
