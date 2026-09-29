library(Seurat)
library(ggplot2)
library(readxl)
#----------------------------------------------------------------------------

## read data
cd4 = readRDS("cd4_all.int.rpca2.ct.rds")

## create treat_celltype in meta.data
cd4@meta.data$treat_celltype = paste0(cd4@meta.data$treatment, "_", cd4@meta.data$celltype)

## subset to hiv+ donors only
Idents(cd4) = "donor"
cd4.hiv = subset(cd4, idents= c("d1", "d2", "d3"))

## create aggr expression objects by

## donor_treat_celltype
cd4.hiv.pb_dtc <- AggregateExpression(cd4.hiv, assays = "RNA", return.seurat = T, group.by = c("donor", "treatment", "celltype"))
## save
saveRDS(object = cd4.hiv.pb_dtc, file= "pbulk_cd4.hiv_dtc.rds")
#saveRDS(object = cd4.hiv.pb_dtc, file= "pbulk_toy_dtc.rds")

## donor_treat
cd4.hiv.pb_dt <- AggregateExpression(cd4.hiv, assays = "RNA", return.seurat = T, group.by = c("donor", "treatment"))
## save
saveRDS(object = cd4.hiv.pb_dt, file= "pbulk_cd4.hiv_dt.rds")
#saveRDS(object = cd4.hiv.pb_dt, file= "pbulk_toy_dt.rds")

## treat_celltype
cd4.hiv.pb_tc <- AggregateExpression(cd4.hiv, assays = "RNA", return.seurat = T, group.by = c("treatment", "celltype"))
#cd4.hiv.pb_tc <- AggregateExpression(cd4, assays = "RNA", return.seurat = T, group.by = c("treatment", "celltype"))
saveRDS(object = cd4.hiv.pb_tc, file= "pbulk_cd4.hiv_tc.rds")


#saveRDS(object=cd4.hiv, file = "toy_hiv.int.ct.rds")
saveRDS(object=cd4.hiv, file = "cd4_hiv.int.rpca2.ct.rds")
