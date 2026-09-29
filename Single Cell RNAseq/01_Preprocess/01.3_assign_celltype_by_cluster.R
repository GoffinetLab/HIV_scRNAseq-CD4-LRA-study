library(Seurat)
library(ggplot2)
library(readxl)
######################################################

## read data

#cd4 = readRDS("toy_int.rds")
cd4 = readRDS("cd4_all.int.rpca2.rds")


## assigned ct to clusters

ct_assign = readxl::read_excel("ct_assign.xlsx")

## pull meta.data
meta = cd4@meta.data

## match celltype to rpca clusters 0.8
meta = plyr::join(meta, ct_assign)

## assign column w/ ct information to meta.data of original object
cd4@meta.data$celltype = meta$celltype

## join layers of RNA assay together to perform DEG analysis

cd4[["RNA"]] = JoinLayers(cd4[["RNA"]])

## save RDS
saveRDS(object = cd4, file = "cd4_all.int.rpca2.ct.rds")
