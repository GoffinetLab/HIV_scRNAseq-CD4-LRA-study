library(Seurat)
library(ggplot2)
library(readxl)
library(dplyr)
library(harmony)

## ------------------------------------------------
options(future.globals.maxSize = 3e+09)
set.seed(42069)

print("_____________ RPCA INTEGRATION _____________")

## load data will all cells

obj = readRDS("cd4_all.rds")

##
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$lib)


## perform standard pre-proc for UNINTEGRATED
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)


## perform SCTransformation
# run sctransform
obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
## Dim reduction

obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 1, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#### PERFORM INTEGRATION

## RPCA

obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = T,
  new.reduction = "integrated.rpca"#,
  #k.weight=50
)


## cluster at multiple resolutions
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "rpca.clusters_0.8")
obj <- FindClusters(obj, resolution = 1, cluster.name = "rpca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")



saveRDS(object = obj, file = "cd4_all.int.rpca.rds")



###############################################

