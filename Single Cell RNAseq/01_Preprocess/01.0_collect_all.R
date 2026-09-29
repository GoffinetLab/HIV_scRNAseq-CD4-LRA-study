library(Seurat)
library(ggplot2)
library(readxl)
library(dplyr)
library(harmony)

## ------------------------------------------------

## read in data from HIV+ donors and healthy control ####

cd4_hiv = Read10X_h5(file  = "/path/to/aggHIVpos/outs/count/filtered_feature_bc_matrix.h5")

cd4_hiv <- CreateSeuratObject(counts = cd4_hiv, project = "cd4.hiv", min.cells = 3, min.features = 200)

# % MT reads
cd4_hiv[["percent.mt"]] <- PercentageFeatureSet(cd4_hiv, pattern = "^MT-")
# Visualize QC metrics as a violin plot
vln = VlnPlot(cd4_hiv, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("qc/vln_qc.png", plot = vln)

# subset according to QC results
cd4_hiv = subset(cd4_hiv, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)


#### read in other independent datasets ####


## read Wang 2020 PBMC dataset - CD4 cells only and 10X PBMC - CD4 only

cd4_wang = readRDS("/path/to/pbmcwang2020.rds")
cd4_10x = readRDS("/path/to/pbmc20k_cd4.rds")


## subset processed data to counts only

# 10X dataset

DefaultAssay(cd4_10x) = "RNA"
cd4_10x_diet = DietSeurat(cd4_10x, layers =NULL, assays=c("RNA"), dimreducs = NULL, graphs = NULL, misc = FALSE)

## remove meta data columns not needed
cd4_10x_diet[["nCount_SCT"]] = NULL
cd4_10x_diet[["nFeature_SCT"]] = NULL
cd4_10x_diet[["SCT_snn_res.0.8"]] = NULL
cd4_10x_diet[["seurat_clusters"]] = NULL

# Wang 2020

DefaultAssay(cd4_wang) = "RNA"
cd4_wang_diet = DietSeurat(cd4_wang, layers =NULL, assays=c("RNA"), dimreducs = NULL, graphs = NULL, misc = FALSE)
## remove meta data columns not needed
cd4_wang_diet[["nCount_SCT"]] = NULL
cd4_wang_diet[["nFeature_SCT"]] = NULL
cd4_wang_diet[["SCT_snn_res.0.8"]] = NULL
cd4_wang_diet[["seurat_clusters"]] = NULL


#### MERGE ALL OBJECTS ####

cd4_all = merge(x = cd4_hiv, y = list(cd4_10x_diet, cd4_wang_diet))
cd4_all[["RNA"]] <- JoinLayers(cd4_all[["RNA"]])


## assign libraries 
cd4_all$bc = rownames(cd4_all@meta.data) 
cd4_all$lib = gsub(".*-(.*)", replacement = "\\1", x=cd4_all$bc)

## read in 
# library assign 
libassign = read.csv(file = "lib_assign.csv", sep=";")
# hiv barcode assign --- read in list of HIV+ cells
hivassign = readxl::read_xlsx(path =  "hivassign.xlsx")

## pull out meta.data and assign lib information and hiv information

meta = cd4_all@meta.data

meta = meta %>%
  plyr::join(., libassign) %>%
  plyr::join(., hivassign) %>%
  mutate(hiv_counts = ifelse(is.na(hiv_counts), 0, hiv_counts)) %>%
rownames(meta) = meta$bc


## join info to meta data
cd4_all@meta.data = plyr::join(cd4_all@meta.data, meta)
## correct rownames to cell barcodes 
rownames(cd4_all@meta.data) = cd4_all$bc


## save merged object

saveRDS(object = cd4_all, file="cd4_all.rds")




