library(Seurat)
library(ggplot2)
library(scCustomize)

##########################################


#read in data
cd4=readRDS("cd4_all.int.rpca2.rds")

# set identity of cluster resolution
Idents(cd4) = "rpca.clusters_0.8"


# create cluster dotplot of celltype markers
Clustered_DotPlot(seurat_object = cd4, 
                  # TN
                  features = c("CD3D", "CD4", "CD8A", "CCR7", "S100A4", 
                               # TCM vs TTM
                               "CD28", "NKG7", "GZMB", "GZMH", "PRF1", "GNLY",
                               # IL2PHA specific TEMs
                               "IFNG", "CCL3", "CCL4",
                               # TREG
                               "IKZF2", "FOXP3", "CTLA4", "IL2RA",
                               # T PROLIF
                               "TOP2A"),
                  exp_color_min=-1, flip = T,
                  cluster_feature=F,
                  cluster_ident=F,
                  show_ident_colors =F,
                  colors_use_exp = viridis_plasma_light_high)+
  # scale_color_viridis_c(direction = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


