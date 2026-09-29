library(Seurat)
library(ggplot2)
library(readxl)
#----------------------------------------------------------------------------


## read data
cd4 = readRDS("cd4_all.int.rpca2.ct.rds")

# Generate UMAP projections 

# FIGURE 1 - A

## UMAP coloured by donor, split by HIV positivity status

param = "donor"

Idents(cd4) = param
umap = DimPlot(cd4, label = FALSE, reduction = "umap.rpca",
               raster=F, split.by = "disease",
               shuffle=T)
ggsave(paste0("umap_alldonors_", param, "_split.HIVposneg.svg"), plot = umap, width=18, height=7)

# FIGURE 1 - B, Figure S1 - A

## UMAP coloured by:
## - celltype
## - cluster

params = c("rpca.clusters_0.8", "celltype")


for (param in params){
  Idents(cd4) = param
  umap = DimPlot(cd4, label = FALSE, reduction = "umap.rpca",
                 raster=F,
                 shuffle=T)
  ggsave(paste0("umap_alldonors_", param, ".svg"), plot = umap, width=9, height=7)
}


# FIGURE 1 - D

## UMAP of samples from HIV+ donors only, coloured by treatment:

### read data - HIV+ only samples
cd4.hiv = readRDS("cd4_hiv.int.rpca2.ct.rds")

params = c("treatment", "donor")


for (param in params){
  Idents(cd4.hiv) = param
  umap = DimPlot(cd4.hiv, label = FALSE, reduction = "umap.rpca",
                 raster=F,
                 shuffle=T)
  ggsave(paste0("umap_hivdonors_", param, ".svg"), plot = umap, width=9, height=7)
}

# FIGURE 1 - G

## extract UMAP emeddings, with HIV+/- detection for each cell

cd4.embed = FetchData(cd4.hiv, vars=c("donor", "treatment", "celltype", "hiv_counts", "umaprpca_1", "umaprpca_2" ))

## UMAP of HIV+ cells, coloured by treatment

cd4.embed_hiv = cd4.embed %>%
  dplyr::filter(hiv_counts > 0)

cd4.embed_neg = cd4.embed %>%
  dplyr::filter(hiv_counts == 0)

umap_hivcells = ggplot( )+
  geom_point(data=cd4.embed_neg, aes(umaprpca_1, umaprpca_2), colour = "grey", size=0.5)+
  geom_point(data=cd4.embed_hiv, aes(umaprpca_1, umaprpca_2, colour=treatment), size=1.5)+
  theme_prism()


ggsave(plot = umap_hivcells, filename = "umap_HIVposCells_by.treatment.svg", 
       width = 9, height=7)
