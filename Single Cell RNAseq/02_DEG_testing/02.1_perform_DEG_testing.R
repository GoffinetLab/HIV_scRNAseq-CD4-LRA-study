library(Seurat)
library(ggplot2)
library(writexl)
library(magrittr)
library(dplyr)
library(DESeq2)

#----------------------------------------------------------------------


## read in pseudobulk results for results by DONOR, TREATMENT, CELLTYPE
cd4 = readRDS("pbulk_cd4.hiv_dtc.rds")

## add treatment celltype meta.data column
cd4$treatment_celltype = paste0(cd4$treatment, "_", cd4$celltype)
cd4@meta.data


## perform pseudobulk DEG testing using DESeq2 algorithm, by treatment/celltype combo
Idents(cd4) = "treatment_celltype"

deg.list_dt = list()

for (treat in c("Panobinostat", "Vorinostat", "IL2/PHA")){
  for (ct in c("TN", "TCM", "TTM", "TREG", "TEM", "TPROLIF")){
    
    
    contrast.name = paste0(treat, ".", ct, "_v_Mock")
    
    print(contrast.name)
    
    id1 = paste0(treat, "_", ct)
    id2 = paste0("Mock", "_", ct)

    
    
    deg <- FindMarkers(object = cd4, 
                       ident.1 = id1, 
                       ident.2 = id2,
                       test.use = "DESeq2")
    
    ## save gene name as column
    deg = deg %>%
      mutate(gene = rownames(.))%>%
      relocate(gene)
    
    
    deg.list_dt[[contrast.name]] = deg
    
  }
}

##rename
names(deg.list_dt) = contrast.name = gsub(pattern = "/", replacement = ".", names(deg.list_dt))

## save results as .rds object for later retrieval

saveRDS(object = deg.list_dt, file = "deg.list_hivdonors_treat.celltype.rds")