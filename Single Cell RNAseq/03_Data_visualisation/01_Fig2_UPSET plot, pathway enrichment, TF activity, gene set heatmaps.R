library(ggplot2)
library(writexl)
library(magrittr)
library(dplyr)
library(ComplexUpset)
library(purrr)
library(ggprism)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(decoupleR)
library(pheatmap)
library(ggplotify)
#---------------------------------------------

# FIGURE 2 - B ####
## UPSET plots of DEGs


## read data
deg_list = readRDS("deg.list_hivdonors_treat.celltype.rds")
deg_list

## for each df in list, add df name as column, then convert to two cols for treat and ct

my_list = deg_list

my_fxn = function(df, i){
  
  #i = 1
  print(i)
  #name_to_add
  df_mod = df %>% 
    mutate( treatment_celltype = i) %>%
    mutate( treatment = gsub("(.*)\\.(.*)_v_Mock", replacement = "\\1", treatment_celltype),
            celltype = gsub("(.*)\\.(.*)_v_Mock", replacement = "\\2", treatment_celltype))
  return(df_mod)
  
}


result = purrr::imap(my_list, my_fxn)

## bind together in mega df for degs
deg_mega = bind_rows(result) %>%
  mutate(treatment = gsub("/", ".", treatment))

## generate UPSET plots

for (updown in c("upregulated", "downregulated")){
  for ( treat in unique(deg_mega$treatment)){
    
    ## print progress #################################
    print(paste0(updown, " - ", treat))
    
    
    if (updown == "upregulated"){
      sub = deg_mega %>%
        dplyr::filter(treatment == treat,
                      avg_log2FC < 0,
                      p_val_adj < 0.01) %>%
        mutate(celltype = factor(celltype, levels = c("TN","TCM","TTM","TREG","TEM","TPROLIF"  )))
    } else {
      sub = deg_mega %>%
        dplyr::filter(treatment == treat,
                      avg_log2FC > 0,
                      p_val_adj < 0.01) %>%
        mutate(celltype = factor(celltype, levels = c("TN","TCM","TTM","TREG","TEM","TPROLIF"  )))
    }
    
    
    ## get table for complexupset input
    tb = table(sub$gene, sub$celltype) %>%
      as.data.frame.matrix() %>%
      mutate_all(as.logical)
    
    ## generate upset plot
    upset.plot = ComplexUpset::upset(tb, names(tb), 
                                     #min_size = 10,
                                     n_intersections = 15, 
                                     keep_empty_groups = T,
                                     width_ratio = 0.2)
    
    ## save plot
    # svg
    ggsave(plot = upset.plot, filename = paste0("upset_",
                                                treat, "_",
                                                updown,
                                                ".svg"),
           width=7, height = 5)
    
  }
}


# FIGURE 2 - C ####
## Pathway enrichment analysis

names(deg_list) = gsub("/", ".", names(deg_list))

## define function to run GO analysis

run_go = function(input, name){
  
  go_res_df = data.frame()
  
  for (contr in names(input)){
    
    print(contr)
    go_input <- input[[contr]] %>%
      #go_input = testdf %>%
      
      #assign ENTREZID ID from symbol
      mutate(entrezID = mapIds(org.Hs.eg.db, keys=.$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")) %>%
      group_by(gene)%>%
      arrange(desc(avg_log2FC)) %>%
      ungroup() %>%
      dplyr::filter(p_val_adj < 0.01) %>%
      dplyr::select(c(entrezID, avg_log2FC)) %>%
      na.omit () %>% # remove unmapped
      arrange(desc(avg_log2FC)) %>%
      deframe()
    
    if (length(go_input) > 0 ){
      ## run GO enrichment
      goenr <- gseGO(geneList     = go_input,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     #minGSSize    = 10,
                     minGSSize    = 20,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     verbose      = T)
      
      
      if(is.null(goenr)){
        goenr_tidy = data.frame()
      } else {
        ## get result
        goenr_tidy <-goenr@result 
      }
      
      
      if (nrow(goenr_tidy) != 0 ){
        
        
        #if (nrow(goenr_tidy) != 0){
        ## clean and determine GeneRatio
        gene_count <- goenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        merged_res <- left_join(goenr_tidy, gene_count, by = "ID") %>% 
          mutate(GeneRatio = count/setSize) %>%
          mutate(zscore = scale(NES)) %>%
          mutate(group = contr)
        
        ## add goRes to dataframe
        go_res_df = rbind(go_res_df, merged_res)
        
      }
    } else {
      print("No DEGs")
    }
    
    
    
  }
  
  ## export to global environment
  assign(x = paste0("gse_GOBP.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("gse_GOBP.",
                                            name, ".rds"))
}

## PLOT PATHWAYS
## define function to plot pathways

## plot pathways on one dotplot
pway_comboDot <- function(contrasts_use = c(unique(pway.mega$group)),
                          num_topPways = 20,
                          minimum_NES = 1,
                          max_padj = 0.001,
                          name = "Pways_ComboDot_AllContrasts",
                          width=8,
                          height=10.5){
  
  ### PWAYS UP # select top n pathways by NES
  pways_up <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    dplyr::filter(p.adjust < max_padj, 
                  NES > minimum_NES,
                  stat == "Significant") %>%
    group_by(group) %>%
    arrange(desc(NES)) %>%
    slice_head(n = num_topPways) %>%
    pull(var= Description)
  
  pways_up <- unique(pways_up) ## get names of pathways to plot
  
  
  
  ## cluster dataframe for better graphical representation
  mat <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    subset(Description %in% pways_up) %>%
    dplyr::select(Description, NES, group) %>%
    pivot_wider(names_from = group, values_from = NES) %>%
    column_to_rownames(var = "Description")
  
  ## scale and cluster
  mat <- scale(mat)
  dist_mat <- dist(mat, method = 'euclidean')
  dist_mat[is.na(dist_mat)] = 0
  hclust_avg <- hclust(dist_mat, method = 'average')
  ## get plotting order
  label_order <- rownames(mat)[hclust_avg$order] 
  
  ## plot pways UP
  plot_up_df <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use, stat == "Significant") %>%
    mutate(group = factor(group, levels = contrasts_use)) %>%
    subset(Description %in% pways_up) %>%
    mutate(Description = factor(Description, levels =label_order ))#%>%
  #arrange(match(rownames(.), label_order))# %>%
  
  plot_up <- ggplot(plot_up_df, aes(group, Description, color=factor(stat, levels = c("NS", "Significant")))) +
    geom_point(aes(size=GeneRatio, fill= NES, stroke=1), shape=21) +
    scale_fill_gradient2(low='#0000FF', mid = "white", high='darkorange')+
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=c("black", "black"), name="Significance") +
    scale_x_discrete(drop=F)+
    #ggtitle("Pathways UP")+
    theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3),
          axis.text.y = element_text(colour = "black", face = "bold", size = 10))
  ## save as PDF and svg
  ggsave(filename = paste0("dotplot_gse.gobp_", name, "_UP.svg"),
         plot = plot_up,
         width=width, height=height)
  
  
  ## PWAYS DOWN
  pways_down <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    dplyr::filter(p.adjust < max_padj, 
                  NES < -minimum_NES,
                  stat == "Significant") %>%
    group_by(group) %>%
    arrange(NES) %>%
    slice_head(n = num_topPways) %>%
    pull(var= Description)
  
  pways_down <- unique(pways_down) ## get names of pathways to plot
  
  ## cluster dataframe for better graphical representation
  mat <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    subset(Description %in% pways_down) %>%
    dplyr::select(Description, NES, group) %>%
    pivot_wider(names_from = group, values_from = NES) %>%
    column_to_rownames(var = "Description")
  
  ## scale and cluster
  mat <- scale(mat)
  dist_mat <- dist(mat, method = 'euclidean')
  dist_mat<- as.dist(dist_mat)
  ## some NAs created
  dist_mat[is.na(dist_mat)] <- 0
  dist_mat[is.nan(dist_mat)] <- 0
  sum(is.infinite(dist_mat))  # THIS SHOULD BE 0
  hclust_avg <- hclust(dist_mat, method = 'average')
  
  ## get plotting order
  label_order <- rownames(mat)[hclust_avg$order] 
  
  ## plot pways down
  plot_down_df <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use, stat == "Significant") %>%
    mutate(group = factor(group, levels = contrasts_use)) %>%
    subset(Description %in% pways_down) %>%
    mutate(Description = factor(Description, levels =label_order ))#%>%
  #arrange(match(rownames(.), label_order))# %>%
  
  plot_down <- ggplot(plot_down_df, aes(group, Description, color=factor(stat, levels = c("NS", "Significant")))) +
    geom_point(aes(size=GeneRatio, fill= NES, stroke=1), shape=21) +
    scale_fill_gradient2(low='#0000FF', mid = "white", high='darkorange')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_color_manual(values=c("black", "black"), name="Significance") +
    ggtitle("Pathways DOWN")+
    theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3),
          axis.text.y = element_text(colour = "black", face = "bold", size = 10))
  
  ## save as PDF and svg
  ggsave(filename = paste0("dotplot_gse.gobp_", name, "_DOWN.svg"),
         plot = plot_down,
         width=width, height=height)
  
}

contrasts = unique(gse_GOBP.hivdonors_treatment.celltype$group)

## plot pathways
pway_comboDot(contrasts_use = contrasts, name = "gse.gobp_hivdonors_treat.celltype", num_topPways = 10, height =9, width =9)

# FIGURE 2 - D

# FIGURE 2 - D ####
## Transcription factor activity analysis


## set col order

col_order = c()

for (treat in c("Mock", "Vorinostat", "Panobinostat", "IL2/PHA")){
  for (ct in c("TN", "TCM", "TTM", "TEM", "TREG", "TPROLIF")){
    
    treat_ct = paste0(treat, "_", ct)
    col_order = c(col_order, treat_ct)
  }
}


# aggregated data by treatment/celltype, hiv pos only
pbulk.hiv.tc = readRDS("pbulk_cd4.hiv_tc.rds")

pbulk.hiv.tc_mtx = AggregateExpression(pbulk.hiv.tc, return.seurat = F, )

#pbulk.hiv.tc_mtx

## convert normalised RNA values to matrix
mat = pbulk.hiv.tc$RNA$data 
## remove columns with Unknown celltype
sub_names = colnames(mat)[!grepl("nonT", colnames(mat))]
# 
pbulk.hiv.tc_mtx = mat[, c(sub_names)]
## rearrange cols
pbulk.hiv.tc_mtx = pbulk.hiv.tc_mtx[, col_order]

## retrieve static file TF activity network
net = read.csv(file = "collectri-human_decouplerPy_260814.csv")
net = net %>%
  rename(mor=weight)

## run TF activity analysis

tf_acts = decoupleR::run_ulm(mat = pbulk.hiv.tc_mtx, 
                             net = net, 
                             .source = 'source', 
                             .target = 'target',
                             .mor = 'mor', 
                             minsize = 5)

## plot TF activity heatmap
n_tfs <- 30

# Transform to wide matrix
sample_acts_mat <- tf_acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- tf_acts %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

sample_acts_mat <- sample_acts_mat[tfs,]

## set hmap params
breaks2 = seq(-2, 2, by=0.1)
cols2 = colorRampPalette(rev(c("#D73027", "#FC8D59", 
                               "#FEE090", "#FFFFBF",  "#91BFDB", "#4575B4", "#235F97")))(n=length(breaks2)-1 )



sample_acts_mat = sample_acts_mat[, col_order]
head(sample_acts_mat)
sample_acts_mat2 = as.matrix(sample_acts_mat)
# Plot
hmap_tf = pheatmap::pheatmap(mat = sample_acts_mat2,
                             color = cols2,
                             #border_color = "white",
                             breaks = breaks2,
                             cellwidth = 15,
                             cellheight = 15,
                             treeheight_row = 20,
                             treeheight_col = 20,
                             scale = "row",
                             cluster_cols  = F)

hmap_tf.gg = as.ggplot(hmap_tf)

ggsave(plot = hmap_tf.gg, filename = paste0("hmap_",
                                            "hiv_donors_TFactivity_top30mostVariableTFs_",
                                            ".svg"),
       width = 8, height =8)




# FIGURE 2 - E ####
## Heatmap for T cell signalling and innate immunity genes

## read in gene set list

genesets = readxl::read_xlsx("GeneSets.xlsx")

## set column order
col_order = c()


for (treat in c("Mock", "Vorinostat", "Panobinostat", "IL2/PHA")){
  for (ct in c("TN", "TCM", "TTM", "TEM", "TREG", "TPROLIF")){
    for (donor in c("d1", "d2", "d3")){
      
      donor_treat_ct = paste0(donor, "_", treat, "_", ct)
      col_order = c(col_order, donor_treat_ct)
    }
  }
}


## read expression data by donor/treatment/celltype
pbulk.hiv.dtc = readRDS("pbulk_cd4.hiv_dtc.rds")

## convert normalised RNA values to matrix
mat = pbulk.hiv.dtc$RNA$data 
## remove columns with Unknown celltype
sub_names = colnames(mat)[!grepl("nonT", colnames(mat))]

pbulk.hiv.dtc_mtx = mat[, c(sub_names)]
## rearrange cols
pbulk.hiv.dtc_mtx = pbulk.hiv.dtc_mtx[, col_order]

## generate heatmap -- innate immune genes

gs1 = genesets %>% 
  dplyr::filter(set == "innate_immunity") %>%
  pull(gene)
mat_gs1 = pbulk.hiv.dtc_mtx[gs1,] %>%
  as.matrix()
mat_gs1

hmap_dtc_ii =pheatmap(mat_gs1,scale = "row", cluster_rows = T, cluster_cols = F,
                      breaks = breaks2,
                      color = cols2,
                      gaps_col = c(18, 36, 54))

hmap_dtc_ii.gg = as.ggplot(hmap_dtc_ii)

ggsave(plot = hmap_dtc_ii.gg, filename = paste0("hmap_",
                                                "hiv_donors_genesets_byDonorCelltypeTreatment_",
                                                "innateimmunity",
                                                ".svg"),
       width = 12, height =6)

## generate heatmap -- t cell signalling

gs2 = genesets %>% 
  dplyr::filter(set == "tcell_signalling") %>%
  pull(gene)
mat_gs2 = pbulk.hiv.dtc_mtx[gs2,] %>%
  as.matrix()
mat_gs2

hmap_dtc_tcs =pheatmap(mat_gs2,scale = "row", cluster_rows = T, cluster_cols = F,
                       breaks = breaks2,
                       color = cols2,
                       gaps_col = c(18, 36, 54))

hmap_dtc_tcs.gg = as.ggplot(hmap_dtc_tcs)

ggsave(plot = hmap_dtc_tcs.gg, filename = paste0("hmap_",
                                                 "hiv_donors_genesets_byDonorCelltypeTreatment_",
                                                 "tcellsignalling",
                                                 ".svg"),
       width = 12, height =6)




