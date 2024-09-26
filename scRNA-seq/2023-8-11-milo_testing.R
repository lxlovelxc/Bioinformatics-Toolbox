# Load dependencies

library(SingleCellExperiment)
library(miloDE)
library(Rtsne)
library(scuttle)
suppressMessages(library(miloR))
suppressMessages(library(uwot))
library(scran)
suppressMessages(library(dplyr))
library(reshape2)
library(scWGCNA)
suppressMessages(library(Seurat))
library(ggplot2)
library(viridis)
library(ggpubr)
library(BiocParallel)

ncores = 6
mcparam = MulticoreParam(workers = ncores)
register(mcparam)

seurat_object <- readRDS("cells_0922.Rdata")
sce_object <- SingleCellExperiment(
  assays = list(counts = seurat_object@assays$RNA@counts),
  colData = seurat_object@meta.data
)
sce_object$label <- seurat_object$label
saveRDS(sce_object, "cells_0922_sce.rds")

all <- readRDS("cells_0922_sce.rds")
dim(counts(all))
all <- logNormCounts(all)
all <- reducedDim(all,PCA)
set.seed(32)
umaps = as.data.frame(uwot::umap(reducedDim(all , "pca.corrected")))

## UMAPs, colored by cluster
# umaps = cbind( as.data.frame(reducedDim(all , "UMAP")) , as.data.frame(colData(all)))
# umaps$celltype[umaps$seurat_clusters == "0"] = "Cytotrophoblast"
# umaps$celltype[umaps$seurat_clusters == "1"] = "Syncytiotrophoblast"
# umaps$celltype[umaps$seurat_clusters == "2"] = "Cluster 2"
# umaps$celltype[umaps$seurat_clusters == "3"] = "Epiblast"
# umaps$celltype[umaps$seurat_clusters == "4"] = "Hypoblast"
# 
# p = ggplot(umaps , aes(x = UMAP_1 , y = UMAP_2 , col = celltype)) + 
#   geom_point() + 
#   theme_bw() +
#   theme(legend.position = "left") + 
#   facet_wrap(~label) +
#   labs(x = "" , y = "") +
#   theme(legend.position = "right" , legend.text=element_text(size=12)) +
#   theme(strip.text = element_text(size = 12) , legend.title=element_text(size=12))
# p


#Assign nhoods
set.seed(15)
stat_k = estimate_neighbourhood_sizes(all, k_grid = seq(10,40,5) ,
                                      order = 2, prop = 0.1 , filtering = TRUE,
                                      reducedDim_name = "PCA" , plot_stat = TRUE)


sce_humanEmbryo = assign_neighbourhoods(all, k = 30, order = 2, 
                                        filtering = TRUE, reducedDim_name = "PCA", verbose = F)

nhoods_sce = nhoods(sce_humanEmbryo)

#######################################################

#
sce_humanEmbryo$celltype[sce_humanEmbryo@colData$seurat_clusters == "0"] = "Cytotrophoblast"
sce_humanEmbryo$celltype[sce_humanEmbryo@colData$seurat_clusters == "1"] = "Syncytiotrophoblast"
sce_humanEmbryo$celltype[sce_humanEmbryo@colData$seurat_clusters == "2"] = "Cluster 2"
sce_humanEmbryo$celltype[sce_humanEmbryo@colData$seurat_clusters == "3"] = "Epiblast"
sce_humanEmbryo$celltype[sce_humanEmbryo@colData$seurat_clusters == "4"] = "Hypoblast"
# 
# 
# add_da = function(sce_milo , reducedDim_name = "PCA" , col_annotate){
#   require(dplyr)
#   require(miloR)
#   sce_milo <- countCells(sce_milo, meta.data = as.data.frame(colData(sce_milo)), samples="embryo")
#   sce_design <- data.frame(colData(sce_milo))[,c("embryo", "label")]
#   sce_design <- distinct(sce_design)
#   rownames(sce_design) <- sce_design$sample
#   nhoodCounts(sce_milo) = as.matrix(nhoodCounts(sce_milo) , "dgCMatrix")
#   da_results <- testNhoods(sce_milo, design = ~ factor(label), design.df = sce_design, fdr.weighting = "graph-overlap", reduced.dim = reducedDim_name)
#   if (!is.null(col_annotate)){
#     da_results = annotateNhoods(sce_milo, da_results, coldata_col = col_annotate)
#   }
#   da_results = da_results[order(da_results$Nhood) , ]
#   return(da_results)
# }
# 
# 
# da_stat = add_da(sce_humanEmbryo , col_annotate = "celltype")
# 
# p2 = plot_milo_by_single_metric(sce_humanEmbryo , nhood_stat = da_stat , colour_by = "logFC" , significance_by = "SpatialFDR", alpha = 0.1, size_range = c(2,5)) + 
#   guides(size=FALSE, width=FALSE ,edge_width=FALSE) +
#   theme(legend.position = "top")


# assign cell types for nhoods
nhood_stat_ct = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_center = colnames(nhoods_sce))
nhood_stat_ct = miloR::annotateNhoods(sce_humanEmbryo , nhood_stat_ct , coldata_col =
                                        "celltype")
nhood_stat_ct = miloR::annotateNhoods(sce_humanEmbryo , nhood_stat_ct , coldata_col =
                                        "label")

p = plot_milo_by_single_metric(sce_humanEmbryo, nhood_stat_ct, colour_by = "label" ,
                               layout = "UMAP" , size_range = c(1.5,3) , edge_width =
                                 c(0.2,0.5))

p


p1 = plot_milo_by_single_metric(sce_humanEmbryo , nhood_stat = nhood_stat_ct , colour_by = "celltype" ) +   
  theme(legend.position = "top")

p1


# #Calculate AUC per neighbourhood
# stat_auc = suppressWarnings(calc_AUC_per_neighbourhood(sce_humanEmbryo , sample_id = "embryo" ,
#                                                        condition_id = "label" ,min_n_cells_per_sample = 1, BPPARAM = mcparam))
# p = plot_milo_by_single_metric(sce_humanEmbryo, stat_auc, colour_by = "auc" ,
#                                layout = "UMAP" , size_range = c(1.5,3) , edge_width =
#                                  c(0.2,0.5)) +
#   scale_fill_viridis(name = "AUC")
# p


#miloDE 
#subset_nhoods = stat_auc$Nhood[!is.na(stat_auc$auc)],

de_stat = de_test_neighbourhoods(sce_humanEmbryo ,
                                 sample_id = "embryo",
                                 design = ~label,
                                 covariates = c("label","batch"),
                                 contrasts = c("labelyoung"),
                                 output_type = "SCE",
                                 plot_summary_stat = TRUE,
                                 layout = "UMAP", BPPARAM = mcparam , 
                                 verbose = T)

stat_de_magnitude = rank_neighbourhoods_by_DE_magnitude(de_stat)

p1 = plot_milo_by_single_metric(sce_humanEmbryo, stat_de_magnitude, colour_by = "n_DE_genes" , 
                                layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# DE genes")
p2 = plot_milo_by_single_metric(sce_humanEmbryo, stat_de_magnitude, colour_by = "n_specific_DE_genes" , 
                                layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# specific\nDE genes" , option = "inferno")
p1
p2

# de_stat@assays@data@listData$logFC[is.na(de_stat@assays@data@listData$logFC)] = 0
# de_stat@assays@data@listData$pval[is.na(de_stat@assays@data@listData$pval)] = 1
# de_stat@assays@data@listData$pval_corrected_across_genes[is.na(de_stat@assays@data@listData$pval_corrected_across_genes)] = 1
# de_stat@assays@data@listData$pval_corrected_across_nhoods[is.na(de_stat@assays@data@listData$pval_corrected_across_nhoods)] = 1
# 
# n_hoods_sig = sum(de_stat@assays@data@listData$pval_corrected_across_nhoods < 0.1 , na.rm = T)
###################################################################################
#WGCNA
get_wgcna_modules = function(de_stat , subset_hoods = NULL , 
                             n_hoods_sig.thresh = 2 ,
                             npcs = 5 ,
                             pval.thresh = 0.1 ){
  require(scWGCNA)
  require(Seurat)
  require(dplyr)
  require(reshape2)
  
  set.seed(32)
  # subset hoods
  if (!is.null(subset_hoods)){
    de_stat = de_stat[de_stat$Nhood %in% subset_hoods , ]
  }
  
  # focus on genes that DE in at least 2 nhoods
  de_stat_per_gene = as.data.frame(de_stat %>% group_by(gene) %>% dplyr::summarise(n_hoods_sig = sum(pval_corrected_across_nhoods < pval.thresh , na.rm = TRUE)))
  genes_sig = de_stat_per_gene$gene[de_stat_per_gene$n_hoods_sig >= n_hoods_sig.thresh]
  
  de_stat = de_stat[de_stat$gene %in% genes_sig, ]
  de_stat = de_stat[order(de_stat$Nhood) , ]
  
  # discard neighbourhoods in which testing was not performed
  de_stat = de_stat[de_stat$test_performed , ]
  
  # for this analysis, set logFC to 0 and pvals to 1 if they are NaN
  de_stat$logFC[is.na(de_stat$logFC)] = 0
  de_stat$pval[is.na(de_stat$pval)] = 1
  de_stat$pval_corrected_across_genes[is.na(de_stat$pval_corrected_across_genes)] = 1
  de_stat$pval_corrected_across_nhoods[is.na(de_stat$pval_corrected_across_nhoods)] = 1
  
  # set logFC to 0 if pval_corrected_across_nhoods > pval.thresh
  de_stat$logFC[de_stat$pval_corrected_across_nhoods >= pval.thresh] = 0
  
  # move the object to Seurat
  de_stat = reshape2::dcast(data = de_stat, formula = gene~Nhood, value.var = "logFC")
  rownames(de_stat) = de_stat$gene
  de_stat = de_stat[,2:ncol(de_stat)]
  
  obj.seurat <- CreateSeuratObject(counts = de_stat)
  DefaultAssay(obj.seurat) <- "RNA"
  obj.seurat = FindVariableFeatures(obj.seurat)
  # scale
  obj.seurat[["RNA"]]@scale.data = as.matrix(obj.seurat[["RNA"]]@data)
  obj.seurat = RunPCA(obj.seurat , npcs = npcs)
  
  # run scwgcna
  clusters_scwgcna = run.scWGCNA(p.cells = obj.seurat, 
                                 s.cells = obj.seurat, 
                                 is.pseudocell = F, 
                                 features = rownames(obj.seurat),
                                 less = TRUE , merging = TRUE)
  # compile stat
  clusters = lapply(1:length(clusters_scwgcna$module.genes) , function(i){
    out = data.frame(cluster = i , gene = clusters_scwgcna$module.genes[[i]] , n_genes = length(clusters_scwgcna$module.genes[[i]]))
    return(out)
  })
  clusters = do.call(rbind , clusters)
  # add colors
  genes_w_colors = clusters_scwgcna$dynamicCols
  genes_w_colors = data.frame(gene = names(genes_w_colors) , cluster_color = genes_w_colors)
  clusters = merge(clusters , genes_w_colors)
  
  return(clusters)
}


# for this vignette, for simplicity we will focus on genes that are DE in at least 4 neighbourhoods
modules_wgcna = suppressMessages(get_wgcna_modules(convert_de_stat(de_stat) , n_hoods_sig.thresh = 2))

# # rearrange modules in the descending order
# tab = sort(table(modules_wgcna$cluster) , decreasing = T)
# tab = data.frame(cluster = c(1:length(tab)) , cluster_old = names(tab))
# colnames(modules_wgcna)[colnames(modules_wgcna) == "cluster"] = "cluster_old"
# modules_wgcna = merge(modules_wgcna , tab , all.x = T)
# n_hoods_sig = sum(pval_corrected_across_nhoods < 0.1 , na.rm = T)
# 
# 
# n_hoods_total = ncol(nhoods(sce_humanEmbryo))
# tab = as.data.frame(table(modules_wgcna$cluster))
# colnames(tab) = c("module" , "n_genes")
# cols_modules = MetBrewer::met.brewer("Homer2" , n =13 )
# p3 = ggplot(tab , aes(x = factor(module) , y = n_genes , fill = factor(module))) +
#   geom_bar(stat = "identity" , position = "dodge") +
#   theme_bw() +
#   scale_fill_manual(values = cols_modules) +
#   theme(legend.position = "none") +
#   labs(x = "Gene module" , y = "Number of genes") +
#   theme(text = element_text(size = 14)) 
# p2 = ggplot(modules_wgcna , aes(x = factor(cluster) , y = n_hoods_sig/n_hoods_total, col = factor(cluster) )) +
#   geom_boxplot() + geom_jitter(width = 0.1 , size = .75) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "Gene module" , y = "Fraction of \n significant hoods") +
#   scale_color_manual(values = cols_modules) +
#   theme(text = element_text(size = 14)) 

# #breakdown by CT
# modules_wgcna$cluster_name = paste0("Module #" , modules_wgcna$cluster)
# modules_wgcna$cluster_name = factor(modules_wgcna$cluster_name , levels = paste0("Module #",levels(modules_wgcna$cluster)))
# 
# 
# plots = lapply(levels(modules_wgcna$cluster_name) , function(cluster_name){
#   genes = modules_wgcna$gene[modules_wgcna$cluster_name == cluster_name]
#   p = plot_beeswarm_gene_set(de_stat , genes = genes , nhoodGroup = "celltype" , size = 3.5) + 
#     scale_color_gradient2(breaks = c(-2,-1,0,1,2) , limits = c(-2.2,2.2) , name = "avg logFC") + 
#     geom_quasirandom(alpha = 1 , size = 3.5 , shape = 1 , colour = "black") +
#     ylim(c(0,1)) + 
#     labs(x = "",y = "") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#     theme(axis.text.y = element_text(face="bold", colour = cols_ct, size=12)) + 
#     ggtitle(cluster_name)
#   return(p)
# })


## Gene module plots

plots = lapply(sort(unique(modules_wgcna$cluster)) , function(cluster){
  p = plot_DE_gene_set(sce_humanEmbryo, de_stat , genes = modules_wgcna$gene[modules_wgcna$cluster == cluster],
                       layout = "UMAP" , size_range = c(0.5,3) ,
                       node_stroke = 0.3, edge_width = c(0.2,0.5)) +
    ggtitle(paste0("Module ", cluster, ", ", length(modules_wgcna$gene[modules_wgcna$cluster == cluster]) , " genes"))
  return(p)
})
p = ggarrange(plotlist = plots)
p

#p = ggarrange(plotlist = plots)
#p

## Module 2
set.seed(1020)
module = 2
n_genes = 6
genes = sample(modules_wgcna$gene[modules_wgcna$cluster == module] , n_genes)


plots = lapply(genes , function(gene){
  p = plot_DE_single_gene(sce_humanEmbryo, de_stat , gene = gene , layout = "UMAP" , set_na_to_0 = TRUE) + 
    ggtitle(gene)
  return(p)
})
p = ggarrange(plotlist = plots , ncol = 2, nrow = 3)
p


##
set.seed(1020)
module = 3
n_genes = 2
genes = sample(modules_wgcna$gene[modules_wgcna$cluster == module] , n_genes)


plots = lapply(genes , function(gene){
  p = plot_DE_single_gene(sce_humanEmbryo, de_stat , gene = gene , layout = "UMAP" , set_na_to_0 = TRUE) + 
    ggtitle(gene)
  return(p)
})
p = ggarrange(plotlist = plots, ncol=2, nrow=3)
p


cluster = 2
gos = get_go_terms(modules_wgcna$gene[modules_wgcna$cluster == cluster] )

genes_cluster_1_selected = c("S100A6" , "SLC7A8", "MDN1" , "ZNF337" , "SLC6A6" , "S100A13")

plot_umap_w_counts = function(gene , sce_milo , umap_name = "UMAP" , size = 0.75){
  umaps = cbind( as.data.frame(reducedDim(sce_milo , umap_name)) , as.data.frame(colData(sce_milo)))
  umaps$counts = as.numeric(logcounts(sce_milo[gene , ]))
  umaps = umaps[order(umaps$counts) , ]
  p = ggplot(umaps , aes(x = UMAP_1 , y = UMAP_2 , col = counts)) +
    geom_point(size=size) +
    scale_color_viridis(discrete = F) +
    theme_bw() +
    facet_wrap(~label) +
    xlim(c(min(umaps$umap_1) , max(umaps$umap_1))) + ylim(c(min(umaps$umap_2) , max(umaps$umap_2))) +
    xlim(c(-2,12)) + ylim(c(-10,0)) +
    ggtitle(gene)
  return(p)
}

plots = lapply(genes_cluster_1_selected, function(gene){
  print(gene)
  p = plot_umap_w_counts(gene , sce_humanEmbryo , size = .75) + 
    labs(x = "" , y = "") +
    theme(legend.position = "top") 
  return(p)
})
p





get_go_terms = function(genes , pval.thresh = 0.1){
  require(enrichR)
  require(stringr)
  gos <- as.data.frame( enrichr(genes, 'GO_Biological_Process_2021')[[1]] )
  gos$n_tested_genes = length(genes)
  gos$n_detected_genes_in_go = NaN
  gos$n_total_genes_in_go = NaN  
  
  gos$Term_id = NaN
  for (i in 1:nrow(gos)){
    current.overlap = gos$Overlap[i]
    current.overlap = str_split(current.overlap, "/")
    gos$n_detected_genes_in_go[i] = as.numeric(current.overlap[[1]][1])
    gos$n_total_genes_in_go[i] = as.numeric(current.overlap[[1]][2])
    
    current.go = gos$Term[i]
    current.go = str_split(current.go, "GO:")
    gos$Term[i] = substring(current.go[[1]][1], 1, nchar(current.go[[1]][1])-2)
    gos$Term_id[i] = substring(current.go[[1]][2], 1, nchar(current.go[[1]][2])-1)
  }
  gos = gos[gos$Adjusted.P.value < pval.thresh, ]
  return(gos)
}