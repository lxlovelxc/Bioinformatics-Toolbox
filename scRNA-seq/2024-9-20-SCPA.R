# library(devtools)
# # install.packages("devtools")
# devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
# devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
# devtools::install_github("jackbibby1/SCPA")

library(SCPA)
library(msigdbr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
data <- readRDS("2024-9-13-with-additional-embryo-raw-data.Rdata")
data1 <- subset(data, subset = nFeature_RNA > 600 )
data1
data1 <- JoinLayers(data1)
data1 <- NormalizeData(data1)
data1 <- FindVariableFeatures(data1, selection.method = "mean.var.plot")
data1 <- ScaleData(data1,model.use = 'linear',do.scale = T,do.center = T)
data1 <- RunPCA(data1, features = NULL)
data1 <- FindNeighbors(data1, dims = 1:15)#20 /25 is the best choice
data1 <- FindClusters(data1, resolution = 0.1)
data1 <- RunUMAP(data1,reduction = "pca",dims = 1:15)
# "C2","CP"

DimPlot(data1) +
  theme(aspect.ratio = 1)

young <- seurat_extract(data1,
                          meta1 = "label", value_meta1 = "young",
                        meta2 = "seurat_clusters",value_meta2 = "1")

old <- seurat_extract(data1,
                            meta1 = "label", value_meta1 = "old",
                            meta2 = "seurat_clusters",value_meta2 = "1")


pathways <- msigdbr("Homo sapiens", "C2","CP") %>%
  format_pathways()



scpa_out <- compare_pathways(samples = list(young, old), 
                             pathways = pathways,
                             parallel = TRUE,
                             cores = 10)
head(scpa_out,10)
plot_rank(scpa_out = scpa_out, 
          pathway = c("ECM","WNT"), 
          base_point_size = 2, 
          highlight_point_size = 5)


plot_heatmap(scpa_out,
             column_names = "Young vs Old",
             show_row_names = T)

#msigdbr_collections()