library(EnhancedVolcano)
library(Seurat)
cell1 <- readRDS("cells_0922.Rdata")
Idents(cell1) <- "label"
de.markers <- FindMarkers(cell1,ident.1 = "old",ident.2 = "young",only.pos = TRUE,
                          test.use = "DESeq2")#default logfc = 0.25
head(de.markers,n = 10)
EnhancedVolcano(de.markers,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = rownames(de.markers),
                FCcutoff = 0.5,
                pCutoff = 0.05,
                title = "genes differential expression")


cluster0 <- subset(cell1,subset = seurat_clusters == "0")
cluster0.de.markers <- FindMarkers(cluster0,ident.1 = "old",ident.2 = "young",only.pos = FALSE,
                                   test.use = "DESeq2",logfc.threshold = 0.25)
cluster0.de.markers <- as.data.frame(cluster0.de.markers)
#p val cutoff = 0.05
cluster0.selected.markers <- cluster0.de.markers[cluster0.de.markers$p_val_adj <= 0.05 &
                                                   abs(cluster0.de.markers$avg_log2FC) >= 0.25,]
#cluster0.selected.markers. <- cluster0.de.markers[cluster0.de.markers$avg_log2FC <= -0.5 & 
#                                                   cluster0.de.markers$p_val_adj <= 0.05,]
cluster0.selected.markers <- na.omit(cluster0.selected.markers)
head(cluster0.selected.markers,n = 10)
write.csv(cluster0.selected.markers,"c0-0.25.csv")
EnhancedVolcano(cluster0.selected.markers,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = rownames(cluster0.selected.markers),
                FCcutoff = 0.25,
                pCutoff = 0.05,
                title = "cluster0 differential expression(log2FC= 0.25)")
VlnPlot(cluster0,"NNMT",group.by = "embryo",split.by = "label")


cluster1 <- subset(cell1,subset = seurat_clusters == "1")
cluster1.de.markers <- FindMarkers(cluster1,ident.1 = "old",ident.2 = "young",only.pos = FALSE,
                                   test.use = "DESeq2",logfc.threshold = 0.25)

cluster1.de.markers <- as.data.frame(cluster1.de.markers)
cluster1.selected.markers <- cluster1.de.markers[cluster1.de.markers$p_val_adj <= 0.05 &
                                                   abs(cluster1.de.markers$avg_log2FC) >= 0.25,]
cluster1.selected.markers <- na.omit(cluster1.selected.markers)
cluster1.selected.markers
write.csv(cluster1.selected.markers,"c1-0.25.csv")
EnhancedVolcano(cluster1.selected.markers,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = rownames(cluster1.selected.markers),
                FCcutoff = 0.25,
                pCutoff = 0.05,
                title = "cluster1 differential expression(log2FC= 0.25)")
VlnPlot(cluster1,"H19",group.by = "embryo",split.by = "label")


cluster2 <- subset(cell1,subset = seurat_clusters == "2")
cluster2.de.markers <- FindMarkers(cluster2,ident.1 = "old",ident.2 = "young",only.pos = FALSE,
                                   test.use = "DESeq2",logfc.threshold = 0.5)
cluster2.de.markers <- as.data.frame(cluster2.de.markers)
cluster2.selected.markers <- cluster2.de.markers[cluster2.de.markers$p_val_adj <= 0.05 &
                                                   abs(cluster2.de.markers$avg_log2FC) >= 0.5,]
cluster2.selected.markers <- na.omit(cluster2.selected.markers)
cluster2.selected.markers
write.csv(cluster2.selected.markers,"c2-0.5.csv")
EnhancedVolcano(cluster2.de.markers,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = rownames(cluster2.de.markers),
                FCcutoff = 0.5,
                pCutoff = 0.05,
                title = "cluster2 differential expression(log2FC= 0.5)")
VlnPlot(cluster2,"RPS26",group.by = "embryo",split.by = "label")


cluster3 <- subset(cell1,subset = seurat_clusters == "3")
cluster3.de.markers <- FindMarkers(cluster3,ident.1 = "old",ident.2 = "young",only.pos = FALSE,
                                   test.use = "DESeq2",logfc.threshold = 0.5)
cluster3.de.markers <- as.data.frame(cluster3.de.markers)
cluster3.selected.markers <- cluster3.de.markers[cluster3.de.markers$p_val_adj <= 0.05 &
                                                   abs(cluster3.de.markers$avg_log2FC) >= 0.5,]
cluster3.selected.markers <- na.omit(cluster3.selected.markers)
cluster3.selected.markers
write.csv(cluster3.selected.markers,"c3-0.5.csv")
EnhancedVolcano(cluster3.de.markers,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = rownames(cluster3.de.markers),
                FCcutoff = 0.5,
                pCutoff = 0.05,
                title = "cluster3 differential expression(log2FC= 0.5)" )
VlnPlot(cluster3,"RPL10",group.by = "embryo",split.by = "label")

####### pseudobulk the counts based on label-celltype
pseudo <- AggregateExpression(cell1,assays = "RNA",return.seurat = T,group.by = c("label","embryo","seurat_clusters")) 

# the metadata for the pseudobulk object is missing, so we need to add it back
pseudo$seurat_clusters <- sapply(strsplit(Cells(pseudo), split = "_"), "[", 3)
pseudo$embryo <- sapply(strsplit(Cells(pseudo), split = "_"), "[", 2)
pseudo$label <- sapply(strsplit(Cells(pseudo), split = "_"), "[", 1)
pseudo$celltype.label <- paste(pseudo$seurat_clusters, pseudo$label, sep = "_")
Idents(pseudo) <- "celltype.label"
bulk.de <- FindMarkers(object = pseudo, 
                      ident.1 = "0_old", 
                      ident.2 = "0_young",
                      test.use = "t",
                      logfc.threshold = 0.5)
head(bulk.de, n = 10)

bulk.de <- as.data.frame(bulk.de)
bulk.de <- bulk.de[bulk.de$p_val_adj <= 0.05,]
bulk.de <- na.omit(bulk.de)
bulk.de

c2 <- subset(cell1,subset = seurat_clusters ==  "2")
VlnPlot(c2,"CAPG",group.by = "embryo",split.by = "label")
VlnPlot(c2,"SERF2",group.by = "embryo",split.by = "label")
VlnPlot(c2,"ENO1",group.by = "embryo",split.by = "label")

c3 <- subset(cell1,subset = seurat_clusters ==  "3")
VlnPlot(c3,c("CAPG","SERF2","ENO1","CTNND1"),group.by = "embryo",split.by = "label")
