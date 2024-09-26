#This is to peform TE analysis after mapping to expression matrix

library(Seurat)
library(readxl)
library(dplyr)
library(patchwork)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ramify)
library(Polychrome)
library(ggplate)

#############
setwd("/home/xi/scTE/csv")
embryo_data_list=c("1-1-1","1-1-3","1-2-2","1-3-1","1-3-2","1-6-1","1-6-2","1-6-3")
data_list = list()
for (n in 1:8){
  data_list[[n]]=read.csv(paste0(embryo_data_list[n] ,".csv"))
  rownames(data_list[[n]]) = data_list[[n]][[1]]
  data_list[[n]] <- data_list[[n]][-1]
  data_list[[n]] <- t(data_list[[n]])
  print(n)
}
R11 <- CreateSeuratObject(counts = data_list[[1]],project = "embryo_human_1101",min.features = 100)
R13 <- CreateSeuratObject(counts = data_list[[2]],project = "embryo_human_1101",min.features = 100)
R22 <- CreateSeuratObject(counts = data_list[[3]],project = "embryo_human_1101",min.features = 100)
R31 <- CreateSeuratObject(counts = data_list[[4]],project = "embryo_human_1101",min.features = 100)
R32 <- CreateSeuratObject(counts = data_list[[5]],project = "embryo_human_1101",min.features = 100)
R61 <- CreateSeuratObject(counts = data_list[[6]],project = "embryo_human_1101",min.features = 100)
R62 <- CreateSeuratObject(counts = data_list[[7]],project = "embryo_human_1101",min.features = 100)
R63 <- CreateSeuratObject(counts = data_list[[8]],project = "embryo_human_1101",min.features = 100)

R11$plate <-"1_1_1"
R13$plate <- "1_1_3"
R22$plate <-"1_2_2_D12_4_1"
R31$plate <- "1_3_1"
R32$plate <- "1_3_2_N12"
R61$plate <- "1_6_1"
R62$plate <- "1_6_2"
R63$plate <- "1_6_3_B3_2_1"

#R22 TO R21 AND R4
R22$barcode <- colnames(R22)
R22.D12 <- read.csv("R22-D12.txt", sep="",header = FALSE)
R22_D12 <- data.frame(R22.D12)
R22_D12
R21_True <- subset(R22,barcode %in% R22_D12$V1)
R21_True



R4RAW <- read_excel ("R4.xlsx",col_names = FALSE)
R4RAW <-data.frame(R4RAW)
R4 <- subset(R22,barcode %in% R4RAW$...1)
R4


#R63 TO R2 AND R6-3
R63$barcode <- colnames(R63)
R6_3 <- read_excel("R6-3.xlsx",col_names = FALSE)
R6_3 <-data.frame(R6_3)
R63_Ture <- subset(R63 , barcode %in% R6_3$...1)
R63_Ture


R22_63 <-  read_excel("R_22.xlsx",col_names = FALSE)
R22_63 <- data.frame(R22_63)
R22_True <- subset(R63, barcode %in% R22_63$...1)
R22_True


R2 <-merge(R21_True,y = R22_True)
R2


R6 <- merge(R61,y = c(R62,R63_Ture))
R6


# merging old and new data
R3 <-merge(R31, y = R32)
R1 <-merge(R11, y = R13)

R1$embryo <- "1-1"
R2$embryo <- "1-2"
R3$embryo <- "1-3"
R4$embryo <- "1-4"
R6$embryo <- "1-6"

young <-merge(R1,y = c(R2,R3))
old <- merge(R4,y = R6)
young$label <- 'young'
old$label <- 'old'
all <- merge(young ,y = old)
all$batch <- "1101"
saveRDS(all,"TE_1101.Rdata")
###################################################################
#1207 data
embryo_data_list=c("2-1-1","2-1-2","2-1-3","2-1-4","2-1-5","2-2-1","2-2-2","2-3-1","2-4-1","2-4-2","2-4-3","2-5-2","2-6")
data_list = list()
for (n in 1:13){
  data_list[[n]]=read.csv(paste0(embryo_data_list[n] ,".csv"))
  rownames(data_list[[n]]) = data_list[[n]][[1]]
  data_list[[n]] <- data_list[[n]][-1]
  data_list[[n]] <- t(data_list[[n]])
  print(n)
}

# Create Seurat Object ,filtering @ 200
R11 <- CreateSeuratObject(counts = data_list[[1]],project = "embryo_human_1207",min.features = 100)
R12 <- CreateSeuratObject(counts = data_list[[2]],project = "embryo_human_1207",min.features = 100)
R13 <- CreateSeuratObject(counts = data_list[[3]],project = "embryo_human_1207",min.features = 100)
R14 <- CreateSeuratObject(counts = data_list[[4]],project = "embryo_human_1207",min.features = 100)
R1551 <- CreateSeuratObject(counts = data_list[[5]],project = "embryo_human_1207",min.features = 100)
R21 <- CreateSeuratObject(counts = data_list[[6]],project = "embryo_human_1207",min.features = 100)
R22 <- CreateSeuratObject(counts = data_list[[7]],project = "embryo_human_1207",min.features = 100)
R31 <- CreateSeuratObject(counts = data_list[[8]],project = "embryo_human_1207",min.features = 100)
R41 <- CreateSeuratObject(counts = data_list[[9]],project = "embryo_human_1207",min.features = 100)
R42 <- CreateSeuratObject(counts = data_list[[10]],project = "embryo_human_1207",min.features = 100)
R43 <- CreateSeuratObject(counts = data_list[[11]],project = "embryo_human_1207",min.features =100)
R52 <- CreateSeuratObject(counts = data_list[[12]],project = "embryo_human_1207",min.features = 100)
R6 <- CreateSeuratObject(counts = data_list[[13]],project = "embryo_human_1207",min.features = 100)

#
R11$LOG10UMI <- log10(R11@meta.data$nCount_RNA)
R12$LOG10UMI <- log10(R12@meta.data$nCount_RNA)
R13$LOG10UMI <- log10(R13@meta.data$nCount_RNA)
R14$LOG10UMI <- log10(R14@meta.data$nCount_RNA)
R1551$LOG10UMI <- log10(R1551@meta.data$nCount_RNA)
R42$LOG10UMI <- log10(R42@meta.data$nCount_RNA)
R43$LOG10UMI <- log10(R43@meta.data$nCount_RNA)
R52$LOG10UMI <- log10(R52@meta.data$nCount_RNA)
R6$LOG10UMI <- log10(R6@meta.data$nCount_RNA)


#
R11$plate <-"2_1_1"
R12$plate <- "2_1_2"
R13$plate <-"2_1_3"
R14$plate <- "2_1_4"
R1551$plate <- "2_1_5-5_1"
R21$plate <- "2_2_1"
R22$plate <- "2_2_2"
R31$plate <- "2_3_1"
R41$plate <- "2_4_1"
R42$plate <- "2_4_2"
R43$plate <- "2_4_3"
R52$plate <- "2_5_2"
R6$plate <- "2_6"


#
R1551$barcode <- colnames(R1551)
R_half <- read.csv("R_A1_D12", sep="",header = FALSE)
R_half <- data.frame(R_half)
R_half
R15_True <- subset(R1551,barcode %in% R_half$V1)
R15_True


R_half_2 <- read.csv("R_E1_H12", sep="",header = FALSE)
R_half_2 <- data.frame(R_half_2)
R_half_2
R51_True <- subset(R1551,barcode %in% R_half_2$V1)
R51_True

E1 <- merge(R11,y = c(R12,R13,R14,R15_True)) 
E1$embryo <- "2-1"
E2 <- merge(R21, y = R22)
E2$embryo <- "2-2"
E3 <- R31
E3$embryo <- "2-3"

E4 <- merge(R41,y = c(R42,R43))
E4$embryo <- "2-4"
E5 <- merge(R52,y = R51_True)
E5$embryo <- "2-5"
E6 <- R6
E6$embryo <- "2-6"



young <- merge(E1,y = c(E2,E3))
young$label <- "young"
old <- merge(E4,y = c(E5,E6))
old$label <- "old"


all_1207 <- merge(young,old)
all_1207$batch <- "1207"

saveRDS(all_1207,"TE_1207.Rdata")

################################################
embryo_data_list=c("3-1-1","3-1-2","3-1-3","3-1-4","3-2-1","3-2-2","3-2-3","3-2-4","3-2-5","3-3-1","3-3-2","3-3-3","3-3-4","3-3-5","3-6-1","3-6-2","3-6-3","3-7")
data_list = list()
for (n in 1:13){
  data_list[[n]]=read.csv(paste0(embryo_data_list[n] ,".csv"))
  rownames(data_list[[n]]) = data_list[[n]][[1]]
  data_list[[n]] <- data_list[[n]][-1]
  data_list[[n]] <- t(data_list[[n]])
  print(n)
}
data_list[[14]]=read.csv(paste0(embryo_data_list[14] ,".csv"))
rownames(data_list[[14]]) <- NULL
rownames(data_list[[14]]) <- data_list[[14]][,1]
data_list[[14]]<- data_list[[14]][,-1]
data_list[[14]] <- t(data_list[[14]])

#data <- read.csv("~/scTE/csv/3-3-5.csv",header = TRUE, row.names = 1)
#View(data)


for (n in 15:18){
  data_list[[n]]=read.csv(paste0(embryo_data_list[n] ,".csv"))
  rownames(data_list[[n]]) = data_list[[n]][[1]]
  data_list[[n]] <- data_list[[n]][-1]
  data_list[[n]] <- t(data_list[[n]])
  print(n)
}


seurat_object = list()
for ( n in 1:18){
  seurat_object[[n]] = CreateSeuratObject(counts = data_list[[n]],project = "embryo_human_0112",min.features = 100)
  print(n)
}


for (n in 1:18){
  seurat_object[[n]]$LOG10UMI <- log10(seurat_object[[n]]@meta.data$nCount_RNA)
  print(n)
}

for (n in 1:18 ){
  seurat_object[[n]]$plate = embryo_data_list[[n]]
  print(n)
}


#split the plate and assign the embryo number
seurat_object[[14]]$barcode <- colnames(seurat_object[[14]])
R51 <- read.csv("A1_E10",sep = "",header = FALSE)
R51 <- data.frame(R51)
R51
r51 <- subset(seurat_object[[14]],barcode %in% R51$V1)

R26 <- read.csv("E10_H12",sep = "", header = FALSE)
R26 <- data.frame(R26)
R26
r26 <- subset(seurat_object[[14]],barcode %in% R26$V1)



E1 <- merge(seurat_object[[1]],y = c(seurat_object[[2]],seurat_object[[3]],seurat_object[[4]]))
E1$embryo <- "3-1"
E2 <- merge(seurat_object[[5]], y = c(seurat_object[[6]],seurat_object[[7]],seurat_object[[8]],seurat_object[[9]],r26))
E2$embryo <- "3-2"
E3 <- merge(seurat_object[[10]], y = c(seurat_object[[11]],seurat_object[[12]],seurat_object[[13]]))
E3$embryo <- "3-3"
E5 <- r51
E5$embryo <- "3-5"
E6 <- merge(seurat_object[[15]],y = c(seurat_object[[16]],seurat_object[[17]]))
E6$embryo <- "3-6"
E7 <- seurat_object[[18]]
E7$embryo <- "3-7"


young_0112 <- merge(E1, y = c(E2,E3))
young_0112$label <- "young"
old_0112 <- merge(E5,y = c(E6,E7))
old_0112$label <- "old"
data_0112 <- merge(young_0112,old_0112)
data_0112$batch <- "0112"
saveRDS(data_0112,"TE_0112.Rdata")
###############################################################

TE_11 <- readRDS("/home/xi/scTE/csv/TE_1101.Rdata")
TE_12 <- readRDS("/home/xi/scTE/csv/TE_1207.Rdata")
TE_01 <- readRDS("/home/xi/scTE/csv/TE_0112.Rdata")


TE_11 <- JoinLayers(TE_11)
TE_12 <- JoinLayers(TE_12)
TE_01 <- JoinLayers(TE_01)
TE_01
sc <- merge(TE_11,c(TE_12,TE_01))
sc <- JoinLayers(sc)
sc
VlnPlot(sc, features = c("nCount_RNA"),log = TRUE)
VlnPlot(sc, features = "nFeature_RNA")

#sc <- readRDS("TE_1101.Rdata")
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
#sc <- SCTransform(sc)
sc <- FindVariableFeatures(sc)
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
VizDimLoadings(sc, dims = 1:2, reduction = "pca")
DimPlot(sc, reduction = "pca")
DimHeatmap(sc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(sc)
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.1)
sc <- RunUMAP(sc, dims = 1:10)
DimPlot(sc, reduction = "umap")
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

sc.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(sc, features = top20$gene) + NoLegend()

saveRDS(sc,"TE_ALL.Rdata")

VlnPlot(sc, features = "TCF7")
VlnPlot(sc, features = "SOX2")
VlnPlot(sc, features = "TFAP2C")

################################################################################
data <- readRDS("TE_ALL.Rdata")

counts <- GetAssayData(object = sc, slot = "counts")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = as.data.frame(sc@meta.data),
                              design = ~group)
#dds <- varianceStabilizingTransformation(dds)
# 执行差异表达分析
dds <- DESeq(dds)

# 获取差异表达基因
res <- results(dds)

# 获取显著差异表达基因
sig_res <- subset(res, padj < 0.005 & abs(log2FoldChange) > 1)

sig_genes <- rownames(sig_res)

# 绘制热图
DoHeatmap(object = sc,
          features = sig_genes,
          group.by = "group")

VlnPlot(object = sc,
        features = sig_genes,
        group.by = "group")


FeaturePlot(object = sc,
            features = sig_genes,
            combine = TRUE)
##################################################
data <- readRDS("embryo_integrated_allembryos_filtered_0517.Rdata")
data2 <- readRDS("/home/xi/scTE/csv/TE_ALL.Rdata")
data1 <- subset(data,subset = Age == "9")
data1$dataset <- "reference_data"
data2$dataset <- "our_data"


y = merge(data1, y = data2)
embryo.list <- SplitObject(y, split.by = "dataset")

embryo.list <- lapply(X = SplitObject(y,split.by = "dataset"), FUN = SCTransform)

embryo.anchors <- FindIntegrationAnchors(object.list = embryo.list, dims = 1:30)
embryo.integrated <- IntegrateData(anchorset = embryo.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(embryo.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
embryo.integrated <- ScaleData(embryo.integrated, verbose = FALSE)
embryo.integrated <- RunPCA(embryo.integrated, npcs = 30, verbose = FALSE)
embryo.integrated <- RunUMAP(embryo.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
DimPlot(embryo.integrated, reduction = "umap", group.by = "dataset",pt.size =0.5)
DimPlot(embryo.integrated, reduction = "umap",group.by = "seurat_clusters",pt.size =0.5,repel = TRUE,label = TRUE)
DimPlot(data1,reduction = "umap")
#DimPlot(all,reduction = "umap")