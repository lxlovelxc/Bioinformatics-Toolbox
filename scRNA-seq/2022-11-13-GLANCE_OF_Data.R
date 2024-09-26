# A GLANCE OF NEW DATA IN THE PAPER


library(dplyr)
library(Seurat)
library(Matrix)
library(gdata)
library(ggplot2)
setwd("C:/Users/luoxi/Desktop")

test <- Read10X(data.dir = 'C:/Users/luoxi/Desktop/test1024/')

counts_per_cell <- Matrix::colSums(test)
genes_per_cell <- Matrix::colSums(test>0)


counts_per_cell_df <- as.data.frame(counts_per_cell)

#counts per cell for 45 and 61 each
hist(log10(counts_per_cell+1),main='counts per cell',col='cadetblue1',
     xlab = "log10(counts_per_cell)",xlim = range(0,6),ylim = range(0:500000),breaks = 30)


hist(log10(counts_per_cell),col = "cadetblue1"
     ,main='log10 of total counts(bin size: 30)',xlab = 'log10(total counts)',
     xlim = range(0,6),ylim = range(0,30000),breaks = 30)


rna <- CreateSeuratObject(counts = test,min.cells = 3, min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna,features = VariableFeatures(object = rna))
#Investigate the intrinsic dimensionality of the data using an elbow plot:
ElbowPlot(rna)

VizDimLoadings(rna,dims = 1:2,reduction = "pca")

DimHeatmap(rna,dims = 1:2,balanced = TRUE)

rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna, resolution = 0.8)
DimPlot(rna,reduction = "pca")


rna <- RunUMAP(rna, dims = 1:12)
DimPlot(rna, reduction = "umap",label = TRUE,pt.size =0.5)
