suppressMessages({
  library(scater)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(Matrix.utils)
  library(edgeR)
  library(dplyr)
  library(magrittr)
  library(Matrix)
  library(purrr)
  library(reshape2)
  library(S4Vectors)
  library(tibble)
  library(SingleCellExperiment)
  library(pheatmap)
  library(apeglm)
  library(png)
  library(DESeq2)
  library(RColorBrewer)
  library(viridis)
  library(EnhancedVolcano)
})
data <- readRDS("2024-9-13-with-additional-embryo-raw-data.Rdata")
data1 <- subset(data, subset = nFeature_RNA > 600 )#& percent.mt < 60
plot(data1$percent.mt,data1$nFeature_RNA,pch = 20)
data1
data1 <- JoinLayers(data1)
data1 <- NormalizeData(data1)
data1 <- FindVariableFeatures(data1, selection.method = "mean.var.plot")
# all.genes <- rownames(data1)
data1 <- ScaleData(data1,model.use = 'linear',do.scale = T,do.center = T)
data1 <- RunPCA(data1, features = NULL)
ElbowPlot(data1)
data1 <- FindNeighbors(data1, dims = 1:15)#20 /25 is the best choice
data1 <- FindClusters(data1, resolution = 0.1)
data1 <- RunUMAP(data1,reduction = "pca",dims = 1:15)
DimPlot(object = data1, reduction = "umap",label= T)

all <- data1
# all <- subset(all, subset = embryo %in% c("1-1","1-3","2-1","3-2","3-3","3-6","3-7","4-1","4-27","5-5.2","5-5.4","5-7.1","5-7.4","5-7.5"))
all$group_id <- all$label
all$sample_id <- all$embryo

#counts <- all@assays$RNA@counts
counts <- GetAssayData(all,slot = "counts")
metadata <- all@meta.data
metadata$cluster_id <- factor(all@active.ident)
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

groups <- colData(sce)[, c("cluster_id", "sample_id")]

assays(sce)
dim(counts(sce))
counts(sce)[1:6,1:6]
dim(colData(sce))
head(colData(sce))


kids <- purrr::set_names(levels(sce$cluster_id))
kids
nk <- length(kids)
nk
sce$sample_id <- as.factor(sce$sample_id)
sids <- purrr::set_names(levels(sce$sample_id))
ns <- length(sids)
ns


n_cells <- as.numeric(table(sce$sample_id))
m <- match(sids, sce$sample_id)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei


groups <- colData(sce)[, c("cluster_id", "sample_id")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
class(pb)
dim(pb)
pb[1:6, 1:6]

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb)
str(pb)

table(sce$cluster_id, sce$sample_id)
colnames(pb[[1]] ) <- levels(sce$sample_id)#cluster0
colnames(pb[[2]] ) <- c("1-1", "1-2", "1-3", "1-4", "1-6", "2-1", "2-2", "2-3", "2-4", "2-5", "2-6", "3-1", "3-2", "3-3", "3-5", "3-6", "3-7", "4-1",
                        "4-27", "5-5.2", "5-5.4", "5-5.5", "5-6.2", "5-6.4", "5-7.1", "5-7.4", "5-7.5")#cluster1
##########################################################################
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}
de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()
# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 
table(all$batch)

metadata
#for all embryos
a <- c(rep("exp1",each = 5),rep("exp2",each = 6),rep("exp3",each = 6),rep("exp4",each = 2),rep("exp5",each= 10),rep("exp1",each = 5),rep("exp2",each = 6),rep("exp3",each = 6),rep("exp4",each = 2),rep("exp5",each=8),
       rep("NA",each = 33))
metadata$experiment <- a
metadata$cluster_id <- as.factor(metadata$cluster_id)
clusters <- levels(metadata$cluster_id)
clusters
############################################################################
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[2]), ]
head(cluster_metadata)
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)
counts <- pb[[clusters[2]]]

cluster_counts <- data.frame(counts[,which(colnames(counts ) %in% rownames(cluster_metadata))])

colnames(cluster_counts) <- rownames(cluster_metadata)

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ experiment + group_id)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA

DESeq2::plotPCA(rld, intgroup = c("group_id","experiment"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id","experiment"), drop=F])


# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)
contrast <- c("group_id","young","old")
# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 coef = 2,
                 type="apeglm")

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl
# write.table(res_tbl,"2024-9-12-cluster0_res_tbl.csv")
# write.table(res_tbl,"2024-9-12-cluster1_res_tbl.csv")
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

# c0_DE <- sig_res
c1_DE <- sig_res
