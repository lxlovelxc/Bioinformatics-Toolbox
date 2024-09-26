# 加载必要的库
suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scater)
  library(muscat)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(cowplot)
  library(UpSetR)
  library(limma)
  library(BiocParallel)
})
# 读取Seurat对象
all <- readRDS("2024-7-8-filtered-data.Rdata")
all <- subset(all, subset = embryo %in% c("1-1","1-3","2-1","3-2","3-3","3-6","3-7","4-1","4-27","5-5.2","5-5.4","5-7.1","5-7.4","5-7.5"))
# 确保meta.data包含所需的列

all$embryo <- all@meta.data$embryo
all$cluster_id <- all@meta.data$seurat_clusters
all$label <- all@meta.data$label



# 将Seurat对象转换为SingleCellExperiment对象
sce <- as.SingleCellExperiment(all)

# 移除未检测到的基因
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

# 计算每个细胞的质量控制（QC）指标
qc <- perCellQCMetrics(sce)

# 移除检测到的基因数过少或过多的细胞
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# 移除低表达基因
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# 计算并标准化库因子
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

# 使用vst方法进一步标准化
library(sctransform)
assays(sce)$vstresiduals <- vst(counts(sce), verbosity = FALSE)$y

# 准备数据以供分析
sce <- prepSCE(sce,
               kid = "cluster_id",  # subpopulation assignments
               gid = "label",      # group IDs (old/young)
               sid = "embryo",          # sample IDs (old.1, young.1, etc.)
               drop = TRUE)         # drop all other colData columns



nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids


# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))
# 计算UMAP
sce <- runUMAP(sce, pca = 20)
sce <- runTSNE(sce)
plotTSNE(sce,colour_by = "cluster_id") 
plotUMAP(sce,colour_by = "cluster_id")

# 数据聚合和降维
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))#sum
assayNames(pb)
pb_mds <- pbMDS(pb)

pb_mds

# 确保对比名称为有效R变量名
levels(pb) <- make.names(levels(pb))

# 差异表达分析
res <- pbDS(pb, verbose = FALSE)

# 检查结果
tbl <- res$table[[1]]
names(tbl)
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))


# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("old-young", levels = mm)

# # run DS analysis
bp <- MulticoreParam(workers = 8)
pbDS(pb, design = mm, contrast = contrast,method = "DESeq2", BPPARAM = bp)


# mm <- mmDS(sce, method = "dream", BPPARAM = bp)
# mm <- mmDS(sce, method = "dream",
#            n_cells = 10, n_samples = 2, min_cells = 20, BPPARAM = bp)
# 
# typeof(mm)
# mm <- as.data.frame(mm)

# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# view top 2 hits in each cluster
top2 <- bind_rows(lapply(tbl_fil, top_n, 2, p_adj.loc))
format(top2[, -ncol(top2)], digits = 2)


li <- as.data.frame(tbl_fil$"0")
li
write.csv(li,"c0_pb_muscat_DE_result.csv",sep = ",",quote = F)
# frq <- calcExprFreqs(sce, assay = "counts", th = 0)
# # one sheet per cluster
# assayNames(frq)
# # expression frequencies in each
# # sample & group; 1st cluster
# t(head(assay(frq), 5))
# 
# # big-table (wide) format; attach CPMs
# resDS(sce, res, bind = "col", cpm = TRUE)
# 
# 
# # compute expression frequencies on the fly
# resDS(sce, res, frq = TRUE)

plotExpression(sce[, sce$cluster_id == "1"],
               features = tbl_fil$"0"$gene[seq_len(1)],
               x = "sample_id", colour_by = "group_id", ncol = 1) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotExpression(sce[, sce$cluster_id == "1"],
               features = c("SNHG17","PTGES2","B3GNT2","FAM81A","PNPLA4"),
               x = "sample_id", colour_by = "group_id", ncol = 1) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pbHeatmap(sce, res, top_n = 5)
pbHeatmap(sce, res, k = "2")
pbHeatmap(sce, mm, g = "SNHG17")

mm <- mm[order(mm$X0.p_adj.glb), ]
write.csv(mm,"c0_mm_muscat_DE_result.csv",sep = ",",quote = F)
