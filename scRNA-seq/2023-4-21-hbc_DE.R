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
all$group_id <- all$label
all$sample_id <- all$embryo
head(all)

counts <- all@assays$RNA@counts 
metadata <- all@meta.data
metadata$cluster_id <- factor(all@active.ident)
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
colnames(pb[[1]] ) <- c("1-1","1-2","1-3","1-4","1-6","2-1","2-2","2-4","2-5","2-6","3-1","3-2","3-3","3-5","3-6","3-7")
colnames(pb[[2]] ) <- levels(sce$sample_id)
colnames(pb[[3]] ) <- c("1-1","1-2","1-4","1-6","2-2","2-4","2-6","3-2","3-3","3-5","3-7")
colnames(pb[[4]] ) <- c("1-2","1-3","1-6","2-2","2-5","2-6","3-3","3-5","3-7")
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

a <- c(rep("exp1 ",each = 5),rep("exp2 ",each = 5),rep("exp3 ",each = 6),rep("exp1 ",each = 5),rep("exp2 ",each = 6),rep("exp3 ",each = 6),
            rep("exp1 ",each = 4),rep("exp2 ",each = 3),rep("exp3 ",each = 4),rep("exp1 ",each = 3),rep("exp2 ",each = 3),rep("exp3 ",each = 3))

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
                              design = ~ group_id + experiment) #+experiment
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])


# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Output results of Wald test for contrast for stim vs ctrl
#levels(as.factor(cluster_metadata$group_id))[2]# young
#levels(as.factor(cluster_metadata$group_id))[1]# old

#contrast <- c("condition", "level_to_compare", "base_level")
#contrast <- c("group_id", levels(as.factor(cluster_metadata$group_id))[2], levels(as.factor(cluster_metadata$group_id))[1])

#contrast <- c("group_id", levels(as.factor(cluster_metadata$group_id))[1], levels(as.factor(cluster_metadata$group_id))[2])
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


# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res
cluster1 <- subset(all,subset = seurat_clusters == "0")
VlnPlot(cluster1,"CRIP1",group.by ="embryo" ,split.by = "label")
# Write significant results to file
write.csv(sig_res,
          paste0("results_0922", clusters[2], "_", "old_vs_young", "_sig_genes_consider_experiment.csv"),
          quote = FALSE, 
          row.names = FALSE)

## ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)



###top40
top40_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)

VlnPlot(all,top40_sig_genes,group.by = "embryo",split.by = "label")
top40_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top40_sig_genes)
colnames(top40_sig_norm) <- c("gene","1-1","1-2","1-3","1-4","1-6","2-1","2-2","2-3","2-4","2-5","2-6","3-1","3-2","3-3","3-5","3-6","3-7")

gathered_top40_sig <- top40_sig_norm %>%
  gather(colnames(top40_sig_norm)[2:length(colnames(top40_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top40_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top40_sig, by = c("sample_id" = "samplename"))

## plot using ggplot2
ggplot(gathered_top40_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = group_id), 
             position=position_jitter(w=0.1,h=0))+
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 40 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

#ggplot(sce) +
#  geom_point(aes(x = rownames(sce[top40_sig_genes,]),
#                 y = counts(sce[top40_sig_genes,])))

#plot(counts(sce[top40_sig_genes,]))



# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)
colnames(sig_norm) <- c("gene","1-1","1-2","1-3","1-4","1-6","2-1","2-2","2-3","2-4","2-5","2-6","3-1","3-2","3-3","3-5","3-6","3-7")
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata[, c("group_id")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of cluster 1 relative to control(consider experiment)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,10)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
