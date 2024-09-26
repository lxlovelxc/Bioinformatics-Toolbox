library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
#GSEA

# c0 <- read.csv("embryo_batch_916-cluster0_res_tbl.csv",sep = "")

# c0 <- read.csv("embryo_batch_916-cluster1_res_tbl.csv",sep = "")
# c0 <- read.csv("donor_batch_916-cluster0_res_tbl.csv",sep = "")
# c0 <- read.csv("donor_batch_916-cluster1_res_tbl.csv",sep = "")
c0 <- read.csv("2024-9-18-cluster1_filter_nfeature_525_tbl.csv",sep = "")
# 将基因符号转换为Entrez ID
gene_list <- bitr(c0$gene, fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# 将log2FC按照Entrez ID进行排序
gene_list <- merge(gene_list, c0, by.x = "SYMBOL", by.y = "gene")
gene_list <- gene_list[order(gene_list$log2FoldChange, decreasing = TRUE), ] # 按log2FC排序
gene_list <- na.omit(gene_list) # 移除缺失值
# 创建一个命名的向量，用于GSEA分析
gene_list_vector <- gene_list$log2FoldChange
names(gene_list_vector) <- gene_list$ENTREZID
gene_list_vector <- sort(gene_list_vector, decreasing = TRUE)

gsea_kegg <- gseKEGG(geneList = gene_list_vector,
                     organism = 'hsa',
                     pvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 2000,
                     verbose = FALSE)
# gseaplot(gsea_kegg, geneSetID = 1, title = "Top Enriched Pathway")
# 绘制 KEGG 通路富集的点图
dotplot(gsea_kegg, showCategory = 10, title = "GSEA KEGG Pathways")

# 绘制 KEGG 富集得分的山脊图
ridgeplot(gsea_kegg, showCategory = 10) + ggtitle("GSEA KEGG Pathway")
####################################################################
# #KEGG
# # 获取 KEGG 基因集
# msigdbr_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
# 
# # 准备 TERM2GENE 数据框，供 GSEA 分析使用
# term2gene <- msigdbr_df %>% dplyr::select(gs_name, entrez_gene)
# # 使用 GSEA 进行 KEGG 基因集分析
# gsea_kegg <- GSEA(geneList = gene_list_vector,   # 排序的基因列表
#                   TERM2GENE = term2gene,        # KEGG 基因集
#                   pvalueCutoff = 0.2,          # P值显著性阈值
#                   minGSSize = 10,               # 基因集最小基因数
#                   maxGSSize = 2000,              # 基因集最大基因数
#                   verbose = FALSE)
# # 绘制 KEGG 通路富集的点图
# dotplot(gsea_kegg, showCategory = 10, title = "GSEA KEGG Pathways")
# 
# # # 绘制 KEGG 富集得分曲线
# # gseaplot(gsea_kegg, geneSetID = 1, title = "Top Enriched KEGG Pathway")
# 
# # 绘制 KEGG 富集得分的山脊图
# ridgeplot(gsea_kegg, showCategory = 10) + ggtitle("GSEA KEGG Pathway")
# # # 查看 GSEA 分析的显著结果
# # head(gsea_kegg)
# # 
# # # 提取显著的通路信息
# # significant_pathways <- gsea_kegg@result[gsea_kegg@result$p.adjust < 0.05, ]
# # 

####################################################################
#msigdbr
# 读取数据
# c0 <- read.csv("embryo_batch_916-cluster1_res_tbl.csv", sep = "")
# c0 <- read.csv("donor_batch_916-cluster1_res_tbl.csv",sep = "")

# 将基因符号转换为Entrez ID
gene_list <- bitr(c0$gene, fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# 合并log2FoldChange和Entrez ID
gene_list <- merge(gene_list, c0, by.x = "SYMBOL", by.y = "gene")
gene_list <- gene_list[order(gene_list$log2FoldChange, decreasing = TRUE), ]  # 按log2FC排序
gene_list <- na.omit(gene_list)  # 移除缺失值

# 创建一个命名的向量，用于GSEA分析
gene_list_vector <- gene_list$log2FoldChange
names(gene_list_vector) <- gene_list$ENTREZID

# 从msigdbr获取基因集（这里以Hallmark基因集为例）
msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")

# 准备TERM2GENE数据框
term2gene <- msigdbr_df %>% dplyr::select(gs_name, entrez_gene)

# 进行GSEA分析
gsea_result <- GSEA(geneList = gene_list_vector,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 0.2,
                    minGSSize = 10,
                    maxGSSize = 2000,
                    verbose = FALSE)

# 将结果转换为可读格式（可选）
gsea_result <- setReadable(gsea_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 绘制点图
dotplot(gsea_result, showCategory = 10, title = "GSEA MSigDB Pathways")
ridgeplot(gsea_result, showCategory = 10)+ ggtitle("GSEA MSigDB Pathways")

# # 绘制GSEA的富集得分曲线
# gseaplot(gsea_result, geneSetID = 1, title = "Top Enriched Pathway")
################################################################################
#Gene Ontology 
library(clusterProfiler)
library(org.Hs.eg.db)

# 创建命名的向量
gene_list_vector <- gene_list$log2FoldChange
names(gene_list_vector) <- gene_list$ENTREZID

# 进行 GSEA 分析
gsea_go <- gseGO(geneList = gene_list_vector,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",  # 可以是 "BP", "CC", "MF" 或 "ALL"
                 pvalueCutoff = 0.2,
                 minGSSize = 10,
                 maxGSSize = 2000,
                 verbose = FALSE)

# 绘制结果
dotplot(gsea_go, showCategory = 10, title = "GSEA GO Terms")
ridgeplot(gsea_go, showCategory = 10) + ggtitle("GSEA GO Terms")
# ridgeplot(gsea_go, showCategory = 3) + ggtitle("GSEA GO Terms")
###################################################################################
library(ReactomePA)
library(org.Hs.eg.db)

# 使用 gene_list_vector

# 进行 GSEA 分析
gsea_reactome <- gsePathway(geneList = gene_list_vector,
                            organism = "human",
                            pvalueCutoff = 0.2,
                            minGSSize = 10,
                            maxGSSize = 2000,
                            verbose = FALSE)

# 绘制结果
dotplot(gsea_reactome, showCategory = 10, title = "GSEA Reactome Pathways")
ridgeplot(gsea_reactome, showCategory = 10)+ ggtitle("GSEA Reactome Pathways")
# ridgeplot(gsea_reactome, showCategory = 6)+ ggtitle("GSEA Reactome Pathways")
###############################################################################
# # 计算带符号的统计量
# gene_list$statistic <- -log10(gene_list$pvalue) * sign(gene_list$log2FoldChange)
# gene_list_vector <- gene_list$statistic
# names(gene_list_vector) <- gene_list$ENTREZID
# gene_list_vector <- sort(gene_list_vector, decreasing = TRUE)
# 
# # 进行GSEA
# gsea_result <- GSEA(geneList = gene_list_vector, 
#                     TERM2GENE = term2gene, 
#                     pvalueCutoff = 0.05,
#                     minGSSize = 10,
#                     maxGSSize = 500,
#                     verbose = FALSE)
# ridgeplot(gsea_reactome, showCategory = 10)