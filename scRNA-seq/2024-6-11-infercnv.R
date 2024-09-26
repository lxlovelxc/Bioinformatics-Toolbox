library(infercnv)
library(Seurat)
#data <- readRDS("2024-6-5-exp1-5_nFeature320_mt50.Rdata")
cou_mat <- GetAssayData(data,slot = "counts")
head(cou_mat)

# 创建注释文件
label <- as.data.frame(data$label)
label$cellID <- rownames(label)
label$celltype <- label$`data$label`
label <- label[, c("cellID", "celltype")]
colnames(label) <- c("cellID", "celltype")

# 确保注释文件中的细胞名称与数据矩阵中的名称匹配
matrix_cells <- colnames(cou_mat)
annotation_cells <- label$cellID
missing_cells <- setdiff(annotation_cells, matrix_cells)

if (length(missing_cells) > 0) {
  stop("The following cells from the annotation file are not found in the data matrix: ", paste(missing_cells, collapse = ", "))
}
# 将注释文件保存为临时文件
temp_annotation_file <- tempfile(fileext = ".txt")
write.table(label, file = temp_annotation_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 创建 infercnv 对象
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = cou_mat,
                                    annotations_file = temp_annotation_file,
                                    delim = "\t",
                                    gene_order_file = "hg38_gencode_v27.txt",
                                    ref_group_names = NULL)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 1,
                             out_dir = "inferCNV_0703",
                             cluster_by_groups = T,
                             denoise = T,
                             HMM = T
)


infercnv::plot_cnv(infercnv_obj,
                   plot_chr_scale =  T,
                   output_filename = "infer_CNV1",
                   output_format = "png")
