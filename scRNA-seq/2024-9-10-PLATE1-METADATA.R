suppressMessages({
  library(ggplot2)   # 加载ggplot2包用于绘图
  library(reshape2)  # 加载reshape2包用于数据重构
  library(Seurat)    # 加载Seurat包用于单细胞RNA测序数据分析
  library(dplyr)     # 加载dplyr包用于数据操作
  library(pheatmap)  # 加载pheatmap包用于绘制热图
  library(ramify)    # 加载ramify包用于矩阵操作
  library(readxl)    # 加载readxl包用于读取Excel文件
  library(readxl)
})

embryo_data_list=c("1-1","1-3","2-2","3-1","3-2","6-1","6-2","6-3")
data_list = list()
for (n in 1:8){
  data_list[[n]]=Read10X(data.dir = paste0("/home/xi/Desktop/2024-5-8-human_embryo/human_1101/",embryo_data_list[n]))
  print(n)
}
#filter cells 
R11 <- CreateSeuratObject(counts = data_list[[1]],project = "embryo_20221101",min.features = 1)
R13 <- CreateSeuratObject(counts = data_list[[2]],project = "embryo_20221101",min.features = 1)
R22 <- CreateSeuratObject(counts = data_list[[3]],project = "embryo_20221101",min.features = 1)
R31 <- CreateSeuratObject(counts = data_list[[4]],project = "embryo_20221101",min.features = 1)
R32 <- CreateSeuratObject(counts = data_list[[5]],project = "embryo_20221101",min.features = 1)
R61 <- CreateSeuratObject(counts = data_list[[6]],project = "embryo_20221101",min.features = 1)
R62 <- CreateSeuratObject(counts = data_list[[7]],project = "embryo_20221101",min.features = 1)
R63 <- CreateSeuratObject(counts = data_list[[8]],project = "embryo_20221101",min.features = 1)

R11$barcode <- colnames(R11)
R13$barcode <- colnames(R13)
R31$barcode <- colnames(R31)
R32$barcode <- colnames(R32)
R61$barcode <- colnames(R61)
R62$barcode <- colnames(R62)


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
R22.D12 <- read.csv("r_related_and_whitelist/R22_D12.txt", sep="",header = FALSE)
R22_D12 <- data.frame(R22.D12)
R22_D12
R21_True <- subset(R22,barcode %in% R22_D12$V1)
R21_True



R4RAW <- read.csv("r_related_and_whitelist/R4.txt",sep="",header = FALSE)
R4RAW <-data.frame(R4RAW)
R4 <- subset(R22,barcode %in% R4RAW$V1)
R4


#R63 TO R2 AND R6-3
R63$barcode <- colnames(R63)
R22_63 <-  read_excel("r_related_and_whitelist/R_22.xlsx",col_names = FALSE)
R22_63 <- data.frame(R22_63)
R22_True <- subset(R63, barcode %in% R22_63$...1)
R22_True
R22_True$embryo <- "1-2"
setdiff(R22_63$...1, Cells(R63))

# 获取 R22_True 中的细胞名（条形码）
R22_cells <- Cells(R22_True)
# 从 R63 中删除 R22_True 中的细胞，保留其余细胞
R63_True <- subset(R63, cells = Cells(R63)[!Cells(R63) %in% R22_cells])
# 查看新的 Seurat 对象
R63_True
# 检查有多少细胞在 R63_True 中
length(Cells(R63_True))
R63_True$embryo <- "1-6"

R63 <- merge(R22_True,R63_True)

R2 <-merge(R21_True,y = R22_True)
R2


R6 <- merge(R61,y = c(R62,R63_True))
R6


# merging old and new data
R3 <-merge(R31, y = R32)
R1 <-merge(R11, y = R13)

R1$embryo <- "1-1"
R2$embryo <- "1-2"
R3$embryo <- "1-3"
R4$embryo <- "1-4"
R6$embryo <- "1-6"

print(table(R6$embryo))



R1$patient_age <- "29"
R2$patient_age <- "29"
R3$patient_age <- "29"
R4$patient_age <- "37"
R6$patient_age <- "37"

R1$patient_id <- "E0008"
R2$patient_id <- "E0008"
R3$patient_id <- "E0008"
R4$patient_id <- "I2100411"
R6$patient_id <- "I2100411"


R1$SequenceTime <- "20221101"
R2$SequenceTime <- "20221101"
R3$SequenceTime <- "20221101"
R4$SequenceTime <- "20221101"
R6$SequenceTime <- "20221101"


R1$PTG_A <- "euploid"
R3$PTG_A <- "euploid"
R2$PTG_A <- "aneuploid"
R4$PTG_A <- "aneuploid"
R6$PTG_A <- "aneuploid"

young <-merge(R1,y = c(R2,R3))
old <- merge(R4,y = R6)
young$label <- 'young'
old$label <- 'old'
all <- merge(young ,y = old)
all$batch <- "20221101"
all$PickingTime <- "2022.11.1"

# 读取 Excel 文件，不将第一行作为列名
data <- read_excel("white_list_for_scRNA.xlsx", col_names = FALSE)

# 给列名命名为 1 到 12
colnames(data) <- as.character(1:12)
rownames(data) <- LETTERS[1:8]

# 创建字典，将行名-列名配对为键，数据值作为值
dict <- list()
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    well <- paste0(rownames(data)[i], colnames(data)[j])
    dict[[well]] <- data[i, j]
  }
}


# # 定义要提取的井位置 (A1到A12, B1到B12, C1到C12, D1到D12)
# wells <- c(
#   paste0("A", 1:12),
#   paste0("B", 1:12),
#   paste0("C", 1:12),
#   paste0("D", 1:12)
# )
# 
# 
# # 从字典中提取对应的条形码
# barcodes <- unlist(lapply(wells, function(well) dict[[well]]))
# 
# # 写入到 txt 文件，每一行一个条形码
# writeLines(barcodes, "R22_D12.txt")



# wells <- c(
#   paste0("E", 1:12),
#   paste0("F", 1:12),
#   paste0("G", 1:12),
#   paste0("H", 1:12)
# )
# # 从字典中提取对应的条形码
# barcodes <- unlist(lapply(wells, function(well) dict[[well]]))
# 
# # 写入到 txt 文件，每一行一个条形码
# writeLines(barcodes, "R4.txt")

# 从 'all' 中提取 metadata
all_metadata <- all@meta.data

# 根据 'plate' 列拆分 metadata
split_metadata <- split(all_metadata, all_metadata$plate)

# 定义包含所有 Seurat 对象和对应 plate 的列表
plates <- list(
  "R11" = "1_1_1",
  "R13" = "1_1_3",
  "R22" = "1_2_2_D12_4_1",
  "R31" = "1_3_1",
  "R32" = "1_3_2_N12",
  "R61" = "1_6_1",
  "R62" = "1_6_2",
  "R63" = "1_6_3_B3_2_1"
)

# 定义一个循环来处理每个 Seurat 对象
for (name in names(plates)) {
  plate <- plates[[name]]

  # 更新 Seurat 对象的 meta.data
  assign(name, get(name))  # 获取 Seurat 对象
  seurat_object <- get(name)
  seurat_object@meta.data <- split_metadata[[plate]]

  # 创建一个空的列表来存储结果
  result_list <- list()

  # 遍历每一行 meta.data，查找条形码
  for (i in 1:nrow(seurat_object@meta.data)) {
    barcode <- seurat_object@meta.data$barcode[i]

    # 在字典中查找与当前条形码匹配的位置（well）
    well <- names(dict)[sapply(dict, function(x) x == barcode)]

    # 如果找到匹配的井位置
    if (length(well) > 0) {
      # 将井的位置和当前行的所有元数据组合在一起
      result_list[[length(result_list) + 1]] <- c(
        well = well,  # 井的位置
        as.list(seurat_object@meta.data[i, ])  # 当前行的所有元数据，转换为列表形式
      )
    }
  }

  # 将结果列表转换为数据框
  result_df <- do.call(rbind, lapply(result_list, function(x) data.frame(t(unlist(x)), stringsAsFactors = FALSE)))

  # 确保列名保持一致
  colnames(result_df) <- c("well", colnames(seurat_object@meta.data))

  # 排序 well 列
  result_df$well_letter <- substr(result_df$well, 1, 1)  # 提取字母部分
  result_df$well_number <- as.numeric(substring(result_df$well, 2))  # 提取数字部分

  # 按照字母部分和数字部分排序
  result_df <- result_df[order(result_df$well_letter, result_df$well_number), ]

  # 删除辅助列
  result_df$well_letter <- NULL
  result_df$well_number <- NULL

  # 删除指定的列
  result_df <- result_df[, !names(result_df) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA")]

  # 打印结果数据框
  print(result_df)

  # 更新处理后的 meta.data 回到 Seurat 对象中
  seurat_object@meta.data <- result_df

  # 重新赋值给对应的 Seurat 对象
  assign(name, seurat_object)

  # 生成 CSV 文件，并以 plate 名称命名文件
  csv_filename <- paste0("PLATE_",plate, ".csv")
  write.csv(result_df, file = csv_filename, row.names = FALSE)
}
##################################################################
