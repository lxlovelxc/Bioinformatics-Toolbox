suppressMessages({
  library(ggplot2)   # 加载ggplot2包用于绘图
  library(reshape2)  # 加载reshape2包用于数据重构
  library(Seurat)    # 加载Seurat包用于单细胞RNA测序数据分析
  library(dplyr)     # 加载dplyr包用于数据操作
  library(pheatmap)  # 加载pheatmap包用于绘制热图
  library(ramify)    # 加载ramify包用于矩阵操作
  library(readxl)    # 加载readxl包用于读取Excel文件
})

# embryo_data_list=c("4-1_1","4-1_2","4-27","4-3_1","4-3_2")
embryo_data_list=c("4-1_1","4-1_2","4-2_1","4-27","4-3_1","4-3_2","4-3_3")

data_list = list()

for (n in 1:7){
  data_list[[n]]=Read10X(data.dir = paste0("/home/xi/Desktop/2024-5-8-human_embryo/human_0503/",embryo_data_list[n]))
  print(n)
}
seurat_object = list()
for ( n in 1:7){
  seurat_object[[n]] = CreateSeuratObject(counts = data_list[[n]],project = "embryo_20240503")
  print(n)
}
for (n in 1:length(seurat_object)) {
  # 检查是否为 Seurat 对象
  if (inherits(seurat_object[[n]], "Seurat")) {
    # 提取并添加条形码（列名）
    seurat_object[[n]]$barcode <- colnames(seurat_object[[n]])
    print(paste("Added barcode to Seurat object:", embryo_data_list[n]))
  } else {
    print(paste("Error: seurat_object[[", n, "]] is not a Seurat object."))
  }
}

seurat_object[[1]]$plate <- "4-1_1"
seurat_object[[2]]$plate <- "4-1_2"
seurat_object[[3]]$plate <- "4-2_1"
seurat_object[[4]]$plate <- "4-27"
seurat_object[[5]]$plate <- "4-3_1"
seurat_object[[6]]$plate <- "4-3_2"
seurat_object[[7]]$plate <- "4-3_3"

E1 <- merge(seurat_object[[1]],seurat_object[[2]])
E1$embryo <- "4-1"
E27 <- seurat_object[[4]]
E27$embryo <- "4-27"
E2 <- seurat_object[[3]]
E2$embryo <- "4-2"
E3 <- merge(seurat_object[[6]],c(seurat_object[[5]],seurat_object[[7]]))
E3$embryo <- "4-3"


E27$patient_id <- "I2100059"
E27$patient_age <- "29"

E1$patient_id <- "I2100512"
E1$patient_age <- "25"

E3$patient_id <- "I2100512"
E3$patient_age <- "25"
E1$label <- "young"
E27$label <- "young"


E3$label <- "young"
# 为E1赋值
E1$patient_id <- "I2100512"
E1$patient_age <- "25"
E1$label <- "young"
E1$Thaw_Time <- "2024.4.11"
E1$SequenceTime <- "2024.5.3"
E1$Donor <- "I2100512"
E1$Age <- 25
E1$ID <- "#1"
E1$Embryo_Day <- "D6"
E1$Embryo_Grade <- "4BA"
E1$PGT_A_Sampling <- "22PR0085-05"
E1$Single_cell_Sampling <- "22PR0048-02"
E1$Genotype <- "46,XX;"
E1$Chromosome_Ploidy <- "Euploidy"
E1$PTG_A <- "euploid"
E1$PickingTime <- "2024.4.15"

# 为E27赋值
E27$patient_id <- "I2100059"
E27$patient_age <- "29"
E27$label <- "young"
E27$Thaw_Time <- "2024.3.11"
E27$SequenceTime <- "2024.5.3"
E27$Donor <- "I2100059"
E27$Age <- 29
E27$ID <- "#1"
E27$Embryo_Day <- "D6"
E27$Embryo_Grade <- "4BB"
E27$PGT_A_Sampling <- "22PR0094-09"
E27$Single_cell_Sampling <- "22BK0023-06"
E27$Genotype <- "46,XX;"
E27$Chromosome_Ploidy <- "Euploidy"
E27$PTG_A <- "euploid"
E27$PickingTime <- "2024.3.15"

# 为E3赋值
E3$patient_id <- "I2100512"
E3$patient_age <- "25"
E3$label <- "young"
E3$Thaw_Time <- "2024.4.11"
E3$SequenceTime <- "2024.5.3"
E3$Donor <- "I2100512"
E3$Age <- 25
E3$ID <- "#3"
E3$Embryo_Day <- "D6"
E3$Embryo_Grade <- "5BB"
E3$PGT_A_Sampling <- "22PR0050-06"
E3$Single_cell_Sampling <- "22BK0049-09"
E3$Genotype <- "46,XY;"
E3$Chromosome_Ploidy <- "Euploidy"
E3$PTG_A <- "euploid"
E3$PickingTime <- "2024.4.15"



E2$patient_id <- "I2100512"
E2$patient_age <- "25"
E2$label <- "young"
E2$Thaw_Time <- "2024.4.11"
E2$SequenceTime <- "2024.5.3"
E2$Donor <- "I2100512"
E2$Age <- 25
E2$ID <- "#2"
E2$Embryo_Day <- "D6"
E2$Embryo_Grade <- "3CB"
E2$PGT_A_Sampling <- "22BK0018-02"
E2$Single_cell_Sampling <- ""
E2$Genotype <- "46,XY;"
E2$Chromosome_Ploidy <- "aneuploid"
E2$PTG_A <- "aneuploid"
E2$PickingTime <- "2024.4.15"

E1$batch <- "20240415"
E2$batch <- "20240415"
E3$batch <- "20240415"
E27$batch <- "20240315"
all_0503 <- merge(E1,c(E27,E3,E2))


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

# 从 'all' 中提取 metadata
all_metadata <- all_0503@meta.data

# 根据 'plate' 列拆分 metadata
split_metadata <- split(all_metadata, all_metadata$plate)

# 将 embryo_data_list 作为 plates 列表
plates <- setNames(embryo_data_list, embryo_data_list)

# 打印 plates 列表以检查其内容
print(plates)

# 定义一个循环来处理每个 Seurat 对象
for (n in 1:length(plates)) {
  
  plate <- plates[[n]]  # 使用 plates 中的 plate 名称
  
  # 获取当前 Seurat 对象，并使用不同变量名避免覆盖 seurat_object 列表
  current_seurat_object <- seurat_object[[n]]
  
  # 更新当前 Seurat 对象的 meta.data
  current_seurat_object@meta.data <- split_metadata[[plate]]
  
  # 创建一个空的列表来存储结果
  result_list <- list()
  
  
  # 遍历每一行 meta.data，查找条形码
  for (i in 1:nrow(current_seurat_object@meta.data)) {
    barcode <- current_seurat_object@meta.data$barcode[i]
    
    # 在字典中查找与当前条形码匹配的位置（well）
    well <- names(dict)[sapply(dict, function(x) x == barcode)]
    
    # 如果找到匹配的井位置
    if (length(well) > 0) {
      # 将井的位置和当前行的所有元数据组合在一起
      result_list[[length(result_list) + 1]] <- c(
        well = well,  # 井的位置
        as.list(current_seurat_object@meta.data[i, ])  # 当前行的所有元数据，转换为列表形式
      )
    }
  }
  
  # 将结果列表转换为数据框
  result_df <- do.call(rbind, lapply(result_list, function(x) data.frame(t(unlist(x)), stringsAsFactors = FALSE)))
  
  # 确保列名保持一致
  colnames(result_df) <- c("well", colnames(current_seurat_object@meta.data))
  
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
  
  # 更新处理后的 meta.data 回到当前 Seurat 对象中
  current_seurat_object@meta.data <- result_df
  
  # 生成 CSV 文件，并以 plate 名称命名文件
  csv_filename <- paste0("PLATE_", plate, ".csv")
  write.csv(result_df, file = csv_filename, row.names = FALSE)
  
  # 将处理后的 Seurat 对象保存回列表中
  seurat_object[[n]] <- current_seurat_object
}