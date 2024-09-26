suppressMessages({
  library(ggplot2)   # 加载ggplot2包用于绘图
  library(reshape2)  # 加载reshape2包用于数据重构
  library(Seurat)    # 加载Seurat包用于单细胞RNA测序数据分析
  library(dplyr)     # 加载dplyr包用于数据操作
  library(pheatmap)  # 加载pheatmap包用于绘制热图
  library(ramify)    # 加载ramify包用于矩阵操作
  library(readxl)    # 加载readxl包用于读取Excel文件
})



# 定义胚胎数据列表
embryo_data_list = c("1-1","1-2","1-3","1-4","2-1","2-2","2-3","2-4","2-5","3-1","3-2","3-3","3-4","5","6-1","6-2","6-3","7")

# 创建空列表用于存储数据
data_list = list()

# 读取数据并存储到 data_list 中
for (n in 1:length(embryo_data_list)) {
  data_path <- paste0("/home/xi/Desktop/2024-5-8-human_embryo/human_0113/", embryo_data_list[n])
  
  # 确保路径存在
  if (dir.exists(data_path)) {
    data_list[[n]] = Read10X(data.dir = data_path)
    
    # 检查数据是否有效
    if (is.null(data_list[[n]])) {
      print(paste("No data found for embryo:", embryo_data_list[n]))
    } else {
      print(paste("Data successfully read for embryo:", embryo_data_list[n]))
    }
  } else {
    print(paste("Directory does not exist for embryo:", embryo_data_list[n]))
  }
}

# 创建 Seurat 对象并存储到 seurat_object 列表中
seurat_object = list()  # 确保是空的列表
for (n in 1:length(embryo_data_list)) {
  if (!is.null(data_list[[n]])) {
    # 使用 Seurat 包中的 CreateSeuratObject() 函数创建 Seurat 对象
    seurat_object[[n]] = CreateSeuratObject(counts = data_list[[n]], project = "embryo_human_0112")
    
    # 确认 Seurat 对象是否成功创建
    if (inherits(seurat_object[[n]], "Seurat")) {
      print(paste("Created Seurat object for embryo:", embryo_data_list[n]))
    } else {
      print(paste("Failed to create Seurat object for embryo:", embryo_data_list[n]))
    }
  } else {
    print(paste("No data to create Seurat object for embryo:", embryo_data_list[n]))
  }
}

# 为每个 Seurat 对象添加 plate 信息
for (n in 1:18) {
  # 检查是否为 Seurat 对象
  if (inherits(seurat_object[[n]], "Seurat")) {
    # 添加 plate 信息
    seurat_object[[n]]$plate = embryo_data_list[[n]]
    print(paste("Added plate to Seurat object:", embryo_data_list[n]))
  } else {
    print(paste("Error: seurat_object[[", n, "]] is not a Seurat object."))
  }
}

# 为每个 Seurat 对象添加 barcode 信息
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

# 读取 R51 数据并确保正确格式
R51 <- readLines("r_related_and_whitelist/A1-E10.txt")
R51 <- data.frame(barcode = R51)

# 根据 R51 中的 barcode 过滤 Seurat 对象
r51 <- subset(seurat_object[[14]], barcode %in% R51$barcode)

# 获取 r51 中的细胞
r51_cells <- Cells(r51)

# 获取不在 r51 中的细胞并创建 r26
r26 <- subset(seurat_object[[14]], cells = Cells(seurat_object[[14]])[!Cells(seurat_object[[14]]) %in% r51_cells])

# 检查 r26 中的细胞数量
print(length(Cells(r26)))

# 合并 Seurat 对象以创建胚胎数据
E1 <- merge(seurat_object[[1]], y = c(seurat_object[[2]], seurat_object[[3]], seurat_object[[4]]))
E1$embryo <- "3-1"

E2 <- merge(seurat_object[[5]], y = c(seurat_object[[6]], seurat_object[[7]], seurat_object[[8]], seurat_object[[9]], r26))
E2$embryo <- "3-2"

E3 <- merge(seurat_object[[10]], y = c(seurat_object[[11]], seurat_object[[12]], seurat_object[[13]]))
E3$embryo <- "3-3"

E5 <- r51
E5$embryo <- "3-5"

E6 <- merge(seurat_object[[15]], y = c(seurat_object[[16]], seurat_object[[17]]))
E6$embryo <- "3-6"

E7 <- seurat_object[[18]]
E7$embryo <- "3-7"

# 添加其他元数据
E1$patient_age  <- "29"
E2$patient_age  <- "29"
E3$patient_age  <- "29"
E5$patient_age  <- "37"
E6$patient_age  <- "37"
E7$patient_age  <- "37"

E1$patient_id  <- "E0008"
E2$patient_id  <- "E0008"
E3$patient_id  <- "E0008"
E5$patient_id  <- "I2100411"
E6$patient_id  <- "I2100411"
E7$patient_id  <- "I2100411"

E1$SequenceTime <- "20230112"
E2$SequenceTime <- "20230112"
E3$SequenceTime <- "20230112"
E5$SequenceTime <- "20230112"
E6$SequenceTime <- "20230112"
E7$SequenceTime <- "20230112"

E1$PTG_A <- "aneuploid"
E2$PTG_A <-"euploid"
E3$PTG_A <- "euploid"
E5$PTG_A <- "aneuploid"
E6$PTG_A <- "euploid"
E7$PTG_A <- "euploid"

# 标签数据并合并为最终数据集
young_0112 <- merge(E1, y = c(E2, E3))
young_0112$label <- "young"

old_0112 <- merge(E5, y = c(E6, E7))
old_0112$label <- "old"

data_0112 <- merge(young_0112, old_0112)
data_0112$batch <- "20230113"
data_0112$PickingTime <- "2023.1.13"

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
all_metadata <- data_0112@meta.data

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
  csv_filename <- paste0("PLATE_3_", plate, ".csv")
  write.csv(result_df, file = csv_filename, row.names = FALSE)
  
  # 将处理后的 Seurat 对象保存回列表中
  seurat_object[[n]] <- current_seurat_object
}


###############################################################################