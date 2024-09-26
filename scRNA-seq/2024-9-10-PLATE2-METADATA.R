suppressMessages({
  library(ggplot2)   # 加载ggplot2包用于绘图
  library(reshape2)  # 加载reshape2包用于数据重构
  library(Seurat)    # 加载Seurat包用于单细胞RNA测序数据分析
  library(dplyr)     # 加载dplyr包用于数据操作
  library(pheatmap)  # 加载pheatmap包用于绘制热图
  library(ramify)    # 加载ramify包用于矩阵操作
  library(readxl)    # 加载readxl包用于读取Excel文件
})

embryo_data_list=c("1-1","1-2","1-3","1-4","1-5","2-1","2-2","3-1","4-1","4-2","4-3","5-2","6")
data_list = list()
for (n in 1:13){
  data_list[[n]]=Read10X(data.dir = paste0("/home/xi/Desktop/2024-5-8-human_embryo/human_1207/",embryo_data_list[n]))
  print(n)
}

# Create Seurat Object ,filtering @ 200
R11 <- CreateSeuratObject(counts = data_list[[1]],project = "embryo_20221207",min.features = 1)
R12 <- CreateSeuratObject(counts = data_list[[2]],project = "embryo_20221207",min.features = 1)
R13 <- CreateSeuratObject(counts = data_list[[3]],project = "embryo_20221207",min.features = 1)
R14 <- CreateSeuratObject(counts = data_list[[4]],project = "embryo_20221207",min.features = 1)
R1551 <- CreateSeuratObject(counts = data_list[[5]],project = "embryo_20221207",min.features = 1)
R21 <- CreateSeuratObject(counts = data_list[[6]],project = "embryo_20221207",min.features = 1)
R22 <- CreateSeuratObject(counts = data_list[[7]],project = "embryo_20221207",min.features = 1)
R31 <- CreateSeuratObject(counts = data_list[[8]],project = "embryo_20221207",min.features = 1)
R41 <- CreateSeuratObject(counts = data_list[[9]],project = "embryo_20221207",min.features = 1)
R42 <- CreateSeuratObject(counts = data_list[[10]],project = "embryo_20221207",min.features = 1)
R43 <- CreateSeuratObject(counts = data_list[[11]],project = "embryo_20221207",min.features =1)
R52 <- CreateSeuratObject(counts = data_list[[12]],project = "embryo_20221207",min.features = 1)
R6 <- CreateSeuratObject(counts = data_list[[13]],project = "embryo_20221207",min.features = 1)

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

R11$barcode <- colnames(R11)
R12$barcode <- colnames(R12)
R13$barcode <- colnames(R13)
R14$barcode <- colnames(R14)
R21$barcode <- colnames(R21)
R22$barcode <- colnames(R22)
R31$barcode <- colnames(R31)
R41$barcode <- colnames(R41)
R42$barcode <- colnames(R42)
R43$barcode <- colnames(R43)
R6$barcode <- colnames(R6)
R52$barcode <- colnames(R52)


#
R1551$barcode <- colnames(R1551)
R_half <- read.csv("r_related_and_whitelist/R22_D12.txt", sep="",header = FALSE)
R_half <- data.frame(R_half)
R_half
R15_True <- subset(R1551,barcode %in% R_half$V1)
R15_True


R_half_2 <- read.csv("r_related_and_whitelist/R4.txt",sep="",header = FALSE)
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


E1$patient_age  <- "29"
E2$patient_age  <- "29"
E3$patient_age  <- "29"
E4$patient_age  <- "37"
E5$patient_age  <- "37"
E6$patient_age  <- "37"


E1$patient_id  <- "E0008"
E2$patient_id  <- "E0008"
E3$patient_id  <- "E0008"
E4$patient_id  <- "I2100411"
E5$patient_id  <- "I2100411"
E6$patient_id  <- "I2100411"


E1$SequenceTime <- "20221207"
E2$SequenceTime <- "20221207"
E3$SequenceTime <- "20221207"
E4$SequenceTime <- "20221207"
E5$SequenceTime <- "20221207"
E6$SequenceTime <- "20221207"


E1$PTG_A <- "euploid"
E2$PTG_A <- "aneuploid"
E3$PTG_A <- "aneuploid"
E4$PTG_A <- "aneuploid"
E5$PTG_A <- "aneuploid"
E6$PTG_A <- "aneuploid"

young <- merge(E1,y = c(E2,E3))
young$label <- "young"
old <- merge(E4,y = c(E5,E6))
old$label <- "old"


all_1207 <- merge(young,old)
all_1207$batch <- "20221206"
all_1207$PickingTime <-"2022.12.6"

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
all_metadata <- all_1207@meta.data

# 根据 'plate' 列拆分 metadata
split_metadata <- split(all_metadata, all_metadata$plate)

# 将 embryo_data_list 作为 plates 列表
plates <- list(
  "R11" = "2_1_1",
  "R12" = "2_1_2",
  "R13" = "2_1_3",
  "R14" = "2_1_4",
  "R1551" = "2_1_5-5_1",
  "R21" = "2_2_1",
  "R22" = "2_2_2",
  "R31" = "2_3_1",
  "R41" = "2_4_1",
  "R42" = "2_4_2",
  "R43" = "2_4_3",
  "R52" = "2_5_2",
  "R6" = "2_6"
)

# 打印 plates 列表以检查其内容
print(plates)

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


###############################################################################