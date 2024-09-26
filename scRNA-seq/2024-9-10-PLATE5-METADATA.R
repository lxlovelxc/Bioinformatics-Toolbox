suppressMessages({
  library(ggplot2)   # 加载ggplot2包用于绘图
  library(reshape2)  # 加载reshape2包用于数据重构
  library(Seurat)    # 加载Seurat包用于单细胞RNA测序数据分析
  library(dplyr)     # 加载dplyr包用于数据操作
  library(pheatmap)  # 加载pheatmap包用于绘制热图
  library(ramify)    # 加载ramify包用于矩阵操作
  library(readxl)    # 加载readxl包用于读取Excel文件
})


embryo_data_list <- c("5A-1","5A-2","5A-3","5A-4","5A-5","5A-6","5A-7","5A-8","5A-9","5A-10","5A-11","5A-12","5A-13","5A-14","5A-15","5A-16","5A-17","5A-18","5B-1","5B-2","5B-3","5B-4","5B-5","5B-6","5B-7","5B-8","5B-9","5B-10","5B-11")

data_list = list()

for (n in 1:29){
  data_list[[n]]=Read10X(data.dir = paste0("/home/xi/Desktop/2024-6-3-human-embryo-exp5/",embryo_data_list[n]))
  print(n)
}

seurat_object = list()

for ( n in 1:29){
  seurat_object[[n]] = CreateSeuratObject(counts = data_list[[n]],project = "embryo_2024-05-31",min.features = 1)
  print(n)
}

for (n in 1:29){
  seurat_object[[n]]$plate = embryo_data_list[[n]]
  print(n)
}


# 为每个 Seurat 对象添加 barcode 信息
for (n in 1:29) {
  # 检查是否为 Seurat 对象
  if (inherits(seurat_object[[n]], "Seurat")) {
    # 提取并添加条形码（列名）
    seurat_object[[n]]$barcode <- colnames(seurat_object[[n]])
    print(paste("Added barcode to Seurat object:", embryo_data_list[n]))
  } else {
    print(paste("Error: seurat_object[[", n, "]] is not a Seurat object."))
  }
}

exp5_1 <- seurat_object[[26]]
exp5_1$label <- "young"
exp5_1$age <- "25"
exp5_1$embryo <- "5-5.1"
exp5_1$Thaw_Time <- "2024.4.18"
exp5_1$SequenceTime <- "2024.5.31"
exp5_1$Donor <- "I2200214"
exp5_1$Age <- 25
exp5_1$ID <- "#1"
exp5_1$Embryo_Day <- "D6"
exp5_1$Embryo_Grade <- "4BB"
exp5_1$PGT_A_Sampling <- ""
exp5_1$Single_cell_Sampling <- "22PR0048-08"
exp5_1$Genotype <- "NA"
exp5_1$Chromosome_Ploidy <- "suspected euploidy"
exp5_1$PickingTime <- "2024.4.22"

exp5_2 <- merge(seurat_object[[1]],c(seurat_object[[2]],seurat_object[[3]],seurat_object[[4]],seurat_object[[10]],seurat_object[[11]],seurat_object[[12]],seurat_object[[13]],seurat_object[[21]]))
exp5_2$label <- "young"
exp5_2$age <- "25"
exp5_2$embryo <- "5-5.2"
#additional information
exp5_2$Thaw_Time <- "2024.4.18"
exp5_2$SequenceTime <- "2024.5.31"
exp5_2$Donor <- "I2200214"
exp5_2$Age <- 25
exp5_2$ID <- "#2"
exp5_2$Embryo_Day <- "D6"
exp5_2$Embryo_Grade <- "5BA"
exp5_2$PGT_A_Sampling <- "22PR0075-06"
exp5_2$Single_cell_Sampling <- "22BK0019-06"
exp5_2$Genotype <- "46,XX;"
exp5_2$Chromosome_Ploidy <- "euploid"
exp5_2$PickingTime <- "2024.4.22"

exp5_3 <-seurat_object[[28]]
exp5_3$label <- "young"
exp5_3$age <- "25"
exp5_3$embryo <- "5-5.3"
exp5_3$Thaw_Time <- "2024.4.18"
exp5_3$SequenceTime <- "2024.5.31"
exp5_3$Donor <- "I2200214"
exp5_3$Age <- 25
exp5_3$ID <- "#3"
exp5_3$Embryo_Day <- "D6"
exp5_3$Embryo_Grade <- "4BB"
exp5_3$PGT_A_Sampling <- "22BK0056-09"
exp5_3$Single_cell_Sampling <- ""
exp5_3$Genotype <- "47,XX"
exp5_3$Chromosome_Ploidy <- "aneuploidy"
exp5_3$PickingTime <- "2024.4.22"


exp5_4 <- merge(seurat_object[[16]],seurat_object[[20]])
exp5_4$label <- "young"
exp5_4$age <- "25"
exp5_4$embryo <- "5-5.4"
#additional information
exp5_4$Thaw_Time <- "2024.4.18"
exp5_4$SequenceTime <- "2024.5.31"
exp5_4$Donor <- "I2200214"
exp5_4$Age <- 25
exp5_4$ID <- "#4"
exp5_4$Embryo_Day <- "D6"
exp5_4$Embryo_Grade <- "5BA"
exp5_4$PGT_A_Sampling <- "22BK0094-12"
exp5_4$Single_cell_Sampling <- "22BK0049-11"
exp5_4$Genotype <- "46,XX;"
exp5_4$Chromosome_Ploidy <- "euploid"
exp5_4$PickingTime <- "2024.4.22"


exp5_5 <- seurat_object[[15]]
exp5_5$label <- "young"
exp5_5$age <- "25"
exp5_5$embryo <- "5-5.5"
#additional information
exp5_5$Thaw_Time <- "2024.4.18"
exp5_5$SequenceTime <- "2024.5.31"
exp5_5$Donor <- "I2200214"
exp5_5$Age <- 25
exp5_5$ID <- "#5"
exp5_5$Embryo_Day <- "D6"
exp5_5$Embryo_Grade <- "4AB"
exp5_5$PGT_A_Sampling <- "22PR0084-08"
exp5_5$Single_cell_Sampling <- "22BK0049-07"
exp5_5$Genotype <- "46,XX"
exp5_5$Chromosome_Ploidy <- "aneuploidy"
exp5_5$PickingTime <- "2024.4.22"

exp5_5

exp5 <- merge(exp5_2,c(exp5_4,exp5_5,exp5_1,exp5_3))
exp5$batch <- "20240422"
exp5

exp6_1 <- seurat_object[[27]]
exp6_1$label <- "old"
exp6_1$age <- "37"
exp6_1$embryo <- "5-6.1"
exp6_1$Thaw_Time <- "2024.4.19"
exp6_1$SequenceTime <- "2024.5.31"
exp6_1$Donor <- "I2200284"
exp6_1$Age <- 37
exp6_1$ID <- "#1"
exp6_1$Embryo_Day <- "D6"
exp6_1$Embryo_Grade <- "4BC"
exp6_1$PGT_A_Sampling <- "22BK0024-09"
exp6_1$Single_cell_Sampling <- "22BK0018-05"
exp6_1$Genotype <- "46,XY;"
exp6_1$Chromosome_Ploidy <- "euploid"
exp6_1$PickingTime <- "2024.4.28"

exp6_1


exp6_2 <- merge(seurat_object[[17]],seurat_object[[18]])
exp6_2$label <- "old"
exp6_2$age <- "41"
exp6_2$embryo <- "5-6.2"
#additional information
exp6_2$Thaw_Time <- "2024.4.19"
exp6_2$SequenceTime <- "2024.5.31"
exp6_2$Donor <- "I2200116"
exp6_2$Age <- 41
exp6_2$ID <- "#2"
exp6_2$Embryo_Day <- "D6"
exp6_2$Embryo_Grade <- "4BB"
exp6_2$PGT_A_Sampling <- "22BK0013-06"
exp6_2$Single_cell_Sampling <- "22PR0085-03"
exp6_2$Genotype <- "46,XY"
exp6_2$Chromosome_Ploidy <- "aneuploidy"
exp6_2$PickingTime <- "2024.4.28"

exp6_2

exp6_3 <- seurat_object[[25]]
exp6_3$label <- "young"
exp6_3$age <- "26"
exp6_3$embryo <- "5-6.3"
exp6_3$Thaw_Time <- "2024.4.19"
exp6_3$SequenceTime <- "2024.5.31"
exp6_3$Donor <- "I2300001"
exp6_3$Age <- 26
exp6_3$ID <- "#3"
exp6_3$Embryo_Day <- "D6"
exp6_3$Embryo_Grade <- "4BB"
exp6_3$PGT_A_Sampling <- "22PR0004-03"
exp6_3$Single_cell_Sampling <- ""
exp6_3$Genotype <- "46,XX"
exp6_3$Chromosome_Ploidy <- "aneuploidy"
exp6_3$PickingTime <- "2024.4.28"

exp6_3



exp6_4 <-seurat_object[[14]]
exp6_4$label <- "young"
exp6_4$age <- "26"
exp6_4$embryo <- "5-6.4"
#additional information
exp6_4$label <- "young"
exp6_4$Thaw_Time <- "2024.4.19"
exp6_4$SequenceTime <- "2024.5.31"
exp6_4$Donor <- "I2300001"
exp6_4$Age <- 26
exp6_4$ID <- "#4"
exp6_4$Embryo_Day <- "D6"
exp6_4$Embryo_Grade <- "5CB"
exp6_4$PGT_A_Sampling <- "22BK0013-08"
exp6_4$Single_cell_Sampling <- "22PR0049-08"
exp6_4$Genotype <- "45,XY"
exp6_4$Chromosome_Ploidy <- "aneuploidy"
exp6_4$PickingTime <- "2024.4.28"

exp6_4

exp6 <- merge(exp6_2,c(exp6_4,exp6_1,exp6_3))
exp6$batch <- "20240428"
exp6



exp7_3 <- seurat_object[[23]]
exp7_3$label <- "old"
exp7_3$embryo <- "5-7.3"
exp7_3$label <- "old"
exp7_3$Thaw_Time <- "2024.5.6"
exp7_3$SequenceTime <- "2024.5.31"
exp7_3$Donor <- "I2200026"
exp7_3$Age <- 37
exp7_3$ID <- "#3"
exp7_3$Embryo_Day <- "D6"
exp7_3$Embryo_Grade <- "4BB"
exp7_3$PGT_A_Sampling <- "22BK0043-09"
exp7_3$Single_cell_Sampling <- "22PR0048-03"
exp7_3$Genotype <- "46,XY"
exp7_3$Chromosome_Ploidy <- "euploid"
exp7_3$PickingTime <- "2024.5.10"
exp7_3

exp7_5 <- merge(seurat_object[[5]],c(seurat_object[[6]],seurat_object[[7]],seurat_object[[8]],seurat_object[[22]]))
exp7_5$label <- "old"
exp7_5$age <- "37"
exp7_5$embryo <- "5-7.5"
#additional information
exp7_5$label <- "old"
exp7_5$Thaw_Time <- "2024.5.6"
exp7_5$SequenceTime <- "2024.5.31"
exp7_5$Donor <- "I2200026"
exp7_5$Age <- 37
exp7_5$ID <- "#5"
exp7_5$Embryo_Day <- "D6"
exp7_5$Embryo_Grade <- "6BA"
exp7_5$PGT_A_Sampling <- "22PR0043-02"
exp7_5$Single_cell_Sampling <- "22BK0018-03"
exp7_5$Genotype <- "NA"
exp7_5$Chromosome_Ploidy <- "suspected euploidy"
exp7_5$PickingTime <- "2024.5.10"


exp7_1 <- merge(seurat_object[[19]],seurat_object[[9]])
exp7_1$label <- "old"
exp7_1$age <- "37"
exp7_1$embryo <- "5-7.1"
#additional information
exp7_1$label <- "old"
exp7_1$Thaw_Time <- "2024.5.6"
exp7_1$SequenceTime <- "2024.5.31"
exp7_1$Donor <- "I2200026"
exp7_1$Age <- 37
exp7_1$ID <- "#1"
exp7_1$Embryo_Day <- "D6"
exp7_1$Embryo_Grade <- "4BB"
exp7_1$PGT_A_Sampling <- "22PR0075-07"
exp7_1$Single_cell_Sampling <- "22BK0018-11"
exp7_1$Genotype <- "46,XX"
exp7_1$Chromosome_Ploidy <- "euploid"
exp7_1$PickingTime <- "2024.5.10"

exp7_1


exp7_2_barcodes <- c("AGTAAACC", "CCGTTTAG", "GACGCCGA")

all_colnames <- colnames(seurat_object[[24]])

exp7_2 <- seurat_object[[24]][, match(exp7_2_barcodes, all_colnames)]

exp7_4 <- merge(seurat_object[[29]],seurat_object[[24]][, -match(exp7_2_barcodes, all_colnames)])


exp7_2$label <- "old"
exp7_2$age <- "37"
exp7_2$embryo <- "5-7.2"
#additional information
exp7_2$label <- "old"
exp7_2$Thaw_Time <- "2024.5.6"
exp7_2$SequenceTime <- "2024.5.31"
exp7_2$Donor <- "I2200026"
exp7_2$Age <- 37
exp7_2$ID <- "#2"
exp7_2$Embryo_Day <- "D6"
exp7_2$Embryo_Grade <- "4BC"
exp7_2$PGT_A_Sampling <- ""
exp7_2$Single_cell_Sampling <- "22PR0048-06"
exp7_2$Genotype <- "NA"
exp7_2$Chromosome_Ploidy <- "suspected euploidy"
exp7_2$PickingTime <- "2024.5.10"

exp7_4$label <- "old"
exp7_4$age <- "37"
exp7_4$embryo <- "5-7.4"
#additional information
exp7_4$label <- "old"
exp7_4$Thaw_Time <- "2024.5.6"
exp7_4$SequenceTime <- "2024.5.31"
exp7_4$Donor <- "I2200026"
exp7_4$Age <- 37
exp7_4$ID <- "#4"
exp7_4$Embryo_Day <- "D6"
exp7_4$Embryo_Grade <- "5BB"
exp7_4$PGT_A_Sampling <- "22BK0016-09"
exp7_4$Single_cell_Sampling <- "22BK0018-12"
exp7_4$Genotype <- "46,XY"
exp7_4$Chromosome_Ploidy <- "euploid"
exp7_4$PickingTime <- "2024.5.10"

exp7_4


exp7 <- merge(exp7_1,c(exp7_5,exp7_2,exp7_4,exp7_3))
exp7$batch <- "20240510"
exp7
experiment5 <-merge(exp5,c(exp6,exp7))

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
all_metadata <- experiment5@meta.data

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
#########################################################################################
# # 读取Excel文件中的数据
# data <- read_excel("white_list_for_scRNA.xlsx", col_names = FALSE)
# 
# # 将列名设置为数字1到12，行名设置为字母A到H
# colnames(data) <- c(1:12)
# rownames(data) <- c("A","B","C","D","E","F","G","H")
# data$Row_Labels <- c("A", "B", "C", "D", "E", "F", "G", "H")
# data <- data[, c(ncol(data), 1:(ncol(data)-1))]
# rownames(data) <- NULL
# data <- as.data.frame(data)
# rownames(data) <- data$Row_Labels
# data <- data[, -1]

# # 处理每个plate
# 
# plates <- c("5A-1","5A-2","5A-3","5A-4","5A-5","5A-6","5A-7","5A-8","5A-9","5A-10","5A-11","5A-12","5A-13","5A-14","5A-15","5A-16","5A-17","5A-18","5B-1","5B-2","5B-3","5B-4","5B-5","5B-6","5B-7","5B-8","5B-9","5B-10","5B-11")
# 
# for (plate_value in plates) {
#   # 选择指定plate的数据
#   data_plate <- subset(experiment5 , plate == plate_value)
#   
#   # 获取data_plate的列名并截取前八个字符
#   col_names <- colnames(data_plate)
#   col_names <- substr(col_names, 1, 8)
#   col_names <- as.data.frame(col_names)
#   
#   # 初始化一个空的列表，用于存储匹配的行列
#   matching_indices <- list()
#   
#   # 遍历col_names中的每个项
#   for (item in col_names$col_name) {
#     # 查找data中与当前项匹配的行列索引
#     match_indices <- which(data == item, arr.ind = TRUE)
#     
#     # 如果找到匹配项，则将其行列名添加到匹配的列表中
#     if (length(match_indices) > 0) {
#       for (i in 1:nrow(match_indices)) {
#         matching_indices[[length(matching_indices) + 1]] <- c(row = rownames(data)[match_indices[i,1]], col = colnames(data)[match_indices[i,2]])
#       }
#     }
#   }
#   
#   # 创建一个新的数据框来存储匹配的行列名
#   matching_df <- do.call(rbind, matching_indices)
#   colnames(matching_df) <- c("col_name", "row_name")
#   
#   # 将匹配的行列名添加回到col_names中
#   col_names <- cbind(col_names, matching_df)
#   
#   # 创建place列，将col_name和row_name合并
#   col_names$place <- paste(col_names$col_name, col_names$row_name, sep = "")
#   
#   # 检查重复项
#   duplicates <- duplicated(col_names$place)
#   
#   # 打印重复项
#   print(col_names$place[duplicates])
#   
#   # 按照col_name和row_name列排序
#   col_names <- col_names[order(col_names$col_name, col_names$row_name), ]
#   
#   # 将row_name转换为因子，并指定顺序
#   col_names$row_name <- factor(col_names$row_name, levels = 1:12)
#   
#   # 重新排列数据框
#   sorted_col_names <- col_names[order(factor(col_names$col_name, levels = LETTERS[1:8]), col_names$row_name), ]
#   
#   # 设置行名并移除不需要的列
#   rownames(sorted_col_names) <- sorted_col_names$place
#   sorted_col_names <- sorted_col_names[, -c(2, 3)]
#   
#   # 获取embryo和label元数据
#   embryo_metadata <- data_plate$embryo
#   label_metadata <- data_plate$label
#   
#   Thaw_Time_metadata <- data_plate$Thaw_Time
#   SequenceTime_metadata <- data_plate$SequenceTime
#   Donor_metadata <- data_plate$Donor
#   Age_metadata <- data_plate$Age
#   ID_metadata <- data_plate$ID
#   Embryo_Day_metadata <- data_plate$Embryo_Day
#   Embryo_Grade_metadata <- data_plate$Embryo_Grade
#   PGT_A_Sampling_metadata <- data_plate$PGT_A_Sampling
#   Single_cell_Sampling_metadata <- data_plate$Single_cell_Sampling
#   Genotype_metadata <- data_plate$Genotype
#   Chromosome_Ploidy_metadata <- data_plate$Chromosome_Ploidy
#   
#   
#   # sorted_col_names <- cbind(sorted_col_names, embryo = embryo_metadata, label = label_metadata,patient_id = id_metadata,patient_age= age_metadata)
#   sorted_col_names <- cbind(
#     sorted_col_names,
#     embryo = embryo_metadata,
#     label = label_metadata,
#     Thaw_Time = Thaw_Time_metadata,
#     SequenceTime = SequenceTime_metadata,
#     Donor = Donor_metadata,
#     Age = Age_metadata,
#     ID = ID_metadata,
#     Embryo_Day = Embryo_Day_metadata,
#     Embryo_Grade = Embryo_Grade_metadata,
#     PGT_A_Sampling = PGT_A_Sampling_metadata,
#     Single_cell_Sampling = Single_cell_Sampling_metadata,
#     Genotype = Genotype_metadata,
#     Chromosome_Ploidy = Chromosome_Ploidy_metadata
#   )
#   colnames(sorted_col_names) <- c(
#     "barcode", "well", "embryo", "label", 
#     "Thaw_Time", "SequenceTime", "Donor", "Age", "ID", "Embryo_Day",
#     "Embryo_Grade", "PGT_A_Sampling", "Single_cell_Sampling",
#     "Genotype", "Chromosome_Ploidy"
#   )
#   # colnames(sorted_col_names) <- c("barcode", "well", "embryo", "label","patient_id","patient_age")  
#   
#   # 在第一列添加与well列相同的数据
#   sorted_col_names <- cbind(sorted_col_names[, "well", drop = FALSE], sorted_col_names)
#   
#   # 将排序后的列名数据写入CSV文件，文件名为plate_3_<plate_value>.csv
#   write.csv(sorted_col_names, file = paste0("plate_", gsub("-", "_", plate_value), ".csv"), row.names = FALSE)
# }
