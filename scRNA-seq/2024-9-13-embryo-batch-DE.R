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



R1$Age <- "29"
R2$Age <- "29"
R3$Age <- "29"
R4$Age <- "37"
R6$Age <- "37"

R1$Donor <- "E0008"
R2$Donor <- "E0008"
R3$Donor <- "E0008"
R4$Donor <- "I2100411"
R6$Donor <- "I2100411"


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
##################################################################################################

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


E1$Age  <- "29"
E2$Age  <- "29"
E3$Age  <- "29"
E4$Age  <- "37"
E5$Age  <- "37"
E6$Age  <- "37"


E1$Donor  <- "E0008"
E2$Donor  <- "E0008"
E3$Donor  <- "E0008"
E4$Donor  <- "I2100411"
E5$Donor  <- "I2100411"
E6$Donor  <- "I2100411"


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
################################################################################
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
E1$Age  <- "29"
E2$Age  <- "29"
E3$Age  <- "29"
E5$Age  <- "37"
E6$Age  <- "37"
E7$Age  <- "37"

E1$Donor  <- "E0008"
E2$Donor  <- "E0008"
E3$Donor  <- "E0008"
E5$Donor  <- "I2100411"
E6$Donor  <- "I2100411"
E7$Donor  <- "I2100411"

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
###############################################################################
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
E1$label <- "young"
E1$Thaw_Time <- "2024.4.11"
E1$SequenceTime <- "20240503"
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
E27$label <- "young"
E27$Thaw_Time <- "2024.3.11"
E27$SequenceTime <- "20240503"
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
E3$label <- "young"
E3$Thaw_Time <- "2024.4.11"
E3$SequenceTime <- "20240503"
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


E2$label <- "young"
E2$Thaw_Time <- "2024.4.11"
E2$SequenceTime <- "20240503"
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
all_0503 <- E1#E27,E3

###############################################################################
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
exp5_1$embryo <- "5-5.1"
exp5_1$Thaw_Time <- "2024.4.18"
exp5_1$SequenceTime <- "20240531"
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
exp5_2$embryo <- "5-5.2"
#additional information
exp5_2$Thaw_Time <- "2024.4.18"
exp5_2$SequenceTime <- "20240531"
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
exp5_3$embryo <- "5-5.3"
exp5_3$Thaw_Time <- "2024.4.18"
exp5_3$SequenceTime <- "20240531"
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
exp5_4$embryo <- "5-5.4"
#additional information
exp5_4$Thaw_Time <- "2024.4.18"
exp5_4$SequenceTime <- "20240531"
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
exp5_5$embryo <- "5-5.5"
#additional information
exp5_5$Thaw_Time <- "2024.4.18"
exp5_5$SequenceTime <- "20240531"
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

exp5 <- merge(exp5_2,c(exp5_4,exp5_5))#,exp5_1,exp5_3
exp5$batch <- "20240422"
exp5

exp6_1 <- seurat_object[[27]]
exp6_1$label <- "old"
exp6_1$embryo <- "5-6.1"
exp6_1$Thaw_Time <- "2024.4.19"
exp6_1$SequenceTime <- "20240531"
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
exp6_2$embryo <- "5-6.2"
#additional information
exp6_2$Thaw_Time <- "2024.4.19"
exp6_2$SequenceTime <- "20240531"
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
exp6_3$embryo <- "5-6.3"
exp6_3$Thaw_Time <- "2024.4.19"
exp6_3$SequenceTime <- "20240531"
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
exp6_4$embryo <- "5-6.4"
#additional information
exp6_4$label <- "young"
exp6_4$Thaw_Time <- "2024.4.19"
exp6_4$SequenceTime <- "20240531"
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

exp6 <- merge(exp6_2,c(exp6_4,exp6_1))#,exp6_3
exp6$batch <- "20240428"
exp6



exp7_3 <- seurat_object[[23]]
exp7_3$label <- "old"
exp7_3$embryo <- "5-7.3"
exp7_3$label <- "old"
exp7_3$Thaw_Time <- "2024.5.6"
exp7_3$SequenceTime <- "20240531"
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
exp7_5$embryo <- "5-7.5"
#additional information
exp7_5$label <- "old"
exp7_5$Thaw_Time <- "2024.5.6"
exp7_5$SequenceTime <- "20240531"
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
exp7_1$embryo <- "5-7.1"
#additional information
exp7_1$label <- "old"
exp7_1$Thaw_Time <- "2024.5.6"
exp7_1$SequenceTime <- "20240531"
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
exp7_2$embryo <- "5-7.2"
#additional information
exp7_2$Thaw_Time <- "2024.5.6"
exp7_2$SequenceTime <- "20240531"
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
exp7_4$embryo <- "5-7.4"
#additional information
exp7_4$Thaw_Time <- "2024.5.6"
exp7_4$SequenceTime <- "20240531"
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


exp7 <- merge(exp7_1,c(exp7_5,exp7_2,exp7_4))#,exp7_3
exp7$batch <- "20240510"
exp7
experiment5 <-merge(exp5,c(exp6,exp7))

###############################################################################
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^MT-")
all_1207[["percent.mt"]] <- PercentageFeatureSet(all_1207, pattern = "^MT-")
data_0112[["percent.mt"]] <- PercentageFeatureSet(data_0112, pattern = "^MT-")
all_0503[["percent.mt"]] <- PercentageFeatureSet(all_0503, pattern = "^MT-")
E27[["percent.mt"]] <- PercentageFeatureSet(E27, pattern = "^MT-")
experiment5[["percent.mt"]] <- PercentageFeatureSet(experiment5, pattern = "^MT-")

data <-merge(all,c(all_1207,all_0503,data_0112,experiment5,E27))

data
# saveRDS(data,"2024-9-13-with-additional-embryo-raw-data.Rdata")
# data <- readRDS("2024-9-13-with-additional-embryo-raw-data.Rdata")
data1 <- subset(data, subset = nFeature_RNA > 600 )#& percent.mt < 25
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
FeaturePlot(data1,"nFeature_RNA")
FeaturePlot(data1,"percent.mt")
VlnPlot(data1,"nFeature_RNA",log = T)
# FeaturePlot(data1,"percent.mt")
# table(data1$seurat_clusters)
# FeaturePlot(data1,c("PEG10","TPM1","FABP5","KRT19"))
# FeaturePlot(data1,c("ATF3","PRR9","LGALS16","ANXA1"))
# FeaturePlot(data1,"percent.mt")
# DimPlot(data1,reduction = "umap",group.by = "batch")
# DimPlot(data1,reduction = "umap",group.by = "SequenceTime")
# DimPlot(data1,reduction = "umap",group.by = "Donor")
# DimPlot(data1,reduction = "umap",group.by = "Age")

#################################################################################
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
############################################################################
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
a <- c(rep("exp1",each = 5),rep("exp2",each = 6),rep("exp3",each = 6),rep("exp5",each = 1),rep("exp4",each= 1),rep("exp6",each=3),rep("exp7",each = 3)
       ,rep("exp8",each = 4),rep("exp1",each = 5),rep("exp2",each = 6),rep("exp3",each = 2),rep("exp5",each = 1),rep("exp4",each=1),rep("exp6",each = 3),
       rep("exp7",each = 2),rep("exp8",each = 3),rep("NA",each = 37))
metadata$experiment <- a
#metadata$byexperiment <- b
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
resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 coef = "group_id_young_vs_old",
                 type="apeglm")

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl
# write.table(res_tbl,"embryo_batch_916-cluster0_res_tbl.csv")
# write.table(res_tbl,"embryo_batch_916-cluster1_res_tbl.csv")
# write.table(res_tbl,"2024-9-12-cluster1_res_tbl.csv")
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res
# c1_res <-sig_res
# write.table(sig_res[sig_res$log2FoldChange < 0,]$gene,col.names = F,row.names = F,quote = F)
# write.table(sig_res[sig_res$log2FoldChange < 0,]$gene,col.names = F,row.names = F,quote = F) 
# write.table(c0_res[c0_res$log2FoldChange > 0,]$gene,col.names = F,row.names = F,quote = F) 
# c0_res <- sig_res
# write.csv(sig_res,"2024-9-13-cluster0_embryo_batch-mt-none.csv")

# write.table(twofive,"2024-9-13-cluster0_DE_genes_table-mt-25.csv")
# write.table(sig_res,"2024-9-13-cluster1_DE_genes_table-mt-25.csv")

########################################################################################################
all <- data1
all <- subset(all, subset = embryo %in% c("1-1","1-3","2-1","3-2","3-3","3-6","3-7","4-1","4-27","5-5.2","5-5.4","5-6.1","5-7.1","5-7.4","5-7.5"))
all$group_id <- all$label
all$sample_id <- all$embryo

#counts <- all@assays$RNA@counts
counts <- GetAssayData(all,slot = "counts")
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

colnames(pb[[1]] ) <- levels(sce$sample_id)#cluster0
colnames(pb[[2]] ) <-c("1-1","1-3","2-1","3-2","3-3","3-6","3-7","4-1","4-27","5-5.2","5-5.4","5-7.1","5-7.4","5-7.5")#
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

metadata
a <- c(rep("exp1 ",each = 2),rep("exp2 ",each = 1),rep("exp3 ",each = 4),rep("exp5 ",each = 1),rep("exp4 ",each= 1),rep("exp6 ",each = 2)
       ,rep("exp7 ",each = 1),rep("exp8 ",each = 3),rep("exp1 ",each = 2),rep("exp2 ",each = 1),rep("exp3 ",each = 4),rep("exp5 ",each = 1),
       rep("exp4 ",each = 1),rep("exp6 ",each = 2),rep("exp8 ",each = 3),rep("NA",each = 17))
metadata$experiment <- a

metadata$cluster_id <- as.factor(metadata$cluster_id)
clusters <- levels(metadata$cluster_id)
clusters
###################################################################################
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)
counts <- pb[[clusters[1]]]

cluster_counts <- data.frame(counts[,which(colnames(counts ) %in% rownames(cluster_metadata))])

colnames(cluster_counts) <- rownames(cluster_metadata)


dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ experiment + group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
# DESeq2::plotPCA(rld, intgroup = c("donor"))
DESeq2::plotPCA(rld, intgroup = c("group_id","experiment"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id","experiment"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
# plotDispEsts(dds)
contrast <- c("group_id","young","old")
resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 coef = "group_id_young_vs_old",
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

# write.table(sig_res,"2024-9-13-cluster1_euploid.csv")