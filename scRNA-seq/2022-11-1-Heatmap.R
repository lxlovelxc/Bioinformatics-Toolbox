# THIS IS THE CODE TO PLOT THE PLATE HEATMAP BASED ON UMI
# THIS IS USED FOR HUMAN EMBRYO SCS 20221101 
# XI LUO

library(ggplot2)
library(reshape2)
library(Seurat)
library(dplyr)
library(pheatmap)
library(ramify)
library(readxl)
library(tidyverse)
library(aplot)


embryo_data_list=c("1-1","1-3","2-2","3-1","3-2","6-1","6-2","6-3")
data_list = list()
for (n in 1:8){
  data_list[[n]]=Read10X(data.dir = paste0("C:/Users/luoxi/Desktop/human_embryo_20221101/",embryo_data_list[n]))
  print(n)
}


white_list_check <- read_excel("white_list_check.xlsx",col_names = FALSE)


ordered_wl = flatten(as.matrix(white_list_check), across = c("rows", "columns"))


for(n in 1:8){
  col_sum  = Matrix::colSums(data_list[[n]])
  col_sum
  matrix_col_sum = as.matrix(col_sum[ordered_wl])
  matrix_col_sum = matrix(matrix_col_sum, nrow = 8, ncol = 12,byrow = TRUE)
  pheatmap(log10(matrix_col_sum+1),main = paste0("plate",n),cluster_row = FALSE,cluster_cols = FALSE,cellwidth = 15, cellheight = 10)
  print(n)
}

(col_sum[ordered_wl])

col_sum = Matrix::colSums(data_list[[1]])
col_sum
matrix_col_sum = as.matrix(col_sum[ordered_wl])
matrix_col_sum = matrix(matrix_col_sum, nrow = 8, ncol = 12,byrow = TRUE)
pheatmap(log10(matrix_col_sum+1),main = paste0("plate",1),cluster_row = FALSE,cluster_cols = FALSE,cellwidth = 15, cellheight = 10)
matrix_col_sum
col_sum

ordered_wl
matrix_col_sum


