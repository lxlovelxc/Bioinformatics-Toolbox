# THIS IS THE CODE TO PLOT THE PLATE HEATMAP BASED ON UMI
# THIS IS USED FOR HUMAN EMBRYO SCS 20221207 
# XI LUO

library(ggplot2)
library(reshape2)
library(Seurat)
library(dplyr)
library(pheatmap)
library(ramify)
library(readxl)
# Load ggplate package
library(ggplate)

setwd("/home/xi/Desktop")
white_list_check <- read_excel("white_list_check.xlsx",col_names = FALSE)

embryo_data_list=c("1-1","1-2","1-3","1-4","1-5_5-1","2-1","2-2","3-1","4-1","4-2","4-3","5-2","6")
data_list = list()
for (n in 1:13){
  data_list[[n]]=Read10X(data.dir = paste0("/home/xi/Desktop/human_embryo_scs_20221207_data/",embryo_data_list[n]))
  print(n)
}

ordered_wl = flatten(as.matrix(white_list_check), across = c("rows", "columns"))


heatmap_plot <- function(plate){
  a <- Matrix::colSums(data_list[[plate]])
  a <- data.frame(a)
  a$barcode <- rownames(a)
  a <- a[ordered_wl,]
  AtoH <- c('A','B','C','D','E','F','G','H')
  a$Y <- rep(AtoH,each = 12)
  a$X <- rep((c(1:12)) , time = 8)
  
  ggplot(data = a) + geom_count(mapping = aes(Y, X, color = log10(a)), size = 8, shape = 16) + 
    scale_color_gradient(low = "white", high = "blue") + theme_bw() + theme_classic() +
    labs(x = NULL , y = paste('plate plot',embryo_data_list[plate]))+
    scale_y_continuous(limits = c(1,12),breaks = c(12,11,10,9,8,7,6,5,4,3,2,1))+
    theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1
    )) +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(axis.text.y = element_text(angle = 90, hjust = 1))+ coord_fixed(ratio=1.5)
  
}


for (i in 1:2){
  print(heatmap_plot(i))
}
