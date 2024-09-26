table <- read.csv("/home/xi/Desktop/table1.csv")
table <- as.data.frame(table)
mean(table$nCount_peaks)
View(table)
mean(table$frip) #55

sample <- read.csv("/home/xi/Desktop/0725-SMART-ATAC/ATAC/sample_info.csv")
sample <- as.data.frame(sample)
View(sample)
# mean(table$frac_open,na.rm = TRUE)

# 在计算前先排除frac_open列中的0和NA项
filtered_data <- sample[sample$frac_open != 0 & !is.na(sample$frac_open), ]
mean_frac_open <- mean(filtered_data$frac_open)
mean_frac_open # 6.3 	  1 - 20

# 在计算前先排除mapping rate列中的0和NA项
filtered_data_mr <- sample[sample$mapping_rate != 0 & !is.na(sample$mapping_rate), ]
mean_mr <- mean(filtered_data_mr$mapping_rate)
mean_mr  # 94  	70 - 99

# 在计算前先排除mt content
filtered_data_mt <- sample[sample$mt_content != 0 & !is.na(sample$mt_content), ]
mean_mt <- mean(filtered_data_mt$mt_content)
mean_mt  # 4.09  0.1 - 90

#uniq fragment
filtered_data_uniq <- sample[sample$uniq_nuc_frags != 0 & !is.na(sample$uniq_nuc_frags), ]
mean_uniq<- mean(filtered_data_uniq$uniq_nuc_frags)
mean_uniq  # 30519  10,000 - 1,000,000

#dup level
filtered_data_dup <- sample[sample$dup_level != 0 & !is.na(sample$dup_level), ]
mean_dup<- mean(filtered_data_dup$dup_level)
mean_dup # 0.458  0.4 - 0.9

#frip
filtered_data_frip <- sample[sample$frip != 0 & !is.na(sample$frip), ]
mean_frip<- mean(filtered_data_frip$frip)
mean_frip  # 7.16  20 - 80

#sequence depth
filtered_data_seq <- sample[sample$sequencing_depth != 0 & !is.na(sample$sequencing_depth), ]
mean_seq<- mean(filtered_data_seq$sequencing_depth)
mean_seq   # 74538  10,000 - 1,000,000

#library_size
filtered_data_lib <- sample[sample$library_size != 0 & !is.na(sample$library_size), ]
mean_lib<- mean(filtered_data_lib$library_size)
mean_lib #73553  10,000 - 1,000,000
library(ggplot2)
# 加载所需包
library(ggplot2)

# 打开PDF设备
pdf("ATAC-quality.pdf")

# 图1：Library Size
p1 <- ggplot(sample, aes(x = cell, y = log10(library_size), fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = ifelse(log10(sample$library_size) >= 5 & log10(sample$library_size) <= 7, "red", "black"))) +
  labs(title = "Library Size for Each Sample", x = "Sample", y = "log10_Library Size") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 7, color = "blue", linetype = "dashed")
print(p1)

# 图2：Mapping Rate
p2 <- ggplot(sample, aes(x = cell, y = mapping_rate, fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = ifelse(sample$mapping_rate >= 70 & sample$mapping_rate <= 99, "red", "black"))) +
  labs(title = "Mapping rate for Each Sample", x = "Sample", y = "Mapping Rate") +
  geom_hline(yintercept = 70, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 99, color = "blue", linetype = "dashed")
print(p2)

# 图3：Duplication Level
p3 <- ggplot(sample, aes(x = cell, y = dup_level, fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = ifelse(sample$dup_level >= 0.4 & sample$dup_level <= 0.9, "red", "black"))) +
  labs(title = "Duplication level for Each Sample", x = "Sample", y = "Duplication level") +
  geom_hline(yintercept = 0.4, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.9, color = "blue", linetype = "dashed")
print(p3)

# 图4：FRIP
p4 <- ggplot(sample, aes(x = cell, y = frip, fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "frip for Each Sample", x = "Sample", y = "frip %") +
  geom_hline(yintercept = 20, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 80, color = "blue", linetype = "dashed")
print(p4)

# 图5：Uniq Nuc Frags
p5 <- ggplot(sample, aes(x = cell, y = log10(uniq_nuc_frags), fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = ifelse(log10(sample$uniq_nuc_frags) >= 5 & log10(sample$uniq_nuc_frags) <= 6, "red", "black"))) +
  labs(title = "uniq_nuc_frags for Each Sample", x = "Sample", y = "log10_uniq_nuc_frags") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 6, color = "blue", linetype = "dashed")
print(p5)

# 图6：Frac Open
p6 <- ggplot(sample, aes(x = cell, y = frac_open, fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = ifelse(sample$frac_open >= 1 & sample$frac_open <= 20, "red", "black"))) +
  labs(title = "frac_open for Each Sample", x = "Sample", y = "frac_open") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 20, color = "blue", linetype = "dashed")
print(p6)

# 关闭PDF设备
dev.off()
