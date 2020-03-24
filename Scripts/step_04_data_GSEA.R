##############################
# 2020_03_22
# zhw101whz
# 准备数据 for GSEA
##############################

rm(list = ls())
options(stringsAsFactors = F)

### 加载包
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(DESeq2)

### 读取数据
exp_df <- fread(file = "TCGA-LUAD.htseq_counts.tsv", header = T)
exp_df[1:3, 1:3]

### id转换
# 去Ensembl_ID版本号
exp_df$Ensembl_ID <- str_split(exp_df$Ensembl_ID, "\\.", simplify = T)[, 1]
exp_df[1:3, 1:3]
# Ensembl_ID -> HUGO_SYMBOL
keytypes(org.Hs.eg.db)
ID <- bitr(exp_df$Ensembl_ID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
# 重复的不要
table(table(ID$ENSEMBL))
table(table(ID$SYMBOL))
ENSEMBL_single <- names(table(ID$ENSEMBL) == 1)[table(ID$ENSEMBL) == 1]
SYMBOL_single <- names(table(ID$SYMBOL) == 1)[table(ID$SYMBOL) == 1]
ID <- ID[ID$ENSEMBL %in% ENSEMBL_single, ]
ID <- ID[ID$SYMBOL %in% SYMBOL_single, ]
exp_df <- exp_df[exp_df$Ensembl_ID %in% ID$ENSEMBL, ]
# 替换
exp_df <- merge(x = exp_df, y = ID, by.x = "Ensembl_ID", by.y = "ENSEMBL")
# 加行名
class(exp_df)
exp_df <- as.data.frame(exp_df)
rownames(exp_df) <- exp_df$SYMBOL
exp_df <- exp_df[, c(-1, -ncol(exp_df))]

### 图方便，只要配对数据(这段写的太乱)
barcode <- data.frame(sample_barcode = colnames(exp_df),
                      patient_barcode = substring(colnames(exp_df), 1, 12),
                      type = substring(colnames(exp_df), 14, 15))
table(table(barcode$patient_barcode))
one <- names(table(barcode$patient_barcode))[table(barcode$patient_barcode) == 1]
barcode$one <- barcode$patient_barcode %in% one
need <- barcode[barcode$one == "FALSE", -ncol(barcode)]
need_n <- need[need$type > 10, ]
table(table(need_n$patient_barcode))
need_t <- need[need$type < 10, ]
need_t <- need_t[substring(need_t$sample_barcode, 16, 16) == "A", ]
need <- merge(need_n, need_t, by.x = "patient_barcode", by.y = "patient_barcode", all.x = T)
need_simple_barcode <- c(need$sample_barcode.x, need$sample_barcode.y)
exp_df <- exp_df[, need_simple_barcode]

### 用DESeq2标准化
# 先还原count值
exp_df <- round(2^exp_df - 1)
# 根据TCGA样本命名分组
group_list <- factor(ifelse(as.numeric(substr(colnames(exp_df), 14, 15)) < 10, "tumor", "normal"))
# 分组信息
colData <- data.frame(row.names = colnames(exp_df), group_list = group_list)
# 整合数据
dds <- DESeqDataSetFromMatrix(countData = exp_df,
                              colData = colData,
                              design = ~ group_list)
# 建模
dds <- DESeq(dds)
normalized_data <- counts(dds, normalized = TRUE)
boxplot(log10(normalized_data[, 1:ncol(normalized_data)] + 1))
# 看下原来的,齐多了，算法决定了均值不一样
boxplot(log10(exp_df + 1))
# 数据按要求写出去，不能有-
normalized_data <- as.data.frame(normalized_data)
colnames(normalized_data) <- gsub("-", "_", colnames(normalized_data))
GSEA_data <- data.frame(NAME = rownames(normalized_data),
                        DESCRIPTION = rep(NA, nrow(normalized_data)))
GSEA_data <- cbind(GSEA_data, normalized_data)
GSEA_data$NAME <- gsub("-", "_", GSEA_data$NAME)
write.table(GSEA_data, file = "GSEA_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 整理phenotype数据，还得自己加两行
phenotype <- as.data.frame(t(as.data.frame(group_list)))
phenotype[1, ] <- apply(phenotype, 2, function(x){
  ifelse(x == "normal", 0, 1)
})
phenotype[1, ] <- as.numeric(phenotype[1, ])
write.table(phenotype, file = "GSEA_pheno.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

###############################################################
# Thank Dr.Jianming Zeng(University of Macau),
# and all the members of his bioinformatics team, biotrainee,
# for generously sharing their experience and codes !
###############################################################