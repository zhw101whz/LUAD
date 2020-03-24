##############################
# 2020_03_23
# zhw101whz
# 准备数据
##############################

rm(list = ls())
options(stringsAsFactors = F)

### 加载包
library(data.table)
library(stringr)
library(biomaRt)

### 读取数据
exp_df <- fread(file = "TCGA-LUAD.htseq_counts.tsv", header = T)
# 简单看一下，数据是log2(count+1)
exp_df[1:4, 1:4]

### 简单处理
# 去Ensembl_ID版本号
exp_df$Ensembl_ID <- str_split(exp_df$Ensembl_ID, "\\.", simplify = T)[, 1]
# 加行名
class(exp_df)
exp_df <- as.data.frame(exp_df)
rownames(exp_df) <- exp_df$Ensembl_ID
exp_df <- exp_df[, -1]
# 图方便，只要配对数据(这一段写的太乱)
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
# 过滤低表达
index <- apply(exp_df[, -1], 1, function(x){
  sum(x < log2(10 + 1)) < length(x)/2
})
exp_df <- exp_df[index, ]

### biomart的id转换，先不用了
if(F){
  ### id转换，只要coding
  my_ensembl_gene_id <- rownames(exp_df)
  # 显示能连接的数据库
  listMarts()
  # 选择数据库
  ensembl <- useMart("ensembl")
  # 查看数据集
  exp_dfasets <- listDatasets(ensembl)
  # 选人类的
  ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
  # 查看注释类型
  filters = listFilters(ensembl)
  # 查看可用属性
  attributes = listAttributes(ensembl)
  # 提取gene_biotype
  symbol_type <- getBM(attributes = c("ensembl_gene_id",
                                      "hgnc_symbol",
                                      "gene_biotype"), 
                       filters= "ensembl_gene_id",
                       values = my_ensembl_gene_id,
                       mart = ensembl)
  symbol <- symbol_type[symbol_type$gene_biotype == "protein_coding", ]
  table(table(symbol$hgnc_symbol))
  symbol <- symbol[symbol$hgnc_symbol != "", ]
  table(table(symbol$ensembl_gene_id))
  names(table(symbol$ensembl_gene_id))[table(symbol$ensembl_gene_id) == 2]
  symbol <- symbol[symbol$ensembl_gene_id != "ENSG00000187510", ]
  symbol <- symbol[symbol$ensembl_gene_id != "ENSG00000276085", ]
  
  exp_df <- exp_df[symbol$ensembl_gene_id, ]
  ensembl_id <- data.frame(ensembl_gene_id = rownames(exp_df))
  order_symbol <- merge(x = ensembl_id, y = symbol, by = "ensembl_gene_id", all.x = T, sort = F)
  
  rownames(exp_df) <- order_symbol$hgnc_symbol
}

### id转换
keytypes(org.Hs.eg.db)
ID <- bitr(rownames(exp_df), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
# 重复的不要
table(table(ID$ENSEMBL))
table(table(ID$SYMBOL))
ENSEMBL_single <- names(table(ID$ENSEMBL) == 1)[table(ID$ENSEMBL) == 1]
SYMBOL_single <- names(table(ID$SYMBOL) == 1)[table(ID$SYMBOL) == 1]
ID <- ID[ID$ENSEMBL %in% ENSEMBL_single, ]
ID <- ID[ID$SYMBOL %in% SYMBOL_single, ]
exp_df <- exp_df[rownames(exp_df) %in% ID$ENSEMBL, ]
need_id <- data.frame(ENSEMBL = rownames(exp_df))
need_id <- merge(need_id, ID, sort = F)
rownames(exp_df) <- need_id$SYMBOL

### 由log2(count+1)还原成count
exp_df <- round(2^exp_df - 1)

save(exp_df, file = "exp_df.RData")

###############################################################
# Thank Dr.Jianming Zeng(University of Macau),
# and all the members of his bioinformatics team, biotrainee,
# for generously sharing their experience and codes !
###############################################################


