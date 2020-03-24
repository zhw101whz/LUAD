##############################
# 2020_03_21
# zhw101whz
# 安装包
##############################

rm(list = ls())
options(stringsAsFactors = F)

### 设置清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

### 安装BioManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 

### 安装包
# CRAN的包
if(!require("data.table"))
  install.packages("data.table")
if(!require("stringr"))
  install.packages("stringr")
if(!require("pheatmap"))
  install.packages("pheatmap")
if(!require("ggplot2"))
  install.packages("ggplot2")
if(!require("ggfortify"))
  install.packages("ggfortify")
if(!require("VennDiagram"))
  install.packages("VennDiagram")

# Bioconductor的包
#if(!require("biomaRt"))
#  BiocManager::install("biomaRt", ask = F, update = F)
if(!require("DESeq2"))
  BiocManager::install("DESeq2", ask = F, update = F)
if(!require("edgeR"))
  BiocManager::install("edgeR", ask = F, update = F)
if(!require("org.Hs.eg.db"))
  BiocManager::install("org.Hs.eg.db", ask = F, update = F)
if(!require("clusterProfiler"))
  BiocManager::install("clusterProfiler", ask = F, update = F)

###############################################################
# Thank Dr.Jianming Zeng(University of Macau),
# and all the members of his bioinformatics team, biotrainee,
# for generously sharing their experience and codes !
###############################################################


