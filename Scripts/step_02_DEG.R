##############################
# 2020_03_20
# zhw101whz
# 差异表达分析
##############################

rm(list = ls())
options(stringsAsFactors = F)

### 加载包
library(DESeq2)
library(edgeR)
library(limma)
library(VennDiagram)

### 加载自定义绘图函数
source("draw_h_p_v.R")

### 载入数据
load("exp_df.RData")

### 根据TCGA样本命名分组
group_list <- factor(ifelse(as.numeric(substr(colnames(exp_df), 14, 15)) < 10, "tumor", "normal"))

### DESeq2
if(T){
  # 构建分组
  colData <- data.frame(row.names = colnames(exp_df), group_list = group_list)
  # 整合数据
  dds = DESeqDataSetFromMatrix(countData = exp_df,
                               colData = colData,
                               design = ~ group_list)
  # 建模，这个很慢，保存一下
  if(!file.exists("TCGA_LUAD_mRNA_DESeq2_dds.RData")){
    dds = DESeq(dds)
    save(dds, file = "TCGA_LUAD_mRNA_DESeq2_dds.RData")
  }
  load(file = "TCGA_LUAD_mRNA_DESeq2_dds.RData")
  # 提取结果
  resultsNames(dds)
  res = results(dds, name = "group_list_tumor_vs_normal")
  # 只要logFC和pvalue
  res_ordered = res[order(res$padj), ]
  head(res_ordered)
  DEG = as.data.frame(res_ordered)
  DESeq2_DEG = na.omit(DEG)
  DESeq2_need_DEG = DESeq2_DEG[, c(2, 6)]
  colnames(DESeq2_need_DEG) = c("log2FoldChange", "pvalue")
  DESeq2_for_enrich <- DESeq2_need_DEG
  DESeq2_for_enrich$change <- ifelse((abs(DESeq2_for_enrich$log2FoldChange) <= 1) | (DESeq2_for_enrich$pvalue >= 0.05), "no",
                                     ifelse(DESeq2_for_enrich$log2FoldChange > 0, "up", "down"))
  
  ### 绘图
  # 改完函数重画时记得source
  draw_h_p_v(exp_df, DESeq2_need_DEG, "DESeq2", group_list, 1)
}

### edgeR
if(T){
  y <- DGEList(counts = exp_df, group = group_list)
  y <- calcNormFactors(y)
  design <- model.matrix(~group_list)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  topTags(qlf)
  edgeR_DEG <- topTags(qlf, n = nrow(y))
  edgeR_DEG <- as.data.frame(edgeR_DEG)
  edgeR_need_DEG <- edgeR_DEG[,c(1,5)]
  colnames(edgeR_need_DEG) <- c('log2FoldChange','pvalue') 
  edgeR_for_enrich <- edgeR_need_DEG
  edgeR_for_enrich$change <- ifelse((abs(edgeR_for_enrich$log2FoldChange) <= 1) | (edgeR_for_enrich$pvalue >= 0.05), "no",
                                    ifelse(edgeR_for_enrich$log2FoldChange > 0, "up", "down"))
  draw_h_p_v(exp_df, edgeR_need_DEG, "edgeR", group_list, 1)
}

### limma + voom
if(T){
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exp_df)
  
  d <- DGEList(counts = exp_df)
  d <- calcNormFactors(d)
  logCPM <- cpm(d, log=TRUE, prior.count = 3)
  
  v = voom(d, design, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend=TRUE)
  
  tempOutput <- topTable(fit, coef = 2, n = Inf)
  limma_voom_DEG <- na.omit(tempOutput)
  limma_need_DEG <- limma_voom_DEG[,c(1,4)]
  colnames(limma_need_DEG) = c('log2FoldChange','pvalue') 
  limma_for_enrich <- limma_need_DEG
  limma_for_enrich$change <- ifelse((abs(limma_for_enrich$log2FoldChange) <= 1) | (limma_for_enrich$pvalue >= 0.05), "no",
                                    ifelse(limma_for_enrich$log2FoldChange > 0, "up", "down"))
  draw_h_p_v(exp_df, limma_need_DEG, "limma", group_list, 1)
}

### 三种方法的相关性
if(T){
  mi = unique(c(rownames(DESeq2_need_DEG), rownames(edgeR_need_DEG), rownames(limma_need_DEG)))
  lf = data.frame(lf1 = DESeq2_need_DEG[mi, 1],
                  lf2 = edgeR_need_DEG[mi, 1],
                  lf3 = limma_need_DEG[mi, 1])
  cor(na.omit(lf))
  # 可以看到采取不同R包，会有不同的归一化算法，这样算到的logFC会稍微有差异。
}

### DEG韦恩图
if(T){
  t1 = rownames(DESeq2_for_enrich[DESeq2_for_enrich$change != "no", ])
  t2 = rownames(edgeR_for_enrich[edgeR_for_enrich$change != "no", ])
  t3 = rownames(limma_for_enrich[limma_for_enrich$change != "no", ])
  venn.diagram(list(DESeq2 = t1, edgeR = t2, "limma + voom" = t3), "DEG3_Venn.png",
               col = "transparent", alpha = c(0.5, 0.5, 0.5),
               fill = c("red", "yellow", "green"), main = "Three DEG_method_result Venn")
}

### 用DESeq2的数据做富集分析
save(DESeq2_for_enrich, file = "DESeq2_for_enrich.RData")

###############################################################
# Thank Dr.Jianming Zeng(University of Macau),
# and all the members of his bioinformatics team, biotrainee,
# for generously sharing their experience and codes !
###############################################################
