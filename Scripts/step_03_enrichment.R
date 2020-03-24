##############################
# 2020_03_21
# zhw101whz
# 富集分析
##############################

rm(list = ls())
options(stringsAsFactors = F)
# 设置科学计数法
options(scipen = 200)

library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")

### 加载数据
load("DESeq2_for_enrich.RData")
### id转换 SYMBOL -> ENTREZID 
symbol_id <- rownames(DESeq2_for_enrich)
ID <- bitr(symbol_id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 重复的不要
table(table(ID$ENTREZID))
table(table(ID$SYMBOL))
ENTREZID_single <- names(table(ID$ENTREZID) == 1)[table(ID$ENTREZID) == 1]
SYMBOL_single <- names(table(ID$SYMBOL) == 1)[table(ID$SYMBOL) == 1]
ID <- ID[ID$ENTREZID %in% ENTREZID_single, ]
ID <- ID[ID$SYMBOL %in% SYMBOL_single, ]
DESeq2_for_enrich <- DESeq2_for_enrich[rownames(DESeq2_for_enrich) %in% ID$SYMBOL, ]
# 改id
need_ID <- data.frame(SYMBOL = rownames(DESeq2_for_enrich))
need_ID <- merge(need_ID, ID, sort = FALSE)
rownames(DESeq2_for_enrich) <- need_ID$ENTREZID
### 分组处理
# 分改变、上调、下调
geneList_total <- DESeq2_for_enrich[, 1]
geneList_diff <- DESeq2_for_enrich[DESeq2_for_enrich$change != "no", 1]
geneList_up <- DESeq2_for_enrich[DESeq2_for_enrich$change == "up", 1]
geneList_down <- DESeq2_for_enrich[DESeq2_for_enrich$change == "down", 1]
# 给名字
names(geneList_total) <- rownames(DESeq2_for_enrich)
names(geneList_diff) <- rownames(DESeq2_for_enrich[DESeq2_for_enrich$change != "no", ])
names(geneList_up) <- rownames(DESeq2_for_enrich[DESeq2_for_enrich$change == "up", ])
names(geneList_down) <- rownames(DESeq2_for_enrich[DESeq2_for_enrich$change == "down", ])
# 排序
geneList_total <- sort(geneList_total, decreasing = TRUE)
geneList_diff <- sort(geneList_diff, decreasing = TRUE)
geneList_up <- sort(geneList_up, decreasing = TRUE)
geneList_down <- sort(geneList_down, decreasing = TRUE)
# 取基因名
gene_total <- names(geneList_total)
gene_diff <- names(geneList_diff)
gene_up <- names(geneList_up)
gene_down <- names(geneList_down)

### GO
# diff
if(T){
  go_CC_diff <- enrichGO(gene          = gene_diff,
                         universe      = gene_total,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_CC_diff, showCategory = 10)
  dotplot(go_CC_diff, showCategory = 10)
  ggsave("go_CC_diff_dotplot.png", width = 10, height = 6)
  cnetplot(go_CC_diff, foldChange = geneList_diff)
  heatplot(go_CC_diff, foldChange = geneList_diff)
  emapplot(go_CC_diff)
  ggsave("go_CC_diff_emapplot.png", width = 10, height = 10)
  
  go_BP_diff <- enrichGO(gene          = gene_diff,
                         universe      = gene_total,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_BP_diff, showCategory = 10)
  dotplot(go_BP_diff, showCategory = 10)
  ggsave("go_BP_diff_dotplot.png", width = 10, height = 6)
  cnetplot(go_BP_diff, foldChange = geneList_diff)
  heatplot(go_BP_diff, foldChange = geneList_diff)
  emapplot(go_BP_diff)
  ggsave("go_BP_diff_emapplot.png", width = 10, height = 10)
  
  go_MF_diff <- enrichGO(gene          = gene_diff,
                         universe      = gene_total,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_MF_diff, showCategory = 10)
  dotplot(go_MF_diff, showCategory = 10)
  ggsave("go_MF_diff_dotplot.png", width = 10, height = 6)
  cnetplot(go_MF_diff, foldChange = geneList_diff)
  heatplot(go_MF_diff, foldChange = geneList_diff)
  emapplot(go_MF_diff)
  ggsave("go_MF_diff_emapplot.png", width = 10, height = 10)
}
# up
if(T){
  go_CC_up <- enrichGO(gene          = gene_up,
                       universe      = gene_total,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_CC_up, showCategory = 10)
  dotplot(go_CC_up, showCategory = 10)
  ggsave("go_CC_up_dotplot.png", width = 10, height = 6)
  cnetplot(go_CC_up, foldChange = geneList_up)
  heatplot(go_CC_up, foldChange = geneList_up)
  emapplot(go_CC_up)
  ggsave("go_CC_up_emapplot.png", width = 10, height = 10)
  
  go_BP_up <- enrichGO(gene          = gene_up,
                       universe      = gene_total,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_BP_up, showCategory = 10)
  dotplot(go_BP_up, showCategory = 10)
  ggsave("go_BP_up_dotplot.png", width = 10, height = 6)
  cnetplot(go_BP_up, foldChange = geneList_up)
  heatplot(go_BP_up, foldChange = geneList_up)
  emapplot(go_BP_up)
  ggsave("go_BP_up_emapplot.png", width = 10, height = 10)
  
  go_MF_up <- enrichGO(gene          = gene_up,
                       universe      = gene_total,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_MF_up, showCategory = 10)
  dotplot(go_MF_up, showCategory = 10)
  ggsave("go_MF_up_dotplot.png", width = 10, height = 6)
  cnetplot(go_MF_up, foldChange = geneList_up)
  heatplot(go_MF_up, foldChange = geneList_up)
  emapplot(go_MF_up)
  ggsave("go_MF_up_emapplot.png", width = 10, height = 10)
}
# down
if(T){
  go_CC_down <- enrichGO(gene          = gene_down,
                         universe      = gene_total,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_CC_down, showCategory = 10)
  dotplot(go_CC_down, showCategory = 10)
  ggsave("go_CC_down_dotplot.png", width = 10, height = 6)
  cnetplot(go_CC_down, foldChange = geneList_down)
  heatplot(go_CC_down, foldChange = geneList_down)
  emapplot(go_CC_down)
  ggsave("go_CC_down_emapplot.png", width = 10, height = 10)
  
  go_BP_down <- enrichGO(gene          = gene_down,
                         universe      = gene_total,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_BP_down, showCategory = 10)
  dotplot(go_BP_down, showCategory = 10)
  ggsave("go_BP_down_dotplot.png", width = 10, height = 6)
  cnetplot(go_BP_down, foldChange = geneList_down)
  heatplot(go_BP_down, foldChange = geneList_down)
  emapplot(go_BP_down)
  ggsave("go_BP_down_emapplot.png", width = 10, height = 10)
  
  go_MF_down <- enrichGO(gene          = gene_down,
                         universe      = gene_total,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  # 可视化，按需保存吧
  barplot(go_MF_down, showCategory = 10)
  dotplot(go_MF_down, showCategory = 10)
  ggsave("go_MF_down_dotplot.png", width = 10, height = 6)
  cnetplot(go_MF_down, foldChange = geneList_down)
  heatplot(go_MF_down, foldChange = geneList_down)
  emapplot(go_MF_down)
  ggsave("go_MF_down_emapplot.png", width = 10, height = 10)
}

### KEGG 这好像还得联网
# diff
if(T){
  kegg_diff <- enrichKEGG(gene         = gene_diff,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05,
                          universe = gene_total)
  # 可视化，按需保存吧
  barplot(kegg_diff, showCategory = 10)
  dotplot(kegg_diff, showCategory = 10)
  ggsave("kegg_diff_dotplot.png", width = 10, height = 6)
  cnetplot(kegg_diff, foldChange = geneList_diff)
  heatplot(kegg_diff, foldChange = geneList_diff)
  emapplot(kegg_diff)
  ggsave("kegg_diff_emapplot.png", width = 10, height = 10)
}
# up
if(T){
  kegg_up <- enrichKEGG(gene         = gene_up,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05,
                        universe = gene_total)
  # 可视化，按需保存吧
  barplot(kegg_up, showCategory = 10)
  dotplot(kegg_up, showCategory = 10)
  ggsave("kegg_up_dotplot.png", width = 10, height = 6)
  cnetplot(kegg_up, foldChange = geneList_up)
  heatplot(kegg_up, foldChange = geneList_up)
  emapplot(kegg_up)
  ggsave("kegg_up_emapplot.png", width = 10, height = 10)
}
# down
if(T){
  kegg_down <- enrichKEGG(gene         = gene_down,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05,
                          universe = gene_total)
  # 可视化，按需保存吧
  barplot(kegg_down, showCategory = 10)
  dotplot(kegg_down, showCategory = 10)
  ggsave("kegg_down_dotplot.png", width = 10, height = 6)
  cnetplot(kegg_down, foldChange = geneList_down)
  heatplot(kegg_down, foldChange = geneList_down)
  emapplot(kegg_down)
  ggsave("kegg_down_emapplot.png", width = 10, height = 10)
}

###############################################################
# Thank Dr.Jianming Zeng(University of Macau),
# and all the members of his bioinformatics team, biotrainee,
# for generously sharing their experience and codes !
###############################################################