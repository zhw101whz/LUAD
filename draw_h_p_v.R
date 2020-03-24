### 绘制DEG的热图、PCA图、火山图
### 需要DEG的log2FoldChange和pvalue

draw_h_p_v = function(exp_df, need_DEG, DEG_method, group_list, logFC_cutoff){
  
  ### 热图
  library(pheatmap)
  choose_gene = head(rownames(need_DEG), 100) ## 50 maybe better
  choose_matrix = exp_df[choose_gene, ]
  choose_matrix = t(scale(t(log2(choose_matrix+1)))) 
  choose_matrix[choose_matrix > 2.5] = 2.5
  choose_matrix[choose_matrix < -2.5] = -2.5
  annotation_col = data.frame(Groups = group_list)
  rownames(annotation_col) = colnames(exp_df)
  
  # 分类排序
  #callback = function(hc, mat){
  #  sv = svd(t(mat))$v[,1]
  #  dend = reorder(as.dendrogram(hc), wts = sv)
  #  as.hclust(dend)
  #}
  
  # 设定annotation颜色
  ann_colors = list(
    Groups = c(normal = "#88F93D", tumor = "#56C1CD")
  )
  
  # 绘图
  pheatmap(choose_matrix, show_colnames = F, annotation_col = annotation_col,
           fontsize = 10, fontsize_row = 4, main = "Top 100 of differentially expressed mRNA",
           filename = paste0(DEG_method,'_DEG_top100_heatmap.png'),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           annotation_colors = ann_colors
           #clustering_callback = callback,
           #cutree_col = 3, cutree_row = 3
  )
  
  ### PCA图
  library(ggfortify)
  df = as.data.frame(t(choose_matrix))
  df$group = group_list
  #png(paste0(DEG_method, "_DEG_top100_PCA.png"), res = 120)
  p = autoplot(prcomp(df[, 1:(ncol(df)-1)]), data = df, colour = 'group') + theme_bw()
  #print(p)
  ggsave(p, filename = paste0(DEG_method, "_DEG_top100_PCA.png"))
  
  ### 火山图
  # 设定
  need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff, "UP", "DOWN"), "NOT")
  )
  this_tile = paste0("Cutoff for logFC is ", round(logFC_cutoff, 3),
                     "\nThe number of up gene is ", nrow(need_DEG[need_DEG$change == "UP", ]),
                     "\nThe number of down gene is ", nrow(need_DEG[need_DEG$change == "DOWN", ]))
  
  # 绘图
  g = ggplot(data = need_DEG, aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
             geom_point(alpha = 0.4, size = 1.75) + 
             theme_set(theme_set(theme_bw(base_size = 20))) +
             xlab("log2 fold change") + ylab("-log10 p-value") +
             ggtitle(this_tile) +
             xlim(-10, 10) +
             theme(plot.title = element_text(size = 15,hjust = 0.5)) +
             scale_colour_manual(values = c("blue", "black", "red"))
  ggsave(g, filename = paste0(DEG_method, "_volcano.png"))
}
