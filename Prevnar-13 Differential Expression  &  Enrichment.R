#Differential Expression & Enrichment

#performing volcano plot to visualize the DEGs"
dge_tophits_significant = dge_tophits[dge_tophits$P.Value < 0.05,]
upreg_genes = dge_tophits_significant[dge_tophits_significant$logFC > 0,]
downreg_genes = dge_tophits_significant[dge_tophits_significant$logFC < 0,]
upreg_genes = upreg_genes[order(upreg_genes$logFC, decreasing = TRUE),][1:5,]
downreg_genes = downreg_genes[order(downreg_genes$logFC, decreasing = FALSE),][1:5,]
upreg_genes$trend = rep("Up-regulated")
downreg_genes$trend = rep("Down-regulated")
dge_info = rbind(upreg_genes, downreg_genes) 

library(ggrepel)
volcano_plot = volcano_plot + 
  +     geom_label_repel(
    +         data = dge_info, 
    +         mapping = aes(x = logFC, y = -log10(P.Value), label = gene_symbol), 
    +         color = "black",
    +         size = 8
     )
ggsave(plot = volcano_plot,
       +        filename = "volcano_plot.jpeg",
       +        dpi = 300,
       +        height = 10,
       +        width = 12,
       +        device = "jpeg",
       +        bg = "white")




#building a heatmap of differentially expressed genes (DEGs)

#Select upregulated genes
#Filter DE results: keeps only genes with positive log fold change (upregulated) and p-value < 0.05.
upreg_geneinfo = dge_tophits[dge_tophits$logFC > 0 & dge_tophits$P.Value < 0.05,] 
upreg_geneinfo = upreg_geneinfo[order(upreg_geneinfo$logFC, decreasing = TRUE),]
upreg_geneinfo = upreg_geneinfo$gene_Symbol[1:5]

#Select downregulated genes
#Filter DE results: keeps only genes with negative log fold change (downregulated) and p-value < 0.05.
downreg_geneinfo = dge_tophits[dge_tophits$logFC < 0 & dge_tophits$P.Value < 0.05,]
downreg_geneinfo = downreg_geneinfo[order(downreg_geneinfo$logFC, decreasing = FALSE),]
downreg_geneinfo = downreg_geneinfo$gene_Symbol[1:5]


#From the whole expression matrix, selects only the rows (genes) corresponding to the top 5 upregulated and top 5 downregulated genes.
clust_cols = hclust(as.dist(1 - cor(diff_exp_genes_subset, method = "spearman")), method = "complete")
clust_rows = hclust(as.dist(1 - cor(t(diff_exp_genes_subset), method = "pearson")), method = "complete")

clust_cols = hclust(as.dist(1 - cor(diff_exp_genes_subset, method = "spearman")), method = "complete")
clust_rows = hclust(as.dist(1 - cor(t(diff_exp_genes_subset), method = "pearson")), method = "complete")

#Assign gene modules
module_assign = cutree(clust_rows, k = 2)
module_colour = rainbow(n = length(unique(module_assign)), start = 0.1, end = 0.9)   
module_colour = module_colour[as.vector(module_assign)]
col_pal = colorRampPalette(colors = c("red", "white", "green"))(10)  #Define heatmap colors

#Annotate rows and samples
my_gene_col = data.frame(cluster = ifelse(test = module_assign == 1, yes = "upregulated", no = "downregulated")) 
my_sample_col = metainfo[,c("ID", "T0", "sample_info")]
my_sample_col = subset(my_sample_col, select = c(T0))

#Draw heatmap
heatmap = pheatmap(
    mat = diff_exp_genes_subset, 
    legend = TRUE,
    RowSideColors = module_colour, 
    col = col_pal, 
    scale = "row", 
    cluster_cols = FALSE, 
    cluster_rows = TRUE,
    labRow = NA, 
    density.info = "none", 
    trace = "none", 
    cexRow = 1, 
    keysize = 2, 
    annotation_row = my_gene_col, 
    annotation_col = my_sample_col, 
    labels_col = rep("", dim(my_sample_col)[1]), 
    cutree_rows = 2
)

expr_mtx = read.delim("GSE298577_Prevnar.D0.D7.RMA.txt")  #Reloads the expression matrix from a file.



#Sort the DEG table
dge_tophits = dge_tophits[order(dge_tophits$adj.P.Val),]
rownames(dge_tophits) = NULL
cerno_results = tmodCERNOtest(dge_tophits$gene_symbol, order.by = "pval") #CERNO test for enrichment
write.table(file="cerno_results.txt", x = cerno_results, sep="\t", col.names=TRUE, quote=FALSE) #Save results

if(dim(cerno_results)[1] >= 2) {
    cerno_results = cerno_results[order(cerno_results$adj.P.Val),]
    cerno_results = cerno_results[1:30,] 

#Identify significant genes
sig_genes = tmodDecideTests(
    g = dge_tophits$gene_symbol, 
    lfc = dge_tophits$logFC, 
    pval = dge_tophits$P.Value,
    lfc.thr = 0.01, 
    pval.thr = 0.05)
names(sig_genes) = "DGE_mods"
#Create enrichment panel plot
tmod_plots = ggPanelplot(res = list(DGE_mods = cerno_results), sgenes = sig_genes)
tmod_plots = tmod_plots + theme(...)
ggsave(plot = tmod_plots, filename = "cerno_resultes", device = "png", width = 14, height = 9, dpi = 300) #Save the plot
}
else{print("cerno test failed - less than 2 modules")}

#View and re-plot
tmod_plots <- ggPanelplot(res = list(DGE_mods = cerno_results), sgenes = sig_genes) + theme(...)
ggsave(plot = tmod_plots, filename = "cerno_results.png", width = 14, height = 9, dpi = 300)













