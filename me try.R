expr_mtx = read.delim("GSE298577_Prevnar.D0.D7.RMA.txt")
gene_annot = read.delim("GPL13158-5065.txt", skip = 11)
gene_annot = gene_annot[,c(1,3)]
gene_annot = gene_annot[!duplicated(gene_annot$Gene.Symbol),]
gene_annot = gene_annot[!is.na(gene_annot$Gene.Symbol),]
gene_annot$Gene.Symbol = unlist(lapply(1:nrow(gene_annot), function(x){strsplit(gene_annot$Gene.Symbol[x], "///")[[1]][1]}))
expr_mtx = expr_mtx[expr_mtx$Pset.ID %in% gene_annot$ID, ]
expr_mtx = merge(expr_mtx, gene_annot, by.x = "Pset.ID", by.y = "ID")
expr_mtx = expr_mtx[!duplicated(expr_mtx$Gene.Symbol),]
expr_mtx = expr_mtx[!is.na(expr_mtx$Gene.Symbol),]
rownames(expr_mtx) = expr_mtx$Gene.Symbol
expr_mtx = expr_mtx[,-c(1,84)]
rm(gene_annot)


#meta data
metainfo = data.frame(ID = colnames(expr_mtx))
metainfo$sample_info = sub("Subj\\.(PNM[A-D][0-9]+_Day[0-9]+).*", "\\1", metainfo$ID)
metainfo$day = unlist(lapply(1:nrow(metainfo), function(x){strsplit(metainfo$sample_info[x], "_")[[1]][2]}))
metainfo = metainfo[order(metainfo$day),]
metainfo = metainfo[,-3]
rownames(metainfo)=NULL
metainfo$T0 = c(rep("1", 41), rep("0", 41))
metainfo$T7 = c(rep("0", 41), rep("1", 41))
rownames(metainfo) = metainfo$sample_info
colnames(expr_mtx) = sub("Subj\\.(PNM[A-D][0-9]+_Day[0-9]+).*", "\\1", colnames(expr_mtx))
expr_mtx = expr_mtx[,match(rownames(metainfo), colnames(expr_mtx))]
form = ~(T7-T0)
dge_dream_fit = dream(exprObj = expr_mtx, formula = form, data = metainfo, BPPARAM = param)
dge_dream_fit_ebayes = eBayes(dge_dream_fit)
ge_tophits = limma::topTable(dge_dream_fit_ebayes, number = 20000, sort.by = "logFC")
dge_tophits$gene_symbol = rownames(dge_tophits)

#volcanoplot

step4_makeVolcanoHeatplots = function(dge_tophits, metainfo, expr_mtx){
       
    +     print("Step 4 - Primary analysis of DGEs through various plots .. ")
       
    +     if(!dir.exists(paste(backup_file, "/dge_plots/", sep = ""))){
      +         
        +         dir.create(paste(backup_file, "/dge_plots/", sep = ""))
      +         
        +         
        +         print(paste("Created a new directory to save the DGE analysis results",sep = ""))
      +         
        +     }else{print(paste("Found backup directory to save the DGE analysis",sep=""))}
 }

upreg_genes_num = dim(dge_tophits[dge_tophits$logFC > 0.001 & dge_tophits$P.Value < 0.05,])[1]
downreg_genes_num = dim(dge_tophits[dge_tophits$logFC < -0.001 & dge_tophits$P.Value < 0.05,])[1]
de_genes = data.frame(genes = c(upreg_genes_num, downreg_genes_num), type = c("Upregulated", "Downregulated"))
de_genes_plot = ggplot(de_genes, aes(y = type, x = as.numeric(genes), fill = type))+
  +     geom_bar(position = "dodge", stat = "identity", width = 0.4)+
  +     ylab("Type")+xlab("Number of genes")+geom_text(aes(label = genes), position = position_dodge(width = 0.9), hjust = -0.25, size = 6)+
  +     scale_fill_manual("Groups", values = c("#0A0DAE", "#F57A00"))+ ggtitle("Number of DEGs (pval < 0.05, lfc > 0.1)")+
  +     theme_bw()+theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title.y = element_blank(), axis.title.x = element_text(size = 16),
                         +                      legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                         +                      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 14))
ggsave(plot = de_genes_plot, filename = paste( "de_genes.jpeg", sep = ""),
         +        dpi = 300, height = 6, width = 15, device = "jpeg", bg = "white")

#performing clustering of samples using PCA"
pc_comp = prcomp(t(expr_mtx))
 pc_var = pc_comp$sdev^2
 View(pc_comp)
 pc_per = round(pc_var/sum(pc_comp$sdev^2)*100, 1)
 pc_comp_df = as_tibble(pc_comp$x)
 groups = metainfo$group
 #dont use ggtitle
 pca_plot = ggplot(pc_comp_df, aes(x=PC1, y=PC2, colour = groups))+
   +   geom_point(size=4, alpha = 0.7)+stat_ellipse(alpha = 0.5)+
   +   theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                        +                    legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                        +                    plot.title = element_text(size = 18, hjust = 0.5))
 ggsave(plot = pca_plot, filename = "pca_plot.jpeg",
          +        dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
 #t-SNE 
 # Load libraries
 library(Rtsne)
 library(ggplot2)
 
 # Run PCA first (recommended for stability before t-SNE)
 pc_comp <- prcomp(t(expr_mtx[,-1]), scale. = TRUE)
 pc_scores <- pc_comp$x[, 1:30]   # take first 30 PCs for t-SNE
 
 # Run t-SNE
 set.seed(123)  # ensures reproducibility
 tsne_out <- Rtsne(pc_scores, dims = 2, perplexity = 30, verbose = TRUE)
 
 # Put results in a dataframe
 tsne_df <- as_tibble(tsne_out$Y)
 colnames(tsne_df) <- c("tSNE1", "tSNE2")
 tsne_df$groups <- metainfo$group
 
 # Plot (similar style to your PCA plot)
 tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, colour = groups)) +
   geom_point(size = 4, alpha = 0.7) +
   stat_ellipse(alpha = 0.5) +
   theme_bw() +
   theme(
     axis.title.x = element_text(size = 16),
     axis.title.y = element_text(size = 16),
     legend.title = element_text(size = 16),
     legend.text  = element_text(size = 14),
     plot.title   = element_text(size = 18, hjust = 0.5)
   )
 
 # Save the plot
 ggsave(
   plot = tsne_plot,
   filename = "tsne_plot.jpeg",
   dpi = 300,
   height = 10,
   width = 12,
   device = "jpeg",
   bg = "white"
 )
 
 
 #umap code
 umap_comp = umap::umap(t(expr_mtx[,-1]))
 
 #extracting the coordinates of the samples in the umap
 umap_layout = umap_comp[["layout"]]
 umap_layout = data.frame(umap_layout)
 
 #adding information on the sampleIDs as a column
 umap_layout$sample_id = rownames(umap_layout)
 umap_layout = subset(umap_layout, select = -c(sample_id))
 groups = metainfo$T0
 #plotting the umap now
 umap_plot = ggplot(umap_layout, aes(x = X1, y = X2, color = groups)) + geom_point(size = 4, alpha = 0.7)+stat_ellipse(alpha = 0.5)+
   xlab("X1")+ylab("X2")+
   scale_colour_manual("Groups", values = c("#EC0B43", "#58355E"))+
   theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                    legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                    plot.title = element_text(size = 18, hjust = 0.5))
 
 #saving the plot
 ggsave(plot = umap_plot, filename = "umap_plot.jpeg",
        dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
 
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
    + )
ggsave(plot = volcano_plot,
       +        filename = "volcano_plot.jpeg",
       +        dpi = 300,
       +        height = 10,
       +        width = 12,
       +        device = "jpeg",
       +        bg = "white")

#performing heatmap plot to visualize the DEGs
upreg_geneinfo = dge_tophits[dge_tophits$logFC > 0 & dge_tophits$P.Value < 0.05,] #setting p-val cutoff
upreg_geneinfo = upreg_geneinfo[order(upreg_geneinfo$logFC, decreasing = TRUE),]
upreg_geneinfo = upreg_geneinfo$gene_Symbol[1:5]
downreg_geneinfo = dge_tophits[dge_tophits$logFC < 0 & dge_tophits$P.Value < 0.05,]
downreg_geneinfo = downreg_geneinfo[order(downreg_geneinfo$logFC, decreasing = FALSE),]
downreg_geneinfo = downreg_geneinfo[order(downreg_geneinfo$logFC, decreasing = FALSE),]
downreg_geneinfo = downreg_geneinfo$gene_Symbol[1:5]
diff_exp_genes_subset = expr_mtx[rownames(expr_mtx) %in% c(upreg_genes$gene_symbol, downreg_genes$gene_symbol),]
 clust_cols = hclust(as.dist(1-cor(diff_exp_genes_subset, method = "spearman")), method = "complete")
 clust_rows = hclust(as.dist(1-cor(t(diff_exp_genes_subset), method = "pearson")), method = "complete")
   module_assign = cutree(clust_rows, k=2)
module_colour = rainbow(n = length(unique(module_assign)), start = 0.1, end = 0.9)   
 module_colour = module_colour[as.vector(module_assign)]   
 col_pal = colorRampPalette(colors = c("red", "white", "green"))(10)   
 my_gene_col = data.frame(cluster = ifelse(test = module_assign == 1, yes = "upregulated", no = "downregulated")) 
 my_sample_col = metainfo[,c("ID", "T0", "sample_info")]
 my_sample_col = subset(my_sample_col, select = c(T0))
 heatmap = pheatmap(mat = diff_exp_genes_subset, 
                    +                    legend = TRUE,
                    +                    RowSideColors = module_colour, 
                    +                    col = col_pal, 
                    +                    scale = "row", 
                    +                    cluster_cols = FALSE, 
                    +                    cluster_rows = TRUE,
                    +                    labRow = NA, 
                    +                    density.info = "none", 
                    +                    trace = "none", 
                    +                    cexRow = 1, 
                    +                    keysize = 2, 
                    +                    annotation_row = my_gene_col, 
                    +                    annotation_col = my_sample_col, 
                    +                    labels_col = rep("", dim(my_sample_col)[1]), 
                    +                    cutree_rows = 2)
 expr_mtx = read.delim("GSE298577_Prevnar.D0.D7.RMA.txt")
 
 
 
 
 ##Identifying DGE mods using tmod
 install.packages("tmod")
 library(tmod)
 library(ggplot2)
 dge_tophits = dge_tophits[order(dge_tophits$adj.P.Val),]
 rownames(dge_tophits) = NULL
 cerno_results = tmodCERNOtest(dge_tophits$gene_symbol, order.by = "pval")
 write.table(file="cerno_results.txt",x = cerno_results,sep = "\t",col.names = TRUE, quote = FALSE)
 if(dim(cerno_results)[1] >= 2) {
         cerno_results = cerno_results[order(cerno_results$adj.P.Val),]
         cerno_results = cerno_results[1:30,]
         sig_genes = tmodDecideTests(g = dge_tophits$gene_symbol, lfc = dge_tophits$logFC, pval = dge_tophits$P.Value,
                                                                        lfc.thr = 0.01, pval.thr = 0.05)
         names(sig_genes) = "DGE_mods"
         tmod_plots = ggPanelplot(res = list(DGE_mods = cerno_results),sgenes = sig_genes)
         tmod_plots = tmod_plots+
             theme(plot.title = element_text(size = 18, hjust = 0.5), axis.text.y =element_text(size = 14), axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
                                  legend.title = element_text(size = 14), legend.text = element_text(size = 12))
         
          #saving the plot
          ggsave(plot = tmod_plots, filename = "cerno_resultes", device = "png", 
                           width = 14, height = 9, dpi = 300)
         
           
         }else{print("cerno test failed - less than 2 modules")}
View(cerno_results)

tmod_plots <- ggPanelplot(res = list(DGE_mods = cerno_results), sgenes = sig_genes) +
   +      theme(
     +        plot.title = element_text(size = 18, hjust = 0.5),
     +        axis.text.y = element_text(size = 14),
     +        axis.text.x = element_text(size = 12),
     +        axis.title.x = element_text(size = 14),
     +        legend.title = element_text(size = 14),
     +        legend.text = element_text(size = 12)    )

ggsave(
   +     plot = tmod_plots, filename = "cerno_results.png",
   +     width = 14, height = 9, dpi = 300)


# performing analysis of enrichment of disease pathways
dge_tophits = dge_tophits[dge_tophits$P.Value < 0.05 & abs(dge_tophits$logFC) > 0.3,]
regulated_entrezids = mapIds(clariomshumantranscriptcluster.db, keys = dge_tophits$gene_symbol,
                             column = 'ENTREZID', keytype = 'SYMBOL')
regulated_entrezids = as.data.frame(regulated_entrezids)
regulated_entrezids$gene_symbol = rownames(regulated_entrezids)
rownames(regulated_entrezids) = NULL
regulated_entrezids = regulated_entrezids[!is.na(regulated_entrezids$regulated_entrezids),]

if(dim(regulated_entrezids)[1] > 1){
  
  dge_tophits = merge(dge_tophits, regulated_entrezids, by = "gene_symbol")
  
  regulated_entrezids_logFC_vector = dge_tophits$logFC
  names(regulated_entrezids_logFC_vector) = dge_tophits$regulated_entrezids

  regulated_entrezids_logFC_vector = sort(regulated_entrezids_logFC_vector, decreasing = TRUE)}


# Making barplot of the distribution of genes in the disease gene networks (Install DOSE)

enrich_disgen = enrichDGN(dge_tophits$regulated_entrezids)
barplot = barplot(enrich_disgen)
ggsave(filename = "barplot_disgen_network_barplot.jpg", plot = barplot,
       dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
dev.off()
 



# Making the network of disease and the associated genes
enrich_design_genenames = setReadable(enrich_disgen, 'org.Hs.eg.db', 'ENTREZID')
cnet_plot = cnetplot(enrich_design_genenames, showCategory = 5, color.params = list(foldChange = regulated_entrezids_logFC_vector, edge = TRUE), 
                     circular = TRUE, layout = "kk")+ scale_color_gradient2(name='Fold Change', low='blue', high='red')

ggsave(filename = "circ_plot_dis_gene_network.jpg", plot = cnet_plot, 
       dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")


#Making the treeplot to cluster the disease condition on the gsea
enrich_pairwise_disease = enrichplot::pairwise_termsim(enrich_design_genenames)
if(length(enrich_design_genenames$ID) > 2){
  
  #making the treeplot
  treeplot = enrichplot::treeplot(enrich_pairwise_disease, geneClusterPanel = "dotplot", cluster.params = list(method = "average"))
  
  #saving with ggsave
  ggsave(filename = paste(backup_file, "/disease_annotate/treeplot_disease_cluster.jpg", sep = ""), plot = treeplot, 
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")}


# Performing gene set enrichment analysis for disease ontology
gene_set_enrich <- gseDO(
  regulated_entrezids_logFC_vector,
  pvalueCutoff = 0.2
)
if(length(gene_set_enrich$pvalue) > 0){
  dotplot = dotplot(gene_set_enrich)
  ggsave(filename = "dotplot.jpg", plot = dotplot, 
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")}





#Performing gene set enrichment analysis
#preaparing the gsea modules
c7_sign = msigdbr(species = "Homo sapiens", category = "C7")
c7_sign_subset = c7_sign[,c("gs_name", "gene_symbol")]
c7_vax_sign = msigdbr(species = "Homo sapiens", category = "C7", subcategory = "VAX")
c7_vax_sign_subset = c7_vax_sign[,c("gs_name", "gene_symbol")]
hallmark_sign = msigdbr(species = "Homo sapiens", category = "H")
hallmark_sign_subset = hallmark_sign[,c("gs_name", "gene_symbol")]
upreg_genes = dge_tophits$gene_symbol[dge_tophits$logFC > 0 & dge_tophits$P.Value < 0.05]
de_genes_df = data.frame(genes = c(upreg_genes,downreg_genes), lfc = dge_tophits$logFC[dge_tophits$gene_symbol %in% c(upreg_genes, downreg_genes)])
de_genes_info = de_genes_df$lfc
names(de_genes_info) = de_genes_df$genes
# performing gsea for c7 (general immunology module)")
c7_gsea = GSEA(de_genes_info, TERM2GENE = c7_sign_subset, pvalueCutoff = 0.05)
c7_gsea_result = c7_gsea@result

if(dim(c7_gsea_result)[1] > 2){
  
  #@print check
  print("(Step 8) - for c7 - found more than 2 gene sets hence making next analysis")
  
  #plotting the gsea results
  c7_gsea_plot = enrichplot::dotplot(object = c7_gsea, showCategory = 5)
  
  #beautifying the plot with ggplot
  c7_gsea_plot = c7_gsea_plot+
    theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                     axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     plot.title = element_text(size = 18, hjust = 0.5))
  
  #to the table, adding the description of the gene set term
  c7_gsea_result$description_entire = unlist(lapply(c7_gsea_result$ID, function(x){unique(c7_sign$gs_description[c7_sign$gs_name == x])}))
  
  #saving the plot
  ggsave(plot = c7_gsea_plot, filename = "c7_gsea_plot.jpeg",
         dpi = 300, height = 10, width = 10, device = "jpeg", bg = "white")
  
  #saving the tabulated results
  write.table(x = c7_gsea_result, file = "c7_gsea_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}



# performing gsea for c7 subcategory VAX

c7_vax_gsea = GSEA(de_genes_info, TERM2GENE = c7_vax_sign_subset, pvalueCutoff = 0.05)
c7_vax_gsea_result = c7_vax_gsea@result

if(dim(c7_vax_gsea_result)[1] > 2){
 
  c7_vax_gsea_plot = enrichplot::dotplot(object = c7_vax_gsea, showCategory = 5)
  
  c7_vax_gsea_plot = c7_vax_gsea_plot+
    theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                     axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     plot.title = element_text(size = 18, hjust = 0.5))
  
  c7_vax_gsea_result$description_entire = unlist(lapply(c7_vax_gsea_result$ID, function(x){unique(c7_vax_sign$gs_description[c7_vax_sign$gs_name == x])}))
  
  ggsave(plot = c7_vax_gsea_plot, filename ="c7_vax_gsea_plot.jpeg" ,
         dpi = 300, height = 12, width = 10, device = "jpeg", bg = "white")
  
  #saving the tabulated results
  write.table(x = c7_vax_gsea_result, file = "c7_vax_gsea_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}


# performing gsea for Hallmark genes

hallmark_gsea = GSEA(de_genes_info, TERM2GENE = hallmark_sign_subset, pvalueCutoff = 0.05)
hallmark_gsea_result = hallmark_gsea@result
if(dim(hallmark_gsea_result)[1] > 2){
  
  hallmark_gsea_plot = enrichplot::dotplot(object = hallmark_gsea, showCategory = 5)
  
  hallmark_gsea_plot = hallmark_gsea_plot+
    theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                     axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     plot.title = element_text(size = 18, hjust = 0.5))
  
 
  hallmark_gsea_result$description_entire = unlist(lapply(hallmark_gsea_result$ID, function(x){unique(hallmark_sign$gs_description[hallmark_sign$gs_name == x])}))
 
  ggsave(plot = hallmark_gsea_plot, filename = "hallmark_gsea_plot.jpeg",
         dpi = 300, height = 12, width = 10, device = "jpeg", bg = "white")

  write.table(x = c7_vax_gsea_result, file = "hallmark_gsea_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

#t-SNE plot
install.packages("Rtsne")
library(Rtsne)
library(ggplot2)
