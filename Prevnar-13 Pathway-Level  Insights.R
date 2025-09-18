# performing analysis of enrichment of disease pathways
#Filter significant DEGs
#This script takes your significant DEGs, converts them from gene symbols to Entrez IDs, merges them with logFC values, and creates a ranked logFC vector ready for pathway/disease enrichment analysis.
#Keeps only genes with: p-value < 0.05 (statistically significant) ; |log fold change| > 0.3 (biologically meaningful effect size)
dge_tophits = dge_tophits[dge_tophits$P.Value < 0.05 & abs(dge_tophits$logFC) > 0.3,]

#Map gene symbols to Entrez IDs
regulated_entrezids = mapIds(
  clariomshumantranscriptcluster.db, 
  keys = dge_tophits$gene_symbol,
  column = 'ENTREZID', 
  keytype = 'SYMBOL'
)

#Clean and reformat Entrez ID table
regulated_entrezids = as.data.frame(regulated_entrezids)
regulated_entrezids$gene_symbol = rownames(regulated_entrezids)
rownames(regulated_entrezids) = NULL
regulated_entrezids = regulated_entrezids[!is.na(regulated_entrezids$regulated_entrezids),]

dge_tophits = merge(dge_tophits, regulated_entrezids, by = "gene_symbol") # Create a named vector for enrichment
regulated_entrezids_logFC_vector = sort(regulated_entrezids_logFC_vector, decreasing = TRUE) #Sort by fold change

enrich_disgen = enrichDGN(dge_tophits$regulated_entrezids) #Run enrichment analysis
barplot = barplot(enrich_disgen)
ggsave(filename = "barplot_disgen_network_barplot.jpg", plot = barplot,
       dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
dev.off()


#building visualizations of disease–gene relationships and clustering disease terms.
enrich_design_genenames = setReadable(enrich_disgen, 'org.Hs.eg.db', 'ENTREZID') #Convert enrichment results to readable gene names

#Build a circular cnetplot
cnet_plot = cnetplot(
    enrich_design_genenames, 
    showCategory = 5, 
    color.params = list(foldChange = regulated_entrezids_logFC_vector, edge = TRUE), 
    circular = TRUE, 
    layout = "kk"
) + scale_color_gradient2(name='Fold Change', low='blue', high='red')

ggsave(filename = "circ_plot_dis_gene_network.jpg", plot = cnet_plot, 
       dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white") #Save the cnetplot

enrich_pairwise_disease = enrichplot::pairwise_termsim(enrich_design_genenames) #Compute similarity between enriched disease terms

#Treeplot (clustering disease terms)
#Only runs if there are at least 3 enriched diseases (otherwise clustering is meaningless).
#Displays genes inside clusters as dotplots.
if(length(enrich_design_genenames$ID) > 2){
  
  treeplot = enrichplot::treeplot(
      enrich_pairwise_disease, 
      geneClusterPanel = "dotplot", 
      cluster.params = list(method = "average")
  )
#Save the treeplot
ggsave(filename = paste(backup_file, "/disease_annotate/treeplot_disease_cluster.jpg", sep = ""), plot = treeplot, 
       dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")



#Run GSEA for Disease Ontology
#This script runs GSEA for Disease Ontology terms, checks if any enrichments are found, then visualizes them in a dotplot and saves the result as a figure
gene_set_enrich <- gseDO(
  regulated_entrezids_logFC_vector,
  pvalueCutoff = 0.2
)
if(length(gene_set_enrich$pvalue) > 0){
dotplot = dotplot(gene_set_enrich)
ggsave(filename = "dotplot.jpg", plot = dotplot, 
       dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")





