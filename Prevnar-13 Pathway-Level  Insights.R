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

#Gene Set Enrichment Analysis (GSEA) with immunology- and vaccine-related signatures
#This pipeline tests whether your ranked gene list (up and down DEGs) is enriched for immune-related gene sets (C7), especially vaccine-related modules (C7-VAX).
#Preparing gene sets (modules)
#C7: Immunologic signatures (derived from immune-related experiments).
#C7-VAX: Subset of immunologic signatures specifically from vaccine studies.
#H is (Hallmark): Core biological processes (e.g., interferon response, proliferation, apoptosis).
#Each is subset to keep only the gene set name and gene symbol, which are required for GSEA input.
c7_sign = msigdbr(species = "Homo sapiens", category = "C7")
c7_sign_subset = c7_sign[,c("gs_name", "gene_symbol")]

c7_vax_sign = msigdbr(species = "Homo sapiens", category = "C7", subcategory = "VAX")
c7_vax_sign_subset = c7_vax_sign[,c("gs_name", "gene_symbol")]

hallmark_sign = msigdbr(species = "Homo sapiens", category = "H")
hallmark_sign_subset = hallmark_sign[,c("gs_name", "gene_symbol")]


#Preparing ranked gene list
upreg_genes = dge_tophits$gene_symbol[dge_tophits$logFC > 0 & dge_tophits$P.Value < 0.05]

de_genes_df = data.frame(
  genes = c(upreg_genes, downreg_genes),
  lfc = dge_tophits$logFC[dge_tophits$gene_symbol %in% c(upreg_genes, downreg_genes)]
)

de_genes_info = de_genes_df$lfc
names(de_genes_info) = de_genes_df$genes

#Runs GSEA using the immunology signatures (C7).
c7_gsea = GSEA(de_genes_info, TERM2GENE = c7_sign_subset, pvalueCutoff = 0.05)
c7_gsea_result = c7_gsea@result

  #Checking number of enriched sets
if(dim(c7_gsea_result)[1] > 2){
  print("(Step 8) - for c7 - found more than 2 gene sets hence making next analysis")
  c7_gsea_plot = enrichplot::dotplot(object = c7_gsea, showCategory = 5) #Makes a dotplot of the top 5 enriched C7 signatures.
c7_gsea_plot = c7_gsea_plot+
    theme_bw()+theme(...)
c7_gsea_result$description_entire = unlist(lapply(
  c7_gsea_result$ID, 
  function(x){unique(c7_sign$gs_description[c7_sign$gs_name == x])}
))
ggsave(plot = c7_gsea_plot, filename = "c7_gsea_plot.jpeg",
       dpi = 300, height = 10, width = 10, device = "jpeg", bg = "white")

  




#vaccine-specific extension of your GSEA workflow

c7_vax_gsea = GSEA(de_genes_info, TERM2GENE = c7_vax_sign_subset, pvalueCutoff = 0.05)
c7_vax_gsea_result = c7_vax_gsea@result #Compares them against C7-VAX modules (gene sets specifically derived from vaccine studies in MSigDB).
  if(dim(c7_vax_gsea_result)[1] > 2){
c7_vax_gsea_plot = enrichplot::dotplot(object = c7_vax_gsea, showCategory = 5)
c7_vax_gsea_plot = c7_vax_gsea_plot + theme_bw() + theme(...)
c7_vax_gsea_result$description_entire = unlist(lapply(
  c7_vax_gsea_result$ID, 
  function(x){unique(c7_vax_sign$gs_description[c7_vax_sign$gs_name == x])}
))
ggsave(plot = c7_vax_gsea_plot, filename ="c7_vax_gsea_plot.jpeg", ...)
write.table(x = c7_vax_gsea_result, file = "c7_vax_gsea_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#This script identifies which known vaccine response gene sets (from MSigDB C7-VAX) are significantly enriched in your DEGs. 
#Then it visualizes the strongest ones in a dotplot and saves both the plot and a detailed table.




#GSEA with Hallmark gene sets
    #This workflow checks whether your differentially expressed genes are enriched in Hallmark pathways (e.g. interferon signaling, TNF-alpha response, apoptosis, glycolysis). 
    #This gives you a high-level overview of biological programs activated or suppressed by vaccination.
    hallmark_gsea = GSEA(de_genes_info, TERM2GENE = hallmark_sign_subset, pvalueCutoff = 0.05)
hallmark_gsea_result = hallmark_gsea@result

    if(dim(hallmark_gsea_result)[1] > 2){
      hallmark_gsea_plot = enrichplot::dotplot(object = hallmark_gsea, showCategory = 5)
hallmark_gsea_plot = hallmark_gsea_plot + theme_bw() + theme(...)
hallmark_gsea_result$description_entire = unlist(lapply(
  hallmark_gsea_result$ID,
  function(x){unique(hallmark_sign$gs_description[hallmark_sign$gs_name == x])} #Add pathway descriptions
))
      ggsave(plot = hallmark_gsea_plot, filename = "hallmark_gsea_plot.jpeg", ...)
write.table(x = c7_vax_gsea_result, file = "hallmark_gsea_table.txt", ...)
write.table(x = hallmark_gsea_result, file = "hallmark_gsea_table.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





