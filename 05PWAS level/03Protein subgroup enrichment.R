library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(openxlsx)
library(aPEAR)
library(data.table)

# Set directories
mcode_dir <- "MCODE/"
output_dir <- "enrichment_results/"

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

# Get all .txt files
mcode_files <- list.files(mcode_dir, pattern = "\\.txt$", full.names = TRUE)

# Process each file
for (file_path in mcode_files) {
  gene_symbols <- readLines(file_path)
  
  # Convert to Entrez IDs
  gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Perform enrichment analyses
  enrich_kegg <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  enrich_go <- enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", 
                        pAdjustMethod = "BH", qvalueCutoff = 0.05)
  enrich_reactome <- enrichPathway(gene = gene_entrez$ENTREZID, organism = "human", 
                                   pAdjustMethod = "BH", qvalueCutoff = 0.05)
  enrich_wiki <- enrichWP(gene = gene_entrez$ENTREZID, organism = "Homo sapiens", 
                          pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  # Extract results
  kegg_result <- enrich_kegg@result
  go_result <- enrich_go@result
  reactome_result <- enrich_reactome@result
  wiki_result <- enrich_wiki@result
  
  # Find common columns
  common_cols <- Reduce(intersect, list(colnames(kegg_result), colnames(go_result), 
                                        colnames(reactome_result), colnames(wiki_result)))
  
  # Retain common columns
  kegg_result <- kegg_result[, common_cols, drop = FALSE]
  go_result <- go_result[, common_cols, drop = FALSE]
  reactome_result <- reactome_result[, common_cols, drop = FALSE]
  wiki_result <- wiki_result[, common_cols, drop = FALSE]
  
  # Merge and filter results
  all_results <- rbind(kegg_result, go_result, reactome_result, wiki_result)
  all_results <- all_results[!duplicated(all_results$Description) & all_results$qvalue < 0.01, ]
  
  # Save to CSV
  base_name <- tools::file_path_sans_ext(basename(file_path))
  output_file <- file.path(output_dir, paste0(base_name, "_enrichment.csv"))
  fwrite(all_results, output_file)
}
