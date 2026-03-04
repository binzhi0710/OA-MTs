library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(openxlsx)
library(parallel)
library(foreach)
library(doParallel)

input_folder <- "Leastplei4"
output_folder <- "TraitEnrichment"
if (!dir.exists(output_folder)) dir.create(output_folder)

perform_enrichment <- function(gene_list) {
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  kegg <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1)
  kegg_results <- as.data.frame(kegg)
  
  go <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1)
  go_results <- as.data.frame(go)
  
  combined_results <- bind_rows(
    kegg_results %>% mutate(Source = "KEGG"),
    go_results %>% mutate(Source = "GO")
  )
  
  combined_results <- combined_results %>%
    arrange(p.adjust, pvalue) %>%
    slice_head(n = 50)
  
  return(combined_results)
}

file_list <- list.files(input_folder, pattern = "\\.txt$", full.names = TRUE)

num_cores <- 40
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(file = file_list, .packages = c("clusterProfiler", "org.Hs.eg.db", "dplyr", "openxlsx")) %dopar% {
  file_name <- tools::file_path_sans_ext(basename(file))
  
  gene_list <- readLines(file)
  
  if (length(gene_list) == 0) {
    return(NULL)
  }
  
  enrichment_results <- perform_enrichment(gene_list)
  
  output_file <- file.path(output_folder, paste0(file_name, ".xlsx"))
  write.xlsx(enrichment_results, output_file, row.names = FALSE)
  
  return(file_name)
}

stopCluster(cl)

cat("Enrichment analysis completed for all files. Results saved in:", output_folder, "\n")