library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(openxlsx)
library(parallel)
library(foreach)
library(doParallel)

if (!dir.exists("Druggenes")) dir.create("Druggenes")
if (!dir.exists("Drugpathway")) dir.create("Drugpathway")

get_genes_for_drug <- function(drug_name, dgidb_drugs, drugcentral_drugs, pharmgkb_drugs) {
  genes <- c(
    dgidb_drugs %>% filter(drug == drug_name) %>% pull(gene),
    drugcentral_drugs %>% filter(drug == drug_name) %>% pull(gene),
    pharmgkb_drugs %>% filter(drug == drug_name) %>% pull(gene)
  )
  unique(genes)
}

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

read_gene_drug_data <- function(folder_path) {
  file_path <- file.path(folder_path, "gene-drug.xlsx")
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  read_excel(file_path) %>% mutate(drug = tolower(drug))
}

dgidb_drugs <- read_gene_drug_data("DGIdb")
drugcentral_drugs <- read_gene_drug_data("Drugcentral")
pharmgkb_drugs <- read_gene_drug_data("PharmGKB")

if (!is.data.frame(dgidb_drugs) || !is.data.frame(drugcentral_drugs) || !is.data.frame(pharmgkb_drugs)) {
  stop("Gene-drug data not loaded as data frames.")
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(i = seq_along(approved_drugs_df$drug), .packages = c("dplyr", "clusterProfiler", "org.Hs.eg.db", "openxlsx")) %dopar% {
  drug_name <- approved_drugs_df$drug[i]
  drug_genes <- get_genes_for_drug(drug_name, dgidb_drugs, drugcentral_drugs, pharmgkb_drugs)
  drug_gene_file <- sprintf("Druggenes/%04d_%s.txt", i, drug_name)
  writeLines(drug_genes, drug_gene_file)
  if (length(drug_genes) > 0) {
    pathway_results <- perform_enrichment(drug_genes)
    pathway_file <- sprintf("Drugpathway/%04d_%s.xlsx", i, drug_name)
    write.xlsx(pathway_results, pathway_file, row.names = FALSE)
  }
  return(NULL)
}

stopCluster(cl)

cat("Gene lists and enrichment analysis results for all drugs completed and saved.\n")