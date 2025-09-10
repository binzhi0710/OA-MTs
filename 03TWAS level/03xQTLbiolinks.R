library(stringr)
library(glue)

# Get GWAS files and tissues
gwas_files <- list.files("D:/07metaboOA/1030MTAGresult", pattern = "\\.txt$", full.names = TRUE)
tissues <- xQTLbiolinks::tissueSiteDetailGTExv8$tissueSiteDetail %>% unique() %>% .[1:54]

# Create output directory for R scripts
output_dir <- "R_scripts/"
dir.create(output_dir, showWarnings = FALSE)

# Generate R scripts
script_number <- 1
for (gwas_file in gwas_files) {
  for (tissue in tissues) {
    gwas_name <- str_replace(basename(gwas_file), "_mtag_meta.txt", "")
    tissue_clean <- str_replace_all(tissue, " ", "-")
    output_file_name <- paste0(gwas_name, "_", tissue_clean, ".tsv")
    script_file_name <- str_replace(output_file_name, "\\.tsv$", ".R")
    
    # Generate R script content
    script_content <- glue("
library(xQTLbiolinks)
library(coloc)
library(hyprcoloc)
library(data.table)
library(stringr)
library(R.utils)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)

output_file_path <- 'result/{output_file_name}'

# Skip if output file exists
if (file.exists(output_file_path)) {{
  message('File {output_file_name} already exists. Skipping...')
  quit(save = 'no')
}}

# Read GWAS data
gwasDF <- fread('{gwas_file}')
gwasDF <- gwasDF[str_detect(SNP, '^rs'), .(rsid = SNP, chrom = CHR, position = BP, pValue = mtag_pval, AF = meta_freq, beta = mtag_beta, se = mtag_se)]
gwasDF$pValue[gwasDF$pValue == 0] <- 1e-300

# Get sentinel SNPs
sentinelSnpDF <- xQTLanalyze_getSentinelSnp(gwasDF, centerRange = 1e6, genomeVersion = 'grch37', grch37To38 = TRUE)

# Get traits data
traitsAll <- xQTLanalyze_getTraits(sentinelSnpDF, detectRange = 1e6, tissueSiteDetail = '{tissue}')

# Get gene information
genesAll <- xQTLquery_gene(unique(traitsAll$gencodeId))

# Set index for faster merging
setindex(gwasDF, rsid)

# Initialize results
colocResultAll <- data.table()

# Process each gene
for (i in 1:nrow(genesAll)) {{
  eQTL_i <- xQTLdownload_eqtlAllAsso(genesAll[i]$gencodeId, geneType = 'gencodeId', tissueLabel = '{tissue}', withB37VariantId = FALSE, data_source = 'liLab')
  gwasDF_i <- gwasDF[rsid %in% eQTL_i$rsid, ]
  if (is.null(eQTL_i) || nrow(gwasDF_i) == 0) {{ next() }}
  
  # Merge eQTL and GWAS data
  eQTL_i <- merge(eQTL_i, gwasDF_i[, .(rsid, chrom, position)], by = 'rsid')[, .(rsid, chrom, position, pValue, maf, beta, se)]
  
  # Perform COLOC analysis
  colocResult_i <- xQTLanalyze_coloc_diy(gwasDF = gwasDF_i, qtlDF = eQTL_i, method = 'Both')
  
  if (!is.null(colocResult_i)) {{ 
    colocResult_i <- colocResult_i$coloc_Out_summary
    colocResult_i <- cbind(genesAll[i, c('geneSymbol', 'gencodeId')], colocResult_i)
    colocResultAll <- rbind(colocResultAll, colocResult_i)
  }}
  message(' == Id: ', i, '/', nrow(genesAll), ' == Gene:', genesAll[i]$gencodeId)
}}

# Merge gene information
outGenes <- xQTLquery_gene(colocResultAll$gencodeId)
outGenes <- merge(colocResultAll, outGenes[, c('geneSymbol', 'gencodeId', 'entrezGeneId', 'geneType')], by = 'gencodeId', sort = FALSE)

# Save results
fwrite(outGenes, file = output_file_path, sep = '\\t')
")
    
    # Save R script
    script_path <- file.path(output_dir, script_file_name)
    writeLines(script_content, script_path)
    
    script_number <- script_number + 1
  }
}

message("All scripts generated successfully!")

# To execute all generated scripts, run in the terminal:
# for script in R_scripts/*.R; do Rscript "$script"; done