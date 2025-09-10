# Load required packages
library(readxl)
library(dplyr)
library(data.table)

# Read Excel file
loci_data <- read_excel("1030metabo_GWAS-PW_loci.xlsx")

# Ensure phen1 column maintains two-digit format
loci_data$phen1 <- sprintf("%02d", loci_data$phen1)

# Group data by source column
grouped_sources <- split(loci_data, loci_data$source)

# Create directory for saving R scripts
dir.create("R_scripts1030", showWarnings = FALSE)

# Iterate through each group and generate corresponding R scripts
for (source_name in names(grouped_sources)) {
  source_data <- grouped_sources[[source_name]]
  
  # Build script content
  script_content <- paste0(
    "# Load required packages\n",
    "library(dplyr)\n",
    "library(data.table)\n\n",
    
    "# Read source data\n",
    "source_data <- data.frame(\n",
    "  chr = c(", paste(source_data$chr, collapse = ", "), "),\n",
    "  start = c(", paste(source_data$start, collapse = ", "), "),\n",
    "  stop = c(", paste(source_data$stop, collapse = ", "), "),\n",
    "  phen1 = '", source_data$phen1[1], "',\n",
    "  phen2 = '", source_data$phen2[1], "'\n",
    ")\n\n",
    
    "# Process phen1\n",
    "phen1_file <- paste0('~/07metabo_OA/1030remake/', source_data$phen1[1], '.txt')\n",
    "phen1 <- fread(phen1_file)\n",
    "phen1 <- phen1 %>% dplyr::select(SNP, chr, pos, Z, se)\n",
    "colnames(phen1) <- c('SNP', 'CHR', 'POS', paste0('Z_', source_data$phen1[1]), paste0('V_', source_data$phen1[1]))\n",
    "phen1$CHR <- paste0('chr', phen1$CHR)\n\n",
    
    "# Process phen2\n",
    "phen2_file <- paste0('~/07metabo_OA/00OA_noMHC/', source_data$phen2[1], '.txt')\n",
    "phen2 <- fread(phen2_file)\n",
    "phen2 <- phen2 %>% dplyr::select(SNP, chr, pos, Z, se)\n",
    "colnames(phen2) <- c('SNP', 'CHR', 'POS', paste0('Z_', source_data$phen2[1]), paste0('V_', source_data$phen2[1]))\n",
    "phen2$CHR <- paste0('chr', phen2$CHR)\n\n",
    
    "# Extract SNPs within loci ranges\n",
    "phen1_filtered <- data.frame()\n",
    "phen2_filtered <- data.frame()\n",
    "for (i in 1:nrow(source_data)) {\n",
    "  chr <- source_data$chr[i]\n",
    "  start <- source_data$start[i]\n",
    "  stop <- source_data$stop[i]\n",
    "  temp_phen1 <- phen1 %>% filter(CHR == paste0('chr', chr), POS >= start, POS <= stop)\n",
    "  phen1_filtered <- rbind(phen1_filtered, temp_phen1)\n",
    "  temp_phen2 <- phen2 %>% filter(CHR == paste0('chr', chr), POS >= start, POS <= stop)\n",
    "  phen2_filtered <- rbind(phen2_filtered, temp_phen2)\n",
    "}\n\n",
    
    "# Merge data for both phenotypes\n",
    "merged_data <- merge(phen1_filtered, phen2_filtered, by = 'SNP', all = TRUE)\n",
    "merged_data <- distinct(merged_data, SNP, .keep_all = TRUE)\n",  # Remove duplicate SNP rows, keeping the first occurrence
    "merged_data <- na.omit(merged_data)\n",
    
    "merged_data <- dplyr::select(merged_data, SNP, CHR.x, POS.x, paste0('Z_', source_data$phen1[1]), paste0('V_', source_data$phen1[1]), paste0('Z_', source_data$phen2[1]), paste0('V_', source_data$phen2[1]))\n",
    "colnames(merged_data) <- c('SNPID', 'CHR', 'POS', paste0('Z_', source_data$phen1[1]), paste0('V_', source_data$phen1[1]), paste0('Z_', source_data$phen2[1]), paste0('V_', source_data$phen2[1]))\n\n",
    
    "# Remove points at the edges of loci\n",
    "for (i in 1:nrow(source_data)) {\n",
    "  chr <- source_data$chr[i]\n",
    "  start <- source_data$start[i]\n",
    "  stop <- source_data$stop[i]\n",
    "  merged_data <- merged_data[!(merged_data$CHR == paste0('chr', chr) & (merged_data$POS == start | merged_data$POS == stop)), ]\n",
    "}\n\n",
    
    "# Sort and save\n",
    "merged_data$CHR_numeric <- as.numeric(gsub('chr', '', merged_data$CHR))\n",
    "merged_data <- merged_data[order(merged_data$CHR_numeric, merged_data$POS),]\n",
    "merged_data <- distinct(merged_data, CHR, POS, .keep_all = TRUE)\n",
    "merged_data$CHR_numeric <- NULL\n",
    "output_file <- paste0('GWASfile1030/', '", source_name, "', '.txt.gz')\n",
    "fwrite(merged_data, output_file, compress = 'gzip', sep = '\\t')  # Use tab as separator\n\n",
    
    "# Generate and sort bed file\n",
    "bed_data <- source_data %>% dplyr::select(chr, start, stop)\n",
    "bed_data$chr <- paste0('chr', bed_data$chr)\n",
    "bed_data$CHR_numeric <- as.numeric(gsub('chr', '', bed_data$chr))\n",
    "bed_data <- bed_data[order(bed_data$CHR_numeric, bed_data$start),]\n",
    "bed_data$CHR_numeric <- NULL\n",
    "bed_output_file <- paste0('bedfile1030/', '", source_name, "', '.bed')\n",
    "fwrite(bed_data, bed_output_file, col.names = FALSE, sep = '\\t')\n"
  )
  
  # Save each R script
  script_file <- paste0("R_scripts1030/", source_name, "_script.R")
  writeLines(script_content, script_file)
}