# Set working directory
setwd("~/S-CARMA/")
library(data.table)
library(magrittr)
library(dplyr)
library(R.utils)
library(arrow)
library(readxl)

# Read loci file
loci <- read_excel("GenomicRiskLoci_CARMASI0729.xlsx")

# Read bim data
bim_data <- fread("EUR2m.bim", select = 2, col.names = "SNP")

# Check column types and convert to numeric if necessary
loci$chr <- as.numeric(loci$chr)
loci$start <- as.numeric(loci$start)
loci$end <- as.numeric(loci$end)

# Iterate through all rows of loci file
for (i in seq_len(nrow(loci))) {
  # Extract information for the current row
  traitpair <- loci$traitpair[i]
  chr <- loci$chr[i]
  start <- loci$start[i]
  end <- loci$end[i]
  genomic_locus <- loci$GenomicLocus[i]
  
  # Dynamically generate MTAG file path
  file_path <- paste0("~/02SI/", traitpair, ".txt")
  data <- fread(file_path)
  
  # Filter SNPs within CHR and BP range
  sumstat <- data[chr == chr & pos >= start - 250000 & pos <= end + 250000]
  
  # Rename columns
  setnames(sumstat, old = c("SNP", "chr", "pos", "A1", "A2", "eaf", "beta", "se", "pval"),
           new = c("SNP", "CHR", "BP", "A1", "A2", "eaf", "beta", "se", "pval"))
  
  # Dynamically select annotation file based on chromosome number
  parquet_file <- paste0("baselineLF2.2.UKB/baselineLF2.2.UKB.", chr, ".annot.parquet")
  parquet_data <- read_parquet(parquet_file)
  
  # Find common SNPs across datasets
  common_snps <- Reduce(intersect, list(sumstat$SNP, bim_data$SNP, parquet_data$SNP))
  
  sumstat_filtered <- sumstat[sumstat$SNP %in% common_snps, ]
  parquet_data_filtered <- parquet_data[parquet_data$SNP %in% common_snps, ]
  
  sumstat_filtered <- sumstat_filtered[order(sumstat_filtered$BP), ]

  snp_order <- sumstat_filtered$SNP
  parquet_data_filtered <- parquet_data_filtered[match(snp_order, parquet_data_filtered$SNP), ]
  
  subset_data <- sumstat_filtered[, c("CHR", "BP", "SNP", "A1", "A2")]
  
  # Dynamically generate base filename based on traitpair and GenomicLocus
  base_filename <- paste0(traitpair, "_locus", genomic_locus)
  

  fwrite(subset_data, file = paste0("./01polyfun/", base_filename, "_polyfun.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  fwrite(parquet_data_filtered, file = paste0("./02annotation/", base_filename, "_annot.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  fwrite(sumstat[, c("SNP", "A1")], file = paste0("./03plink/", base_filename, "_plink.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Print progress
  cat("Processed:", traitpair, "for Genomic Locus", genomic_locus, "\n")
}