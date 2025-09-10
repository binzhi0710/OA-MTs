library(data.table)

# Set current working directory
current_dir <- getwd()

# Read sample sizes from file
sample_sizes <- fread("samplesize.txt")

# Get GWAS files
gwas_files <- list.files(path = "~/S-MTAG/MTAGresult/", pattern = "_mtag_meta.txt$", full.names = TRUE)

# Generate R scripts for each GWAS file
for (gwas_file in gwas_files) {
  file_name <- basename(gwas_file)
  base_name <- sub("_mtag_meta.txt$", "", file_name)
  
  # Get sample size for Traitpair
  gwas_n <- sample_sizes$N[sample_sizes$Traitpair == base_name]
  
  # Create output directories
  folder_path <- file.path(current_dir, base_name)
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  cor_dir <- file.path(folder_path, "cor_matrix")
  dir.create(cor_dir, showWarnings = FALSE)
  
  # Generate R script content
  script_content <- sprintf('
library(ctwas)
library(data.table)
library(EnsDb.Hsapiens.v86)

gwas <- fread("%s")
z_snp <- gwas[, .(id = SNP, A1 = A1, A2 = A2, z = mtag_z)]

gwas_n <- %d

# Load region information
region_file <- "~/S-cTWAS/extdata/ldetect/EUR.b38.ldetect.regions.RDS"
region_info <- readRDS(region_file)
genome_version <- "b38"

# Set LD and SNP file paths
LD_dir <- "~/S-cTWAS/LD/"
LD_filestem <- sprintf("ukb_%%s_0.1_chr%%s.R_snp.%%s_%%s", genome_version, region_info$chrom, region_info$start, region_info$stop)
region_metatable <- region_info
region_metatable$LD_file <- file.path(LD_dir, paste0(LD_filestem, ".RDS"))
region_metatable$SNP_file <- file.path(LD_dir, paste0(LD_filestem, ".Rvar"))

# Create SNP and LD maps
res <- create_snp_LD_map(region_metatable)
region_info <- res$region_info
snp_map <- res$snp_map
LD_map <- res$LD_map

# Preprocess SNP data
z_snp <- preprocess_z_snp(z_snp, snp_map, 
                          drop_multiallelic = TRUE, 
                          drop_strand_ambig = TRUE,
                          varID_converter_fun = convert_to_ukb_varIDs)

# Load and preprocess weights
db_files <- list.files(path = "~/S-cTWAS/weights/eqtl/mashr/", pattern = "\\\\.db$", full.names = TRUE)
all_weights <- list()
for (db_file in db_files) {
  weight_name <- sub(".*mashr_(.*)\\\\.db$", "\\\\1_expression", db_file)
  context <- sub(".*mashr_(.*)\\\\.db$", "\\\\1", db_file)
  weight <- preprocess_weights(db_file,
                              region_info,
                              gwas_snp_ids = z_snp$id,
                              snp_map = snp_map,
                              type = "expression",
                              context = context,
                              weight_name = weight_name,
                              weight_format = "PredictDB",
                              drop_strand_ambig = TRUE,
                              scale_predictdb_weights = TRUE,
                              load_predictdb_LD = TRUE,
                              filter_protein_coding_genes = TRUE,
                              varID_converter_fun = convert_to_ukb_varIDs,
                              ncore = 6)
  all_weights[[weight_name]] <- weight
}
weights <- do.call(c, all_weights)

# Run cTWAS analysis
ctwas_res <- ctwas_sumstats(z_snp, 
                            weights, 
                            region_info, 
                            LD_map, 
                            snp_map, 
                            thin = 0.1,
                            group_prior_var_structure = "shared_type", 
                            filter_L = TRUE,
                            filter_nonSNP_PIP = FALSE,
                            min_nonSNP_PIP = 0.5,
                            min_abs_corr = 0.1, 
                            ncore = 6, 
                            ncore_LD = 4,
                            save_cor = TRUE,
                            cor_dir = "%s",
                            force_compute_cor = FALSE)

# Save results
saveRDS(ctwas_res, file = file.path("%s", "ctwas_res.RDS"))
', gwas_file, gwas_n, cor_dir, folder_path)
  
  # Save R script
  script_file <- file.path(folder_path, paste0(base_name, ".R"))
  writeLines(script_content, con = script_file)
}

# To execute all generated scripts, run in the terminal:
# for dir in */; do if [ -f "$dir"/*.R ]; then Rscript "$dir"/*.R; fi; done