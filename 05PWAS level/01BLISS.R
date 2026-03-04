library(data.table)
library(BEDMatrix)

# Set directories
sumstats_dir <- "MTAGGWAS/"
output_dir <- "results/"

# Get sumstats files
sumstats_files <- list.files(sumstats_dir)

# Generate R scripts for each chromosome
for (i in 1:22) {
  script_content <- paste0(
    'source("BLISSAssociation_Support.R")\n',
    'sumstats_dir <- "', sumstats_dir, '"\n',
    'output_dir <- "', output_dir, '"\n',
    'sumstats_files <- list.files(sumstats_dir)\n',
    'weights_models_list <- c("UKB", "deCODE", "ARIC")\n',
    'sumstatfile <- sumstats_files[', i, ']\n',
    'for (weights_model in weights_models_list) {\n',
    '  args <- list(\n',
    '    sumstats = sumstatfile,\n',
    '    sumstats_dir = sumstats_dir,\n',
    '    weights_models = weights_model,\n',
    '    output_dir = output_dir,\n',
    '    output_name = paste0(tools::file_path_sans_ext(sumstatfile), "_", weights_model, ".txt")\n',
    '  )\n',
    '  sumstatfile <- args$sumstats\n',
    '  gwassum.dir <- args$sumstats_dir\n',
    '  weights.models <- args$weights_models\n',
    '  out.dir <- args$output_dir\n',
    '  output.name <- args$output_name\n',
    '  CHR.list <- c(1:22)\n',
    '  if (weights.models == "UKB") {\n',
    '    load("models/UKB_weights_info.RData")\n',
    '    ldref <- "1000G/1000G.EUR.usedSNP.QC.CHR"\n',
    '    weights.dir <- "models/UKB/"\n',
    '  } else if (weights.models == "deCODE") {\n',
    '    load("models/deCODE_weights_info.RData")\n',
    '    ldref <- "1000G/1000G.EUR.usedSNP.QC.CHR"\n',
    '    weights.dir <- "models/deCODE/"\n',
    '  } else if (weights.models == "ARIC") {\n',
    '    load("models/ARIC_weights_info.RData")\n',
    '    ldref <- "1000G/1000G.EUR.usedSNP.QC.CHR"\n',
    '    weights.dir <- "models/ARIC/"\n',
    '  }\n',
    '  if (is.na(sumstatfile) || is.na(gwassum.dir)) stop("Both --sumstats and --sumstats_dir must be specified.")\n',
    '  sumstats.org <- fread(paste0(gwassum.dir, "/", sumstatfile)) %>% as.data.frame()\n',
    '  header.inner <- tolower(colnames(sumstats.org))\n',
    '  header.inner[header.inner %in% c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps", "rsids")] <- "SNP"\n',
    '  header.inner[header.inner %in% c("zscore", "z-score", "gc_zscore", "z")] <- "Z"\n',
    '  header.inner[header.inner %in% c("a1", "allele1", "allele_1", "nea", "non_effect_allele")] <- "A1"\n',
    '  header.inner[header.inner %in% c("a2", "allele2", "allele_2", "effect_allele", "ea")] <- "A2"\n',
    '  header.inner[header.inner %in% c("chrom", "ch", "chr", "chromosome", "#chrom")] <- "CHROMOSOME"\n',
    '  header.inner[header.inner %in% c("n", "samplesize", "num_samples", "sample")] <- "N"\n',
    '  colnames(sumstats.org) <- header.inner\n',
    '  if (sum(grepl("chr", sumstats.org$CHROMOSOME)) <= 100) sumstats.org$CHROMOSOME <- paste0("chr", sumstats.org$CHROMOSOME)\n',
    '  proteins <- list.files(weights.dir)\n',
    '  outres <- as.data.frame(matrix(NA, 7000, 12))\n',
    '  colnames(outres) <- c("chr", "p0", "p1", "gene", "R2", "Zscore.classic", "p.classic", "beta_BLISS", "se_BLISS", "p_BLISS", "n_used_snp", "n_snp")\n',
    '  outindx <- 1\n',
    '  for (chromosome in CHR.list) {\n',
    '    reference.bim <- fread(paste0(ldref, chromosome, ".bim"), data.table = FALSE)\n',
    '    reference.bed <- BEDMatrix(paste0(ldref, chromosome), simple_names = TRUE)\n',
    '    sumstats <- sumstats.org[sumstats.org$CHROMOSOME == paste0("chr", chromosome), ]\n',
    '    sumstats <- sumstats[rowSums(is.na(sumstats)) == 0, ]\n',
    '    used.lookup <- lookup[lookup$chr == chromosome, ]\n',
    '    used.proteins <- proteins[proteins %in% paste0(used.lookup$gene, ".RData")]\n',
    '    reference.dup <- duplicated(reference.bim$V2)\n',
    '    if (sum(reference.dup) != 0) {\n',
    '      reference.bim <- reference.bim[!reference.dup, ]\n',
    '      reference.bed <- reference.bed[, !reference.dup]\n',
    '    }\n',
    '    for (j in 1:length(used.proteins)) {\n',
    '      load(paste0(weights.dir, used.proteins[j]))\n',
    '      protein <- gsub(".RData", "", used.proteins[j])\n',
    '      outres[outindx, 1:4] <- used.lookup[used.lookup$gene == protein, c("chr", "p.0", "p.1", "gene")]\n',
    '      outres[outindx, 5] <- R2\n',
    '      outres[outindx, 12] <- sum(weight != 0)\n',
    '      tryCatch({\n',
    '        outres[outindx, 6:11] <- TestAssociation(sumstats, weight, ldref, n.sumstats, reference.bim, reference.bed, skip.robust = TRUE)\n',
    '      }, error = function(e) {\n',
    '        cat("Warning: No overlapping SNPs found in protein:", protein, "\\n")\n',
    '        cat("Error Details:", conditionMessage(e), "\\n")\n',
    '      })\n',
    '      outindx <- outindx + 1\n',
    '      if (outindx %% 20 == 0) cat("Finished Index", outindx, "\\n")\n',
    '    }\n',
    '    cat("Finished CHR", chromosome, "\\n")\n',
    '  }\n',
    '  outres <- outres[!is.na(outres$chr), ]\n',
    '  write.table(outres, file = paste0(out.dir, "/", output.name), quote = FALSE, row.names = FALSE, sep = "\\t")\n',
    '}\n'
  )
  
  # Save script
  script_filename <- sprintf("%03d.R", i)
  writeLines(script_content, con = script_filename)
}

# To execute all generated scripts, run in the terminal:
# for script in *.R; do Rscript "$script"; done