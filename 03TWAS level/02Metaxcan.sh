#!/bin/bash

GWAS_TOOLS=~/S-metaxcan/MetaXcan-master/software/summary-gwas-imputation-master/src
DATA=~/S-metaxcan/data
GWAS=~/S-MTAG/MTAGresult
OUTPUT=~/S-metaxcan
METAXCAN=~/S-metaxcan/MetaXcan-master/software

# Create scripts directory
mkdir -p "$OUTPUT/scripts"

# Process each GWAS file
for gwas_file in "$GWAS"/*.txt; do
    base_name=$(basename "$gwas_file" .txt)
    script_file="$OUTPUT/scripts/${base_name}_analysis.sh"

    # Write analysis script
    cat << EOF > "$script_file"
#!/bin/bash

#SBATCH -J metaxcan_${base_name}
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=15G
#SBATCH -p cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

conda activate MetaXcan

# GWAS Harmonization
python "$GWAS_TOOLS/gwas_parsing.py" \\
    -gwas_file "$GWAS/${base_name}.txt" \\
    -liftover "$DATA/liftover/hg19ToHg38.over.chain.gz" \\
    -snp_reference_metadata "$DATA/reference_panel_1000G/variant_metadata.txt.gz" METADATA \\
    -output_column_map SNP variant_id \\
    -output_column_map A2 non_effect_allele \\
    -output_column_map A1 effect_allele \\
    -output_column_map mtag_beta effect_size \\
    -output_column_map mtag_pval pvalue \\
    -output_column_map CHR chromosome \\
    --chromosome_format \\
    -output_column_map BP position \\
    -output_column_map meta_freq frequency \\
    --insert_value sample_size 184305 \\
    -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size \\
    -output "$OUTPUT/harmonized_gwas/${base_name}.txt.gz"

# GWAS Imputation (22 chromosomes × 10 batches)
for chr in {1..22}; do
    for batch in {0..9}; do
        python "$GWAS_TOOLS/gwas_summary_imputation.py" \\
        -by_region_file "$DATA/eur_ld.bed.gz" \\
        -gwas_file "$OUTPUT/harmonized_gwas/${base_name}.txt.gz" \\
        -parquet_genotype "$DATA/reference_panel_1000G/chr\${chr}.variants.parquet" \\
        -parquet_genotype_metadata "$DATA/reference_panel_1000G/variant_metadata.parquet" \\
        -window 100000 \\
        -parsimony 7 \\
        -chromosome \${chr} \\
        -regularization 0.1 \\
        -frequency_filter 0.01 \\
        -sub_batches 10 \\
        -sub_batch \${batch} \\
        --standardise_dosages \\
        -output "$OUTPUT/summary_imputation/${base_name}_chr\${chr}_sb\${batch}_reg0.1_ff0.01_by_region.txt.gz"
    done
done

# Combine Imputation Results
python "$GWAS_TOOLS/gwas_summary_imputation_postprocess.py" \\
    -gwas_file "$OUTPUT/harmonized_gwas/${base_name}.txt.gz" \\
    -folder "$OUTPUT/summary_imputation" \\
    -pattern "${base_name}.*" \\
    -parsimony 7 \\
    -output "$OUTPUT/processed_summary_imputation/${base_name}.txt.gz"

# PrediXcan Analysis (across 49 tissues)
for db_file in \$DATA/models/eqtl/mashr/*.db; do
    tissue_name=\$(basename "\$db_file" .db)
    cov_file="\$DATA/models/eqtl/mashr/\${tissue_name}.txt.gz"
    python "$METAXCAN/SPrediXcan.py" \\
        --gwas_file "$OUTPUT/processed_summary_imputation/${base_name}.txt.gz" \\
        --snp_column panel_variant_id \\
        --effect_allele_column effect_allele \\
        --non_effect_allele_column non_effect_allele \\
        --zscore_column zscore \\
        --model_db_path "\$db_file" \\
        --covariance "\$cov_file" \\
        --keep_non_rsid \\
        --additional_output \\
        --model_db_snp_key varID \\
        --throw \\
        --output_file "$OUTPUT/spredixcan/eqtl/${base_name}__PM__\${tissue_name}.csv"
    echo "Finished processing tissue \${tissue_name} for ${base_name}"
done

# MultiXcan Analysis
python "$METAXCAN/SMulTiXcan.py" \\
    --models_folder "$DATA/models/eqtl/mashr" \\
    --models_name_pattern "mashr_(.*).db" \\
    --snp_covariance "$DATA/models/gtex_v8_expression_mashr_snp_covariance.txt.gz" \\
    --metaxcan_folder "$OUTPUT/spredixcan/eqtl/" \\
    --metaxcan_filter "${base_name}_(.*).csv" \\
    --metaxcan_file_name_parse_pattern "(.*)__PM__(.*).csv" \\
    --gwas_file "$OUTPUT/processed_summary_imputation/${base_name}.txt.gz" \\
    --snp_column panel_variant_id \\
    --effect_allele_column effect_allele \\
    --non_effect_allele_column non_effect_allele \\
    --zscore_column zscore \\
    --keep_non_rsid \\
    --model_db_snp_key varID \\
    --cutoff_condition_number 30 \\
    --verbosity 7 \\
    --throw \\
    --output "$OUTPUT/smultixcan/eqtl/${base_name}_smultixcan.txt"

echo "Finished processing ${base_name}"
EOF

    # Make script executable
    chmod +x "$script_file"
    echo "Generated script for $base_name"
done

# To submit all generated scripts, run:
# for script in ~/S-metaxcan/scripts/*_analysis.sh; do sbatch "$script"; done