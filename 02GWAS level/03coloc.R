# Load required libraries
library(dplyr)
library(coloc)
library(data.table)
library(readxl)

setDTthreads(32)

# Read MTAGcolocloci.xlsx file
fuma_data <- read_excel("MTAGcolocloci.xlsx")

# Create an empty data frame to store results
results <- data.frame(
  PSY = character(),
  osteoarthritis = character(),
  chr = numeric(),
  start = numeric(),
  end = numeric(),
  df1_file = character(),
  df2_file = character(),
  PP.H0.abf = numeric(),
  PP.H1.abf = numeric(),
  PP.H2.abf = numeric(),
  PP.H3.abf = numeric(),
  PP.H4.abf = numeric(),
  stringsAsFactors = FALSE
)

# Iterate through each row to read corresponding files and perform analysis
for (i in 1:nrow(fuma_data)) {
  PSY_id <- fuma_data$PSY[i]
  osteoarthritis_id <- fuma_data$osteoarthritis[i]
  
  # Get chr, start, and end for the current row
  chr <- as.numeric(fuma_data$chr[i])
  start <- as.numeric(fuma_data$start[i])
  end <- as.numeric(fuma_data$end[i])
  
  # Construct file paths
  PSY_file <- file.path("00PSYGWAS_noMHC", paste0(PSY_id, ".txt"))
  osteoarthritis_file <- file.path("00OAGWAS_noMHC", paste0(osteoarthritis_id, ".txt"))
  
  # Read files
  df1 <- fread(PSY_file, header = TRUE)
  df2 <- fread(osteoarthritis_file, header = TRUE)
  
  # Preprocess df1
  df1 <- data.frame(SNP = df1$SNP, chrom = df1$chr, pos = df1$pos, A1 = df1$A1, A2 = df1$A2, beta = df1$beta, se = df1$se, eaf = df1$eaf)
  df1$MAF <- ifelse(df1$eaf < 0.5, df1$eaf, 1 - df1$eaf)
  
  # Filter region for colocalization
  df1 <- subset(df1, chrom == chr & pos >= start & pos <= end)
  
  # Preprocess df2
  df2 <- data.frame(SNP = df2$SNP, A1 = df2$A1, A2 = df2$A2, beta = df2$beta, se = df2$se)
  
  # Merge data frames
  dfall <- merge(df1, df2, by = "SNP")
  dfall <- dfall[!duplicated(dfall$SNP), ]
  
  # Convert A1.y and A2.y to uppercase
  dfall$A1.y <- toupper(dfall$A1.y)
  dfall$A2.y <- toupper(dfall$A2.y)
  
  # Filter and adjust beta.y
  dfall <- dfall %>%
    filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y)) %>%
    mutate(beta.y = ifelse(A1.x == A1.y, beta.y, -beta.y))
  
  dfall$VAR1 <- dfall$se.x^2
  dfall$VAR2 <- dfall$se.y^2
  dfall <- dfall[dfall$VAR1 != 0 & dfall$VAR2 != 0, ]
  
  # Remove rows with NA values
  dfall <- na.omit(dfall)
  
  cdf1 <- data.frame(beta = dfall$beta.x, varbeta = dfall$VAR1, snp = dfall$SNP, MAF = dfall$MAF)
  cdf2 <- data.frame(beta = dfall$beta.y, varbeta = dfall$VAR2, snp = dfall$SNP)
  
  cdf1 <- as.list(cdf1)
  cdf2 <- as.list(cdf2)
  
  cdf1$type <- "cc"  # Continuous variable
  cdf2$type <- "cc"  # Binary variable
  
  colocresult <- coloc.abf(cdf1, cdf2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  
  # Extract colocalization results and store
  result_row <- data.frame(
    PSY = PSY_id,
    osteoarthritis = osteoarthritis_id,
    chr = chr,
    start = start,
    end = end,
    df1_file = PSY_file,
    df2_file = osteoarthritis_file,
    PP.H0.abf = colocresult$summary[["PP.H0.abf"]],
    PP.H1.abf = colocresult$summary[["PP.H1.abf"]],
    PP.H2.abf = colocresult$summary[["PP.H2.abf"]],
    PP.H3.abf = colocresult$summary[["PP.H3.abf"]],
    PP.H4.abf = colocresult$summary[["PP.H4.abf"]]
  )
  
  results <- rbind(results, result_row)
}

# Save results data frame to CSV file
write.csv(results, "MTAGlocicoloc_results.csv", row.names = FALSE)