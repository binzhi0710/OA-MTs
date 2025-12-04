library(readxl)
library(dplyr)
library(data.table)

loci_data <- read_excel("1030metabo_GWAS-PW_loci.xlsx")
loci_data$phen1 <- sprintf("%02d", loci_data$phen1)
grouped_sources <- split(loci_data, loci_data$source)

dir.create("R_scripts1030", showWarnings = FALSE)

for (source_name in names(grouped_sources)) {
  
  source_data <- grouped_sources[[source_name]]
  
  script_content <- paste0(
    "library(dplyr)\n",
    "library(data.table)\n",
    "loci <- data.frame(chr = c(", paste(source_data$chr, collapse = ", "), "), start = c(", paste(source_data$start, collapse = ", "), "), stop = c(", paste(source_data$stop, collapse = ", "), "), phen1 = '", source_data$phen1[1], "', phen2 = '", source_data$phen2[1], "')\n",
    "f1 <- paste0('~/07metabo_OA/1030remake/', loci$phen1[1], '.txt')\n",
    "d1 <- fread(f1)\n",
    "d1 <- d1 %>% select(SNP, chr, pos, Z, se)\n",
    "colnames(d1) <- c('SNP','CHR','POS',paste0('Z_',loci$phen1[1]),paste0('V_',loci$phen1[1]))\n",
    "d1$CHR <- paste0('chr', d1$CHR)\n",
    "f2 <- paste0('~/07metabo_OA/00OA_noMHC/', loci$phen2[1], '.txt')\n",
    "d2 <- fread(f2)\n",
    "d2 <- d2 %>% select(SNP, chr, pos, Z, se)\n",
    "colnames(d2) <- c('SNP','CHR','POS',paste0('Z_',loci$phen2[1]),paste0('V_',loci$phen2[1]))\n",
    "d2$CHR <- paste0('chr', d2$CHR)\n",
    "x1 <- data.frame()\n",
    "x2 <- data.frame()\n",
    "for (i in 1:nrow(loci)) {\n",
    "  c <- loci$chr[i]\n",
    "  s <- loci$start[i]\n",
    "  e <- loci$stop[i]\n",
    "  y1 <- d1 %>% filter(CHR == paste0('chr',c), POS >= s, POS <= e)\n",
    "  y2 <- d2 %>% filter(CHR == paste0('chr',c), POS >= s, POS <= e)\n",
    "  x1 <- rbind(x1,y1)\n",
    "  x2 <- rbind(x2,y2)\n",
    "}\n",
    "m <- merge(x1,x2,by='SNP',all=TRUE)\n",
    "m <- distinct(m,SNP,.keep_all=TRUE)\n",
    "m <- na.omit(m)\n",
    "m <- dplyr::select(m,SNP,CHR.x,POS.x,paste0('Z_',loci$phen1[1]),paste0('V_',loci$phen1[1]),paste0('Z_',loci$phen2[1]),paste0('V_',loci$phen2[1]))\n",
    "colnames(m) <- c('SNPID','CHR','POS',paste0('Z_',loci$phen1[1]),paste0('V_',loci$phen1[1]),paste0('Z_',loci$phen2[1]),paste0('V_',loci$phen2[1]))\n",
    "for (i in 1:nrow(loci)) {\n",
    "  c <- loci$chr[i]\n",
    "  s <- loci$start[i]\n",
    "  e <- loci$stop[i]\n",
    "  m <- m[!(m$CHR == paste0('chr',c) & (m$POS == s | m$POS == e)),]\n",
    "}\n",
    "m$CHR_numeric <- as.numeric(gsub('chr','',m$CHR))\n",
    "m <- m[order(m$CHR_numeric,m$POS),]\n",
    "m <- distinct(m,CHR,POS,.keep_all=TRUE)\n",
    "m$CHR_numeric <- NULL\n",
    "out1 <- paste0('GWASfile1030/','", source_name, "','.txt.gz')\n",
    "fwrite(m,out1,compress='gzip',sep='\\t')\n",
    "b <- loci %>% select(chr,start,stop)\n",
    "b$chr <- paste0('chr',b$chr)\n",
    "b$CHR_numeric <- as.numeric(gsub('chr','',b$chr))\n",
    "b <- b[order(b$CHR_numeric,b$start),]\n",
    "b$CHR_numeric <- NULL\n",
    "out2 <- paste0('bedfile1030/','", source_name, "','.bed')\n",
    "fwrite(b,out2,col.names=FALSE,sep='\\t')\n"
  )
  
  script_file <- paste0("R_scripts1030/", source_name, "_script.R")
  writeLines(script_content, script_file)
}
