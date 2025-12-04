library(dplyr)
library(coloc)
library(data.table)
library(readxl)

setDTthreads(32)

fuma_data <- read_excel("MTAGcolocloci.xlsx")

results <- data.frame(
  PSY=character(),
  osteoarthritis=character(),
  chr=numeric(),
  start=numeric(),
  end=numeric(),
  df1_file=character(),
  df2_file=character(),
  PP.H0.abf=numeric(),
  PP.H1.abf=numeric(),
  PP.H2.abf=numeric(),
  PP.H3.abf=numeric(),
  PP.H4.abf=numeric(),
  stringsAsFactors=FALSE
)

for (i in 1:nrow(fuma_data)) {
  x1 <- fuma_data$PSY[i]
  x2 <- fuma_data$osteoarthritis[i]
  c1 <- as.numeric(fuma_data$chr[i])
  s1 <- as.numeric(fuma_data$start[i])
  e1 <- as.numeric(fuma_data$end[i])
  
  f1 <- file.path("00PSYGWAS_noMHC", paste0(x1, ".txt"))
  f2 <- file.path("00OAGWAS_noMHC", paste0(x2, ".txt"))
  
  d1 <- fread(f1)
  d2 <- fread(f2)
  
  d1 <- data.frame(SNP=d1$SNP, chrom=d1$chr, pos=d1$pos, A1=d1$A1, A2=d1$A2, beta=d1$beta, se=d1$se, eaf=d1$eaf)
  d1$MAF <- ifelse(d1$eaf < 0.5, d1$eaf, 1 - d1$eaf)
  d1 <- subset(d1, chrom == c1 & pos >= s1 & pos <= e1)
  
  d2 <- data.frame(SNP=d2$SNP, A1=d2$A1, A2=d2$A2, beta=d2$beta, se=d2$se)
  
  m <- merge(d1, d2, by="SNP")
  m <- m[!duplicated(m$SNP),]
  m$A1.y <- toupper(m$A1.y)
  m$A2.y <- toupper(m$A2.y)
  
  m <- m %>%
    filter((A1.x==A1.y & A2.x==A2.y) | (A1.x==A2.y & A2.x==A1.y)) %>%
    mutate(beta.y = ifelse(A1.x==A1.y, beta.y, -beta.y))
  
  m$VAR1 <- m$se.x^2
  m$VAR2 <- m$se.y^2
  m <- m[m$VAR1 != 0 & m$VAR2 != 0,]
  m <- na.omit(m)
  
  c1d <- data.frame(beta=m$beta.x, varbeta=m$VAR1, snp=m$SNP, MAF=m$MAF)
  c2d <- data.frame(beta=m$beta.y, varbeta=m$VAR2, snp=m$SNP)
  c1d <- as.list(c1d)
  c2d <- as.list(c2d)
  
  c1d$type <- "cc"
  c2d$type <- "cc"
  
  r <- coloc.abf(c1d, c2d, p1=1e-4, p2=1e-4, p12=1e-5)
  
  results <- rbind(
    results,
    data.frame(
      PSY=x1,
      osteoarthritis=x2,
      chr=c1,
      start=s1,
      end=e1,
      df1_file=f1,
      df2_file=f2,
      PP.H0.abf=r$summary[["PP.H0.abf"]],
      PP.H1.abf=r$summary[["PP.H1.abf"]],
      PP.H2.abf=r$summary[["PP.H2.abf"]],
      PP.H3.abf=r$summary[["PP.H3.abf"]],
      PP.H4.abf=r$summary[["PP.H4.abf"]]
    )
  )
}

write.csv(results, "MTAGlocicoloc_results.csv", row.names=FALSE)
