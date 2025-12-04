library(dplyr)
library(stringr)
library(data.table)
library(future)
setDTthreads(32)
source("./CPASSOC_FunctionSet.R")
snplist <- fread(input = "plink_pruning.prune.in", header = FALSE, data.table = FALSE)
trait_pairs <- fread("1030traitpair.txt")

out_path <- "./data_cpassoc"
dir.create(out_path, showWarnings = FALSE)
out_res_path <- "./result_cpassoc"
dir.create(out_res_path, showWarnings = FALSE)

for (i in 1:nrow(trait_pairs)) {
  
  t1 <- sprintf("%02d", as.numeric(trait_pairs$trait1[i]))
  t2 <- trait_pairs$trait2[i]
  
  p1 <- file.path("./00metabo/dataprepared0.01maf/", paste0(t1, ".txt"))
  p2 <- file.path("./00OA_noMHC/maf0.01/", paste0(t2, ".txt"))
  
  d1 <- fread(p1) %>%
    rename(SNP=SNP,effect_allele=A1,other_allele=A2,beta=beta,se=se,eaf=eaf,pval=pval,samplesize=N,z=Z) %>%
    mutate(chr=as.numeric(chr),pos=as.numeric(pos)) %>%
    na.omit()
  
  d2 <- fread(p2) %>%
    rename(SNP=SNP,effect_allele=A1,other_allele=A2,beta=beta,se=se,eaf=eaf,pval=pval,samplesize=N,z=Z) %>%
    mutate(chr=as.numeric(chr),pos=as.numeric(pos),effect_allele=toupper(effect_allele),other_allele=toupper(other_allele)) %>%
    na.omit()
  
  f1 <- paste0(out_path, "/", t1, "_formatted.txt")
  f2 <- paste0(out_path, "/", t2, "_formatted.txt")
  
  fwrite(d1, file=f1, sep="\t", quote=FALSE, row.names=FALSE)
  fwrite(d2, file=f2, sep="\t", quote=FALSE, row.names=FALSE)
  
  d <- cpassoc.combine_assoc(files_path=c(f1,f2), traits=c(t1,t2))
  
  fout1 <- paste0(out_path, "/", t1, "_", t2, ".txt")
  fout2 <- paste0(out_res_path, "/", t1, "_", t2, "_CPASSOC_SHet_SHom.txt")
  
  fwrite(d, file=fout1, quote=FALSE, sep="\t", row.names=FALSE)
  
  cm <- d %>%
    dplyr::filter(SNP %in% snplist$V1) %>%
    dplyr::filter(dplyr::if_all(tidyselect::starts_with("z"), ~ abs(.x) < 1.96)) %>%
    dplyr::select(tidyselect::starts_with("z")) %>%
    as.matrix() %>%
    cor()
  
  options(future.globals.maxSize = 4 * 1024^3)
  
  out <- cpassoc.shet_shom_test(
    gwas_data = d,
    SHet_test = TRUE,
    SHom_test = TRUE,
    workers = 24,
    mvrnorm_n = 1e6,
    set_seed = 666
  )
  
  fwrite(out, file=fout2, quote=FALSE, sep="\t", row.names=FALSE)
}
