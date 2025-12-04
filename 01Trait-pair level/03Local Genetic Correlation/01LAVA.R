gwas_pairs <- read.table("LAVAGWAS.txt", header = TRUE)

if (!dir.exists("R_scripts")) dir.create("R_scripts")

for (i in 1:152) {
  
  phenos <- c(gwas_pairs$GWAS1[i], gwas_pairs$GWAS2[i])
  
  script_content <- paste0(
    'library(LAVA)\n',
    'loci <- read.loci("LAVA.loci")\n',
    'n.loc = nrow(loci)\n',
    'phenos <- c("', phenos[1], '", "', phenos[2], '")\n',
    'input <- process.input(input.info.file = "lavainput.txt", ',
    'sample.overlap.file = "sample.overlap.txt", ',
    'ref.prefix = "g1000_eur", phenos = phenos)\n',
    'progress = ceiling(quantile(1:n.loc, seq(.0005,1,.0005)))\n',
    'univ.p.thresh = 0.05\n',
    'u=b=list()\n',
    'for (i in 1:n.loc) {\n',
    '  if (i %in% progress) {}\n',
    '  locus = process.locus(loci[i,], input)\n',
    '  if (!is.null(locus)) {\n',
    '    loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)\n',
    '    loc.out = run.univ.bivar(locus, univ.thresh = univ.p.thresh)\n',
    '    u[[i]] = cbind(loc.info, loc.out$univ)\n',
    '    if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)\n',
    '  }\n',
    '}\n',
    'out.fname = paste0("LAVAresult/", phenos[1], "_", phenos[2])\n',
    'write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T,sep="\\t")\n',
    'write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T,sep="\\t")\n'
  )
  
  script_filename <- sprintf("R_scripts/%03d.R", i)
  writeLines(script_content, con = script_filename)
}

if (!dir.exists("sh_scripts")) dir.create("sh_scripts")

for (i in 1:152) {
  
  sh_content <- paste0(
    '#!/bin/bash\n',
    '#SBATCH -J jobx_', i, '\n',
    '#SBATCH -N 1\n',
    '#SBATCH -n 4\n',
    '#SBATCH --mem=15G\n',
    '#SBATCH -p node_a,node_b,node_c\n',
    '#SBATCH -t 15-00:00:00\n',
    '#SBATCH --comment=runx\n',
    'module load apps/R/4.3.1\n',
    'Rscript ./R_scripts/', sprintf("%03d.R", i), '\n'
  )
  
  sh_filename <- sprintf("sh_scripts/%03d.sh", i)
  writeLines(sh_content, con = sh_filename)
}

for (i in 1:152) {
  sh_filename <- sprintf("sh_scripts/%03d.sh", i)
  system(paste("sbatch", sh_filename))
}
