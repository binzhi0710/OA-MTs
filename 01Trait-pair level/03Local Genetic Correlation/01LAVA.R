# Read trait pairs file
gwas_pairs <- read.table("LAVAGWAS.txt", header = TRUE)

# Ensure "R_scripts" folder exists in working directory
if (!dir.exists("R_scripts")) {
  dir.create("R_scripts")
}

# Generate 152 R scripts
for (i in 1:152) {
  
  # Extract trait pair for the i-th row
  phenos <- c(gwas_pairs$GWAS1[i], gwas_pairs$GWAS2[i])
  
  # Script content
  script_content <- paste0(
    'library(LAVA)\n\n',
    '# Read loci file\n',
    'loci <- read.loci("LAVA.loci")\n',
    'n.loc = nrow(loci)\n\n',
    '# Extract trait pair for row ', i, '\n',
    'phenos <- c("', phenos[1], '", "', phenos[2], '")\n\n',
    '# Read input information file\n',
    'input <- process.input(input.info.file = "lavainput.txt", \n',
    '                       sample.overlap.file = "sample.overlap.txt",\n',
    '                       ref.prefix = "g1000_eur", \n',
    '                       phenos = phenos)\n\n',
    'progress = ceiling(quantile(1:n.loc, seq(.0005,1,.0005)))\n\n',
    'univ.p.thresh = 0.05\n\n',
    'u=b=list()\n',
    'for (i in 1:n.loc) {\n',
    '  if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))    \n',
    '  locus = process.locus(loci[i,], input)                                          \n',
    '  if (!is.null(locus)) {\n',
    '    loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)\n',
    '    loc.out = run.univ.bivar(locus, univ.thresh = univ.p.thresh)\n',
    '    u[[i]] = cbind(loc.info, loc.out$univ)\n',
    '    if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)\n',
    '  }\n',
    '}\n\n',
    '# Specify output file name and path\n',
    'out.fname = paste0("LAVAresult/", phenos[1], "_", phenos[2])\n\n',
    '# Save output results to specified folder\n',
    'write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T,sep = "\\t")\n',
    'write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T,sep = "\\t")\n\n',
    'print(paste0("Done! Analysis output written to LAVAresult folder as ",out.fname,".*.lava"))\n'
  )
  
  # Write script file with three-digit naming format
  script_filename <- sprintf("R_scripts/%03d.R", i)
  writeLines(script_content, con = script_filename)
}

print("All R scripts have been generated.")

# Ensure "sh_scripts" folder exists in working directory
if (!dir.exists("sh_scripts")) {
  dir.create("sh_scripts")
}

# Generate 152 .sh scripts
for (i in 1:152) {
  
  # Generate script content
  sh_content <- paste0(
    '#!/bin/bash\n\n',
    '#SBATCH -J rscript\n',
    '#SBATCH -N 1\n',
    '#SBATCH -n 4\n',
    '#SBATCH --mem=15G\n',
    '#SBATCH -p cydell_1,cydell_2,cydell_3,cyhpc_1\n',
    '#SBATCH -t 20-12:00:00\n',
    '#SBATCH --comment=test\n\n',
    'module load apps/R/4.3.1\n\n',
    'Rscript ./R_scripts/', sprintf("%03d.R", i), '\n'
  )
  
  # Write script file with three-digit naming format
  sh_filename <- sprintf("sh_scripts/%03d.sh", i)
  writeLines(sh_content, con = sh_filename)
}

print("All .sh scripts have been generated.")

# Submit all .sh scripts to the cluster using sbatch
for (i in 1:152) {
  sh_filename <- sprintf("sh_scripts/%03d.sh", i)
  system(paste("sbatch", sh_filename))
  cat("Submitted:", sh_filename, "\n")
}

print("All .sh scripts have been submitted to the cluster.")