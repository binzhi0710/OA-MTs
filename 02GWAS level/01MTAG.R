library(data.table)

trait_pairs <- fread("1030traitpair.txt")

base_script <- "#!/bin/bash
#SBATCH -J job_x9
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=15G
#SBATCH -p node_a,node_b,node_c
#SBATCH -t 15-00:00:00
#SBATCH --comment=runx
/usr/bin/python2.7 mtag.py \\
"

for (i in seq_len(nrow(trait_pairs))) {
  t1 <- sprintf('%02d', as.numeric(trait_pairs$trait1[i]))
  t2 <- trait_pairs$trait2[i]
  
  script_content <- paste0(
    base_script,
    "  --sumstats ./1030meta_prepare/", t1, ".txt,./OAprepare/", t2, ".txt \\\n",
    "  --out ./1030MTAGresult/", t1, "_", t2, " \\\n",
    "  --stream_stdout \\\n",
    "  --perfect_gencov \\\n",
    "  --equal_h2 \\\n",
    "  --force\n"
  )
  
  script_name <- sprintf("%03d.sh", i)
  writeLines(script_content, script_name)
}
