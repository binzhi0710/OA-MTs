#!/bin/bash

#SBATCH -J rscript
#SBATCH -N 2
#SBATCH -n 50
#SBATCH --mem=300G
#SBATCH -p hpc_smp,cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

module load apps/R/4.3.1

DATADIR=./garfield-data
ANNOTDIR=$DATADIR/annotation
ANNOTLINKFILE=$ANNOTDIR/link_file.txt
PTHRESH=1e-5,1e-8
BINNING=m5,n5,t5
CONDITION=0
SUBSET="1-1005"

export DATADIR ANNOTDIR ANNOTLINKFILE PTHRESH BINNING CONDITION SUBSET

process_input() {
  INPUTNAME=$1
  PRUNETAGSDIR=$DATADIR/tags/r01
  CLUMPTAGSDIR=$DATADIR/tags/r08
  MAFTSSDDIR=$DATADIR/maftssd
  PVALDIR=$DATADIR/pval/$INPUTNAME
  OUTDIR=$DATADIR/output/$INPUTNAME
  mkdir -p $OUTDIR

  F1=$OUTDIR/garfield.prep.$INPUTNAME.out
  F0=$OUTDIR/garfield.Meff.$INPUTNAME.out

  echo 'Prune and Clump for' $INPUTNAME
  echo -n > $F1
  for CHR in $(seq 1 22) # X
  do
    echo 'CHR' $CHR
    ./garfield-prep-chr -ptags $PRUNETAGSDIR/chr$CHR -ctags $CLUMPTAGSDIR/chr$CHR -maftss $MAFTSSDDIR/chr$CHR -pval $PVALDIR/chr$CHR -ann $ANNOTDIR/chr$CHR -excl 895,975,976,977,978,979,980 -chr $CHR -o $F1 || { echo 'Failure!'; }
  done

  echo 'Calculate effective number of annotations for' $INPUTNAME
  Rscript garfield-Meff-Padj.R -i $F1 -o $F0
  NEA=$(head -1 $F0 | awk '{print $2}')
  Padj=$(tail -1 $F0 | awk '{print $2}')

  echo 'Calculate Enrichment and Significance for' $INPUTNAME
  F2=$OUTDIR/garfield.test.$INPUTNAME.out
  Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -s $SUBSET -c $CONDITION
  echo 'GARFIELD single annotation analysis complete for' $INPUTNAME

  echo 'Create Plots for' $INPUTNAME
  Rscript garfield-plot.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

  echo 'GARFIELD Analysis Complete for' $INPUTNAME
}

export -f process_input

ls -d $DATADIR/pval/*/ | xargs -n 1 basename | parallel -j 10 process_input {}