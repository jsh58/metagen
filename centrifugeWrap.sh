#!/bin/bash -e

#SBATCH -p general,bos-info
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem 196000
#SBATCH -t 2-00:00

if [ $# -lt 1 ]; then
  echo 'Usage: bash -e centrifugeWrap.sh  <fol>  [<lane>]'
  echo '  <fol>   Folder of sequencing run'
  echo '  <lane>  Lane to analyze (def. Lane1)'
  exit -1
fi

fol=$1
mkdir -p $fol
lane="Lane1"
if [ $# -gt 1 ]; then
  lane=$2
fi
lanes=( Lane1 Lane2 Lane3 Lane4 Lane5 Lane6 Lane7 Lane8 )

for f in /n/seqcfs/sequencing/analysis_finished/$fol/$lane.*/Fastq/*R1.fastq.gz; do

  # skip undetermined
  if [[ $f =~ Undetermined ]]; then
    continue
  fi

  # check for R2
  f2=None
  r2=${f%.R1.fastq.gz}.R2.fastq.gz
  if [ -f $r2 ]; then
    f2=$r2
  fi

  # skip if output already exists
  dir=$(dirname $f)
  base=$(basename ${f%.R1.fastq.gz})
  if [ -f $fol/$base.html ]; then
    continue
  fi
  echo $base

  # check for files in other lanes
  for ln in ${lanes[@]}; do

    if [ $ln == $lane ]; then
      continue
    fi

    r1=${dir/$lane/$ln}/$base.R1.fastq.gz
    r2=${dir/$lane/$ln}/$base.R2.fastq.gz
    if [ -f $r1 ]; then
      f=$f','$r1
    fi
    if [ -f $r2 ]; then
      f2=$f2','$r2
    fi

  done

  # run centrifuge
  bash -e centrifuge.sh $f $f2 $fol/$base.html

done
