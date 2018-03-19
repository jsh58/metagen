#!/bin/bash -e

#SBATCH -p bigmem,bos-info
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem 500000
#SBATCH -t 3-00:00

module load centrifuge

# index database
db=nt.fa
pre=nt
rm -f $pre.1.cf $pre.2.cf $pre.3.cf $pre.4.cf
centrifuge-build \
  -p16 \
  --conversion-table acc2taxid.txt \
  --taxonomy-tree nodes.dmp \
  --name-table names.dmp \
  $db \
  $pre \
  &> log2.txt
