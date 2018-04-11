#!/bin/bash -e

#SBATCH -p bigmem,bos-info
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 500000
#SBATCH -t 3-00:00

module load centrifuge

# append adapters to nt
cat adapters2.fa >> nt.fa
cat nodes2.dmp >> nodes.dmp
cat names2.dmp >> names.dmp
cat acc2taxid2.txt >> acc2taxid.txt

# index database
db=nt.fa
pre=nt
rm -f $pre.1.cf $pre.2.cf $pre.3.cf $pre.4.cf $pre.tree
centrifuge-build \
  -p4 \
  --conversion-table acc2taxid.txt \
  --taxonomy-tree nodes.dmp \
  --name-table names.dmp \
  $db \
  $pre \
  &> log2.txt

# produce taxonomy tree
centrifuge-inspect \
  --taxonomy-tree \
  $pre \
  > $pre.tree
