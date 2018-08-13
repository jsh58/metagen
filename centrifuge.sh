#!/bin/bash -e

#SBATCH -p bos-info,general
#SBATCH -N 1
#SBATCH -n 9
#SBATCH --mem 196000
#SBATCH -t 1-00:00

module load centrifuge

if [ $# -lt 3 ]; then
  echo 'Usage: bash -e centrifuge.sh  <R1>  <R2>  <out>  [<idx>]  [<proc>]  ["mm"]'
  echo '  Required arguments:'
  echo '    <R1>      Input fastq file of R1 reads'
  echo '    <R2>      Input fastq file of R2 reads; use "None" for SE reads'
  echo '    <out>     Output html file'
  echo '  Optional arguments:'
  echo '    <idx>     Centrifuge index (def. /n/regal/informatics_public/metagen/nt)'
  echo '    <proc>    Number of processors to use with centrifuge (def. 8)'
  echo '    "mm"      Use memory-mapping option (--mm) with centrifuge'
  exit -1
fi

# determine if inputs are SE or PE
if [ "$2" == "None" ]; then
  reads="-U ${1/ /,}"
  echo 'Analyzing SE reads:'
  echo '  '$1
else
  reads="-1 ${1/ /,} -2 ${2/ /,}"
  echo 'Analyzing PE reads:'
  echo '  '$1
  echo '  '$2
fi

# save optional args
idx=/n/regal/informatics_public/metagen/nt
if [ $# -gt 3 ]; then
  idx=$4
fi
proc=8
if [[ $# -gt 4 && $5 -gt 0 ]]; then
  proc=$5
fi
mm=""
if [[ $# -gt 5 && $6 = "mm" ]]; then
  mm="--mm"
fi

# classify reads
centrifuge \
  -p $proc \
  -x $idx \
  $reads \
  $mm \
  --no-abundance \
  --report-file /dev/null \
  | centrifuge-kreport \
  -x $idx \
  --no-lca \
  > "$3.raw"

# produce taxonomy tree (if necessary)
tree=$idx.tree
if [ ! -f $tree ]; then

  centrifuge-inspect \
    --taxonomy-tree \
    $idx \
    > $tree.tmp

  # add taxonomic summary to tree
  if [[ ! -f acc2taxid.txt || ! -f $idx.fa ]]; then
    echo 'Error! Cannot add summary to taxonomy tree'
    exit -1
  fi
  python /n/regal/informatics_public/metagen/ntSumm.py \
    acc2taxid.txt \
    $tree.tmp \
    $idx.fa \
    $tree

  rm $tree.tmp
fi

# determine versions of centrifuge, nt
base=$(dirname $(which centrifuge))
if [ -f $base/VERSION ]; then
  version=$(cat $base/VERSION)
fi
base2=$(dirname $idx)
if [ -f $base2/DATE ]; then
  date=$(cat $base2/DATE)
fi

# select top 20 hits
python /n/regal/informatics_public/metagen/centSumm3.py \
  "$3.raw"  $tree  "$3"  20  $version  "$date"
echo 'Output file: '$3
