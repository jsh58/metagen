#!/bin/bash -e

#SBATCH -p bos-info,general
#SBATCH -N 1
#SBATCH -n 9
#SBATCH --mem 196000
#SBATCH -t 1-00:00

module load centrifuge

idx=/n/regal/informatics_public/metagen/nt

if [ $# -lt 3 ]; then
  echo 'Usage: bash -e centrifuge.sh  <R1>  <R2>  <out>'
  echo '   OR: bash -e centrifuge.sh  <R1>  None  <out>  # for SE'
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

# classify reads
centrifuge \
  -p8 \
  -x $idx \
  $reads \
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
    > $tree
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
