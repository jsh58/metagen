#!/bin/bash -e

#SBATCH -p bigmem,bos-info
#SBATCH -N 1
#SBATCH -n 9
#SBATCH --mem 392000
#SBATCH -t 1-00:00

module load centrifuge

cent=centrifuge
cent_kr=centrifuge-kreport
cent_insp=centrifuge-inspect

idx=nt

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
$cent \
  -p8 \
  -x $idx \
  $reads \
  --report-file /dev/null \
  | $cent_kr \
  -x $idx \
  > "$3.raw"

# produce taxonomy tree (if necessary)
tree=$idx.tree
if [ ! -f $tree ]; then
  $cent_insp \
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
python centSumm2.py \
  "$3.raw"  $tree  "$3"  20  $version  "$date"
echo 'Output file: '$3
