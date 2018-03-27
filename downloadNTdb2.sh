#!/bin/bash -e

#SBATCH -p general,bos-info
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 24000
#SBATCH -t 3-00:00

module load blast
module load bedtools2

# usage (two options):
# - download nt; filter; mask regions in adapters2.bed:
#   $ bash downloadNTdb2.sh
# - download nt; filter; blast given <fasta> against nt;
#   create bed regions of matches; mask those regions:
#   $ bash downloadNTdb2.sh  <fasta>

# download nt database, check md5
rm -f nt.gz nt.gz.md5
date +'%F %r' > DATE  # record date/time of download
wget -q ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
wget -q ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz.md5
md5sum -c nt.gz.md5

# filter out sequences that are pure Ns or shorter than 30bp;
#   mask subsequences that match adapters (adapters2.bed)
if [[ $# -eq 0 && -f adapters2.bed ]]; then
  mask="adapters2.bed"
fi
python filterNT2.py nt.gz nt.fa 30 $mask

# search for matching adapter sequences
if [[ $# -gt 0 && -f $1 ]]; then
  rm -f ${1%.*}.blast ${1%.*}.bed

  # create blast db
  rm -f nt.nal nt.*.nhr nt.*.nin nt.*.nsq
  makeblastdb \
    -dbtype nucl \
    -in nt.fa \
    -out nt \
    -logfile makeblastdb.log

  # megablast
  blastn \
    -query $1 \
    -db nt \
    -word_size 18 -ungapped \
    -evalue 250 \
    -max_target_seqs 10000 \
    -out ${1%.*}.blast \
    -outfmt 6

  # convert output to BED
  awk 'OFS="\t" {if ($9 < $10) print $2, $9-1, $10; else print $2, $10-1, $9}' \
    ${1%.*}.blast \
    | sort -k 1,1 -k 2,2n - \
    | bedtools merge -i - > ${1%.*}.bed

  # mask sequences
  python filterNT2.py nt.fa nt.fa.tmp 30 ${1%.*}.bed
  mv nt.fa.tmp nt.fa

fi
