#!/bin/bash -e

#SBATCH -p general,bos-info
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 10000
#SBATCH -t 0-04:00

# download taxonomy files, check md5
rm -f taxdump.tar.gz taxdump.tar.gz.md5 nodes.dmp names.dmp merged.dmp delnodes.dmp
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
md5sum -c taxdump.tar.gz.md5
tar xzf taxdump.tar.gz  # creates many files, including 'nodes.dmp' and 'names.dmp' (required by centrifuge-build)
                        #   and 'merged.dmp' and 'delnodes.dmp' (required by updateTaxID2.py)

# download accession-to-taxID files, check md5s
rm -f nucl_gb.accession2taxid.gz nucl_gb.accession2taxid.gz.md5
wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5
md5sum -c nucl_gb.accession2taxid.gz.md5

rm -f nucl_wgs.accession2taxid.gz nucl_wgs.accession2taxid.gz.md5
wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5
md5sum -c nucl_wgs.accession2taxid.gz.md5

#rm -f nucl_est.accession2taxid.gz nucl_est.accession2taxid.gz.md5
#wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz
#wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz.md5
#md5sum -c nucl_est.accession2taxid.gz.md5

#rm -f nucl_gss.accession2taxid.gz nucl_gss.accession2taxid.gz.md5
#wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz
#wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz.md5
#md5sum -c nucl_gss.accession2taxid.gz.md5

# update merged and deleted taxIDs from accession2taxid files, and then combine them
python updateTaxID2.py \
  merged.dmp delnodes.dmp \
  acc2taxid.txt \
  nucl_gb.accession2taxid.gz \
  nucl_wgs.accession2taxid.gz
#  nucl_est.accession2taxid.gz \
#  nucl_gss.accession2taxid.gz
