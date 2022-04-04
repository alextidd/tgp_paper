#!/bin/bash
# run on hpcapp01 head node: #
cd "$(dirname "$0")"/.. ; pwd

GWASC_metadata=data/GWASCatalog/gwas-catalog-v1.0.3-studies-r2022-03-08.tsv
FTPPath=ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
roundup() { awk -v n=$1 -v d=$2 'BEGIN{print int((n+d-1/2)/d) * d}'; }

# get GWASCatalog associations
wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative \
  -O data/GWASCatalog/gwas_catalog_v1.0-associations_e105_r2022-03-08.tsv

# get GWASCatalog studies info
wget https://www.ebi.ac.uk/gwas/api/search/downloads/studies_new \
  -O $GWASC_metadata

# get summary stats for studies of interest
cat output/metadata.tsv |
awk -F'\t' '$6 == "summstats"' |
cut -f1-5 |
while IFS=$'\t' read -r trait trait_info PMID first_author year ; do
  Trait_info=${trait_info^} ; echo $Trait_info
  GCST=$(awk -F'\t' -v OFS='\t' \
          -v PMID="$PMID" -v Trait_info="$Trait_info" \
          '$2 == PMID && $8 == Trait_info {print $15}' \
          $GWASC_metadata)
  GCST_n=$(echo ${GCST#GCST} | sed -e 's/^[0]*//')
  GCST_upper=$(roundup $GCST_n 1000})
  GCST_lower="$(($GCST_upper-999))"
  GCST_bin=$(printf 'GCST%06d-GCST%06d' $GCST_lower $GCST_upper)

  url="$FTPPath/$GCST_bin/$GCST/"
  summstatsDir=data/traits/$trait/variants/${first_author}${year}_summstats/
  echo $url > $summstatsDir/README.txt
  wget -r -nH --no-parent --cut-dirs 6 -P $summstatsDir $url/*
done

