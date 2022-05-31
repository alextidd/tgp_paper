#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Chahal2016 (BCC) LD #
mkdir output/BCC_Chahal2016_LD/
( cat data/Chahal2016/Table2.tsv | 
  sed 1d | cut -f1 ;
  cat data/Chahal2016/SupplementaryTable8.tsv |
  sed 1d | cut -f1 ;
) | cat > output/BCC_Chahal2016_LD/index_SNPs.txt
# run SNiPA...
SNiPA_to_variants output/BCC_Chahal2016_LD/

# GWASC assoc LD #
mkdir output/BCC_GWASC_assoc/
cat data/GWASCatalog/gwas_catalog_v1.0-associations_e105_r2022-03-08.tsv |
awk -F'\t' -vOFS='\'t '$35=="basal cell carcinoma" {print substr($21, 1, length($21)-2)}' \
> output/BCC_GWASC_assoc/index_SNPs.txt
SNiPA_to_variants output/BCC_GWASC_assoc/

##### KNOWN GENES #####

# Bonilla2016 #
get_kgs BCC \
<(  cat data/Bonilla2016/SupplementaryTable5.tsv |
    sed 1d | cut -f1 | uniq )


