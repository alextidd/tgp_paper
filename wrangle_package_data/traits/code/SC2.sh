#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Chahal2016 (BCC) and Sarin2020 (SCC) LD #
mkdir output/SC_Chahal2016andSarin2020_LD/
( cat data/Chahal2016/Table2.tsv | 
  sed 1d | cut -f1 ;
  cat data/Chahal2016/SupplementaryTable8.tsv |
  sed 1d | cut -f1 ;
  cat data/Sarin2020/SupplementaryTable6.tsv  | 
  sed 1d | awk '$3>0 {print $4}' ;
) | cat > output/SC_Chahal2016andSarin2020_LD/index_SNPs.txt
# run SNiPA...
SNiPA_to_variants output/SC_Chahal2016andSarin2020_LD/

##### KNOWN GENES #####

# Bonilla2016 and Chang2021 #

get_kgs SC2 \
<(  (
      cat data/Bonilla2016/SupplementaryTable5.tsv |
      sed 1d | cut -f1 | uniq ;
      cat data/Chang2021/Figure4B.tsv ;
    ) | cat | uniq )