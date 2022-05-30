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

# Dietlein2020 #
get_kgs SC \
<(  cat data/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
    awk -F'\t' -vOFS='\t' '$2 == "Skin" && $4 != "level D" {print $1}' )
