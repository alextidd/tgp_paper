#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Went2018 (MM) and Law2017 (CLL) LD #
mkdir output/BloodCancer_Law2017andWent2018_LD/
( cat data/Went2018/Table1.tsv | 
  sed 1d | cut -f1 ;
  cat data/Law2017/41467_2017_BFncomms14175_MOESM1173_ESM_SupplementaryTable3.tsv  | 
  sed 1d | cut -f3 | awk '$0!=""'
) | cat > output/BloodCancer_Law2017andWent2018_LD/index_SNPs.txt
# run SNiPA...
SNiPA_to_variants output/BloodCancer_Law2017andWent2018_LD/

##### KNOWN GENES #####

# Dietlein2020 #
get_kgs BloodCancer \
<(  cat data/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
    awk -F'\t' -vOFS='\t' '$2 == "Blood" && $4 != "level D" {print $1}' )
