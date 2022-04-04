#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Sud2018 LD #
mkdir output/Lymphoma_Sud2018_LD/
cat data/Sud2018/Table2.tsv |
sed 1d | cut -f2 \
> output/Lymphoma_Sud2018_LD/index_SNPs.txt
# run SNiPA...
SNiPA_to_variants output/Lymphoma_Sud2018_LD/

##### KNOWN GENES #####

# Dietlein2020 #
get_kgs Lymphoma \
<(  cat data/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
    awk -F'\t' -vOFS='\t' '$2 == "Blood" && $4 != "level D" {print $1}' )

