#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Sarin2020 (SCC) LD #
mkdir output/SCC_Sarin2020_LD/
( cat data/Sarin2020/SupplementaryTable6.tsv  | 
  sed 1d | awk '$3>0 {print $4}' ;
) | cat > output/SCC_Sarin2020_LD/index_SNPs.txt
# run SNiPA...
SNiPA_to_variants output/SCC_Sarin2020_LD/

##### KNOWN GENES #####

# Chang2021 #
get_kgs SCC \
<(  cat data/Chang2021/Figure4B.tsv )
