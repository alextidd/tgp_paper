#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# McKay2017_LD #
mkdir output/LC_McKay2017_LD/
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  # LD variants
  cat data/McKay2017/SupplementaryTable7.tsv | 
  awk -F'\t' -vOFS='\t' 'NR > 1 {print "chr"$2,$3-1,$3,$1,$4}' |
  sort -k1,1 -k2,2n ) |
cat > output/LC_McKay2017_LD/variants.tsv
sed 1d output/LC_McKay2017_LD/variants.tsv > output/LC_McKay2017_LD/variants.bed

##### KNOWN GENES #####

# Dietlein2020 #
get_kgs LC \
<(  cat data/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
    awk -F'\t' -vOFS='\t' '$2 ~ "Lung" && $4 != "level D" {print $1}' )


