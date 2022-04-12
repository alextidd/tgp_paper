#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Jones2020 FM #
mkdir output/EOC_Jones2020_FM/
( echo "chrom\tstart\tend\tvariant\tcs" ;
  cat data/Jones2020/mmc2_ST1.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 {print "chr"$8,$9-1,$9,$7,$2}' |
  sort -u ;
) > output/EOC_Jones2020_FM/variants.tsv
sed 1d output/EOC_Jones2020_FM/variants.tsv > output/EOC_Jones2020_FM/variants.bed

##### KNOWN GENES #####

# Yamulla2020 #
get_kgs EOC \
<(  cat data/Yamulla2020/Table1.tsv |
    sed 1d | cut -f1 | tr a-z A-Z )
    
    

# # Phelan2017 LD #
# mkdir output/EOC_Phelan2017_LD/
# # ( cd data/Phelan2017/ ; unzip Phelan_Archive.zip )
# cat data/Phelan2017/Summary_chr*.txt | 
# awk -F, -vOFS='\t' '$2~"rs" && $12<0.000000005 && $12!="-99" {print $2}' |
# sed 's/:.*//g' \
# > output/EOC_Phelan2017_LD/index_SNPs.txt
# # run SNiPA...
# SNiPA_to_variants output/EOC_Phelan2017_LD/
