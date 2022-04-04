#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Law2019 LD #
mkdir output/CRC_Law2019_LD/
# 1.  get index SNPs from CSs
cat data/Law2019/41467_2019_9775_MOESM4_ESM.tsv |
sed -n '/Principally Asian GWAS/q;p' |
awk -F'\t' -vOFS='\t' 'NR > 6 && $2 != "" {print $2}' |
sed 's/\s\*\*//g' \
> output/CRC_Law2019_LD/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
#     and save in $outDir/CRC_index_proxies/SNiPA/
# 3.  generate variant list
SNiPA_to_variants output/CRC_Law2019_LD/

##### KNOWN GENES #####

# Bailey2018_and_Dietlein2020 #
get_kgs CRC \
<(  (cat data/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
    awk -F'\t' -vOFS='\t' '$2 == "Colorectal" && $4 != "level D" {print $1}' ; 
    cat data/Bailey2018/NIHMS948705-supplement-8_Table_S1.tsv |
    awk -F'\t' -vOFS='\t' '$2 == "COADREAD" {print $1}' ) |
    cat | sort -u )