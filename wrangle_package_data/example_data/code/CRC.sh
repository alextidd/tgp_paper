#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/traits/CRC/ ; mkdir $dataDir
outDir=output/traits/CRC/ ; mkdir -p $outDir/{variants,known_genes}
. code/example_data/utils.sh

##### KNOWN GENES #####

# Dietlein2020 #
mkdir $outDir/known_genes/Dietlein2020/
cat $dataDir/known_genes/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
awk -F'\t' -vOFS='\t' '$2 == "Colorectal" && $4 != "level D" {print $1}' \
> $outDir/known_genes/Dietlein2020/known_genes.txt

# Bailey2018 #
mkdir $outDir/known_genes/Bailey2018/
cat $dataDir/known_genes/Bailey2018/NIHMS948705-supplement-8_Table_S1.tsv |
awk -F'\t' -vOFS='\t' '$2 == "COADREAD" {print $1}' \
> $outDir/known_genes/Bailey2018/known_genes.txt

# Bailey2018_and_Dietlein2020 #
mkdir $outDir/known_genes/Bailey2018_and_Dietlein2020/
cat $outDir/known_genes/{Bailey2018,Dietlein2020}/known_genes.txt |
sort -u \
> $outDir/known_genes/Bailey2018_and_Dietlein2020/known_genes.txt

##### VARIANTS #####

# Law2019_proxies #
mkdir -p $outDir/variants/Law2019_proxies/SNiPA/
# 1.  get index SNPs from CSs
cat $dataDir/variants/Law2019/41467_2019_9775_MOESM4_ESM.tsv |
sed -n '/Principally Asian GWAS/q;p' |
awk -F'\t' -vOFS='\t' 'NR > 6 && $2 != "" {print $2}' |
sed 's/\s\*\*//g' \
> $outDir/variants/Law2019_proxies/SNiPA/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
#     and save in $outDir/CRC_index_proxies/SNiPA/
# 3.  generate variant list
SNiPA_to_variants CRC Law2019_proxies
