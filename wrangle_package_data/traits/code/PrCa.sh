#!/bin/bash
. code/utils.sh 

##### VARIANTS #####

# Benafif2018_LD #
mkdir output/PrCa_Benafif2018_LD/
# 1. get index SNPs 
cat data/Benafif2018/Table1.tsv |
sed 1d | cut -f1 \
> output/PrCa_Benafif2018_LD/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants output/PrCa_Benafif2018_LD/

# L2GxSchumacher2018_FM #
mkdir output/PrCa_L2GxSchumacher2018_FM/
L2G_to_variants PrCa Schumacher2018

# cS2GxSchumacher0271_FM #
mkdir output/PrCa_cS2GxSchumacher2018_assoc/
cS2G_to_variants PrCa Schumacher2018

##### KNOWN GENES #####

# Giambartolomei2021 #
get_kgs PrCa \
<(  cat data/Giambartolomei2021/PrCa_GeneList_Used.csv |
    awk -F, 'NR > 1 {print $1}'  |
    sort -u )


##### GRAVEYARD #####

# # Schumacher2018_summstats #
# mkdir output/PrCa_Schumacher2018_summstats/
# ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
#   zcat data/GWASCatalog/Schumacher2018/harmonised/29892016-GCST006085-EFO_0001663-build37.f.tsv.gz |
#   awk -F'\t' -v OFS='\t' 'NR>1 && $9<0.000000005 {print "chr"$2,$3-1,$3,$1,$1}' ) |
# cat > output/PrCa_Schumacher2018_summstats/variants.tsv
# sed 1d output/PrCa_Schumacher2018_summstats/variants.tsv > output/PrCa_Schumacher2018_summstats/variants.bed

# # Schumacher2018_LD #
# mkdir output/PrCa_Schumacher2018_LD/
# # 1.  get index SNPs from CSs
# cat /working/lab_georgiat/jonathB/PROJECTS/fredwards/multicancer_captureHiC/data/assocsignal_dist/SNPs/prostate.bed |
# cut -f4 | sort -u \
# > output/PrCa_Schumacher2018_LD/index_SNPs.txt
# # 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
# #     https://snipa.helmholtz-muenchen.de/snipa3/
# # 3.  generate variant list
# SNiPA_to_variants output/PrCa_Schumacher2018_LD/

# # GWASCatalog_LD #
# # 1.  get index SNPs (p<5e-8)
# mkdir output/PrCa_GWASCatalog_LD/
# awk -F'\t' -vOFS='\t' \
#   '$8 == "Prostate cancer" && ($28 + 0) <= 5E-8 {print $22}' \
#   data/GWASCatalog/gwas_catalog_v1.0-associations_e105_r2022-03-08.tsv |
# sort -u |
# sed 's/;\s/\n/g' \
# > output/PrCa_GWASCatalog_LD/index_SNPs.txt
# # 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
# #     https://snipa.helmholtz-muenchen.de/snipa3/
# # 3.  generate variant list
# SNiPA_to_variants output/PrCa_GWASCatalog_LD/

# # Giambartolomei2021_FM #
# mkdir output/PrCa_Giambartolomei2021_FM/
# ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
#   zcat data/Giambartolomei2021/credsets95_paintor_1causals_137regions.txt.gz | 
#   # chr23 -> chrX
#   sed 's/chr23/chrX/g' | 
#   # get coords
#   awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$15,$14}' | 
#   sort -k1,1 -k2,2n ) |
# cat > output/PrCa_Giambartolomei2021_FM/variants.tsv
# sed 1d output/PrCa_Giambartolomei2021_FM/variants.tsv > output/PrCa_Giambartolomei2021_FM/variants.bed
# 
# # Giambartolomei2021_LD #
# mkdir output/PrCa_Giambartolomei2021_LD/
# # 1.  get index SNPs from CSs
# zcat data/Giambartolomei2021/credsets95_paintor_1causals_137regions.txt.gz |
# sed 1d | cut -f17 | sort -u \
# > output/PrCa_Giambartolomei2021_LD/index_SNPs.txt
# # 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
# #     https://snipa.helmholtz-muenchen.de/snipa3/
# # 3.  generate variant list
# SNiPA_to_variants output/PrCa_Giambartolomei2021_LD/
# 
# # Dadaev2018_FM #
# mkdir output/PrCa_Dadaev2018_FM/
# ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
#   # FM variants
#   cat data/Dadaev2018/Supplementary_Data_1a.tsv | 
#   awk -F'\t' -v OFS='\t' 'NR > 1 {print "chr"$5,$6-1,$6,$4,$1,$17}' |
#   # chr23 -> chrX
#   sed 's/chr23/chrX/g' | 
#   # sort
#   sort -k1,1 -k2,2n ) |
# cat > output/PrCa_Dadaev2018_FM/variants.tsv
# sed 1d output/PrCa_Dadaev2018_FM/variants.tsv > output/PrCa_Dadaev2018_FM/variants.bed
