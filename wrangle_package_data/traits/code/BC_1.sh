#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Michailidou2017_FM #
# col 4 = snp_onco = the SNP ID from the oncoarray chip
# col 5 = snp_icogs = the SNP ID from the icogs chip
# col 6 = rsID = the dbSNP rsID.
# The reason for the difference is that there are 6 SNPs without rsIDs, 
# so have unknownrsID in col6. There are also 2 ambiguous rsIDs where the 
# variant is an indel and another where there are G/A/C alternate alleles. 
mkdir output/BC_Michailidou2017_FM/
BCACFM_variants=/working/lab_georgiat/jonathB/bcac_oncoarray/BCACFM.CCRV.STRONG.variants.bed
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  awk -F'\t' -v OFS='\t' 'NR==FNR{a[$0];next} { 
    if ($6 in a) 
      print $1,$2,$3,$4,$7,$8 ; 
    else 
      print $1,$2,$3,$6,$7,$8 ;
    }' \
    <(cat $BCACFM_variants | cut -f6 | uniq -d) \
    <(cat $BCACFM_variants) |
  sed_ichav |
  sort -k1,1 -k2,2n ) |
cat > output/BC_Michailidou2017_FM/variants.tsv
sed 1d output/BC_Michailidou2017_FM/variants.tsv > output/BC_Michailidou2017_FM/variants.bed

# Michailidou2017_LD #
mkdir output/BC_Michailidou2017_LD/
# 1. get rsIDs
(
# 1a. get index rsIDs for existing loci (ST14)
  cat data/Michailidou2017/Supplementary_Table_14.tsv |
  awk '$1 == $8 && $7 != "NULL" {print $7}' ;
# 1b. get index rsIDs for novel loci (ST7)
  cat data/Michailidou2017/Supplementary_Table_7.tsv |
  sed 1d | cut -f4 ;
) | cat | sort -u > output/BC_Michailidou2017_LD/index_SNPs.txt
# 1c. get missing rsIDs
module load R/4.0.2 ; Rscript code/BC_2.R ; 
# 1d. input missing_index_rsIDs into UCSC Table Browser > found_index_rsIDs...
cat output/BC_Michailidou2017_LD/found_index_rsIDs.tsv |
cut -f4 | sort -u >> output/BC_Michailidou2017_LD/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants output/BC_Michailidou2017_LD/

# cS2GxMichailidou2017_assoc # - treat each variant as a credible set
mkdir output/BC_cS2GxMichailidou2017_assoc/
cS2G_to_variants BC Michailidou2017

# L2GxMichailidou2017_FM #
mkdir output/BC_L2GxMichailidou2017_FM/
L2G_to_variants BC Michailidou2017

# # Michailidou2017_LD - # the old index SNPs
# cut -f4 /working/lab_georgiat/jonathB/PROJECTS/trench_lab/BC_riskSNP_targets/output/nearest_gene/BCrisk/BCAC_196signals_bestPval.bed 

# # Michailidou2017_summstats #
# mkdir output/BC_Michailidou2017_summstats/
# ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
#   zcat data/GWASCatalog/Michailidou2017/harmonised/29059683-GCST004988-EFO_0000305-build37.f.tsv.gz |
#   awk -F'\t' -v OFS='\t' 'NR>1 && $9<0.000000005 {print "chr"$2,$3-1,$3,$1,$1}' ) |
# cat > output/BC_Michailidou2017_summstats/variants.tsv
# sed 1d output/BC_Michailidou2017_summstats/variants.tsv > output/BC_Michailidou2017_summstats/variants.bed

##### KNOWN GENES #####

# Fachal2020 #
get_kgs BC \
<(  cat /working/lab_georgiat/alexandT/target_gene_prediction/data/BC/drivers/breast_cancer_drivers_2021.txt )



#########
# # best SNP per CS
# bestSNPs=/working/lab_georgiat/jonathB/PROJECTS/trench_lab/BC_riskSNP_targets/output/nearest_gene/BCrisk/BCAC_196signals_bestPval.bed
# ( echo -e "chrom\tBestSNPPos\tBestSNP\tCredibleSet" ;
#   cat $bestSNPs |
#   sed_ichav |
#   cut -f1,3,4,5 ) |
# cat > $outDir/BC/best_SNP_per_CS.tsv
# 
# # credible set list
# module load R/4.0.2
# Rscript code/CredibleSetList.R BC

# # expanded variants
# ( cat data/BC/BCACFM.CCRV.STRONG.variants.bed | 
#   awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$7,$8,"STRONG"}'|
#   sed_ichav  ;
#   cat /working/lab_georgiat/jonathB/bcac_oncoarray/BCAC.secondarySignalSNPs.9col.bed | 
#   awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$8,$9,"SECONDARY"}'  |
#   sed_ichav ; ) | 
# cat | sort -k1,1 -k2,2n > $expandedOutDir/VariantList.bed
# ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tSignal" ;
#   cat $expandedOutDir/VariantList.bed ) |
# cat > $expandedOutDir/VariantList.tsv
