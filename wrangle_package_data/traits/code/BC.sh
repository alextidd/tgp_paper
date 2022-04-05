#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Michailidou2017 FM #
mkdir output/BC_Michailidou2017_FM/
( echo -e "chrom\tstart\tend\tvariant\tcs" ;
  cat /working/lab_georgiat/jonathB/bcac_oncoarray/BCACFM.CCRV.STRONG.variants.bed |
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$7,$8}'  |
  sed_ichav |
  sort -k1,1 -k2,2n ) |
cat > output/BC_Michailidou2017_FM/variants.tsv
sed 1d output/BC_Michailidou2017_FM/variants.tsv > output/BC_Michailidou2017_FM/variants.bed

# Michailidou2017 LD #
mkdir output/BC_Michailidou2017_LD/
# 1.  get index/best SNPs from CSs
cat /working/lab_georgiat/jonathB/PROJECTS/trench_lab/BC_riskSNP_targets/output/nearest_gene/BCrisk/BCAC_196signals_bestPval.bed | 
cut -f4 |
sed 's/\:.*//g' \
> output/BC_Michailidou2017_LD/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants output/BC_Michailidou2017_LD/

# # Michailidou2017 summstats #
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
