#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
outDir=output/traits/BC/ ; mkdir -p $outDir/{variants,known_genes}/
. code/example_data/utils.sh
function sed_ichav () { cat $1 | sed 's/\tichav/\./g' | sed 's/\tCIMBA/\.CIMBA/g' ; }

##### KNOWN GENES #####

# Fachal2020 #
mkdir $outDir/known_genes/Fachal2020/
cp /working/lab_georgiat/alexandT/target_gene_prediction/data/BC/drivers/breast_cancer_drivers_2021.txt \
  $outDir/known_genes/Fachal2020/known_genes.txt

##### VARIANTS #####

# Michailidou2017 #
mkdir $outDir/variants/Michailidou2017
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  cat /working/lab_georgiat/jonathB/bcac_oncoarray/BCACFM.CCRV.STRONG.variants.bed |
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$7,$8}'  |
  sed_ichav |
  sort -k1,1 -k2,2n ) |
cat > $outDir/variants/Michailidou2017/variants.tsv
sed 1d $outDir/variants/Michailidou2017/variants.tsv > $outDir/variants/Michailidou2017/variants.bed

# Michailidou2017_proxies #
mkdir -p $outDir/variants/Michailidou2017_proxies/SNiPA/
# 1.  get index/best SNPs from CSs
cat /working/lab_georgiat/jonathB/PROJECTS/trench_lab/BC_riskSNP_targets/output/nearest_gene/BCrisk/BCAC_196signals_bestPval.bed | 
cut -f4 |
sed 's/\:.*//g' \
> $outDir/variants/Michailidou2017_proxies/SNiPA/best_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants BC Michailidou2017_proxies

# Michailidou2017_summstats #
mkdir $outDir/variants/Michailidou2017_summstats/
summstats_to_variants BC Michailidou2017_summstats






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
# ( cat data/Traits/BC/BCACFM.CCRV.STRONG.variants.bed | 
#   awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$7,$8,"STRONG"}'|
#   sed_ichav  ;
#   cat /working/lab_georgiat/jonathB/bcac_oncoarray/BCAC.secondarySignalSNPs.9col.bed | 
#   awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$8,$9,"SECONDARY"}'  |
#   sed_ichav ; ) | 
# cat | sort -k1,1 -k2,2n > $expandedOutDir/VariantList.bed
# ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tSignal" ;
#   cat $expandedOutDir/VariantList.bed ) |
# cat > $expandedOutDir/VariantList.tsv
