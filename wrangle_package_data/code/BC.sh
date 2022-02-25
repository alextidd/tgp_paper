#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/Traits/BC/ ; mkdir -p $dataDir
outDir=output/Traits/BC/ ; mkdir -p $outDir

# fix CS names function
function sed_ichav () {
  cat $1 | sed 's/\tichav/\./g' | sed 's/\tCIMBA/\.CIMBA/g'
}

# variant list
cp \
  /working/lab_georgiat/jonathB/bcac_oncoarray/BCACFM.CCRV.STRONG.variants.bed \
  $dataDir
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  cat ${dataDir}/BCACFM.CCRV.STRONG.variants.bed |
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$7,$8}'  |
  sed_ichav |
  sort -k1,1 -k2,2n ) |
cat > $outDir/BC.VariantList.tsv
sed 1d $outDir/BC.VariantList.tsv > $outDir/BC.VariantList.bed

# expanded variants
expandedOutputDir=output/Traits/BC_expanded/ ; mkdir -p $expandedOutputDir
( cat data/Traits/BC/BCACFM.CCRV.STRONG.variants.bed | 
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$7,$8,"STRONG"}'|
  sed_ichav  ;
  cat /working/lab_georgiat/jonathB/bcac_oncoarray/BCAC.secondarySignalSNPs.9col.bed | 
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$8,$9,"SECONDARY"}'  |
  sed_ichav ; ) | 
cat | sort -k1,1 -k2,2n > $expandedOutputDir/BC_expanded.VariantList.bed
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tSignal" ;
  cat $expandedOutputDir/BC_expanded.VariantList.bed ) |
cat > $expandedOutputDir/BC_expanded.VariantList.tsv
  
# genes list
cp \
  /working/lab_georgiat/alexandT/target_gene_prediction/data/BC/drivers/breast_cancer_drivers_2021.txt \
  $dataDir
cat ${dataDir}/breast_cancer_drivers_2021.txt \
> $outDir/BC.Drivers.txt
cat ${dataDir}/breast_cancer_drivers_2021.txt \
> $expandedOutputDir/BC_expanded.Drivers.txt

# best SNP per CS
bestSNPs=/working/lab_georgiat/jonathB/PROJECTS/trench_lab/BC_riskSNP_targets/output/nearest_gene/BCrisk/BCAC_196signals_bestPval.bed
( echo -e "BestSNPPos\tBestSNP\tCredibleSet" ;
  cat $bestSNPs |
  sed_ichav |
  cut -f3,4,5 ) |
cat > $outDir/BC.best_SNP_per_CS.tsv

# credible set list
module load R/4.0.2
Rscript code/CredibleSetList.R BC
