#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/Traits/PrCa/ ; mkdir -p $dataDir
outDir=output/Traits/PrCa/ ; mkdir -p $outDir

# from Giambortolomei et al., 2021
# https://github.com/bogdanlab/hichip/blob/main/additional/credsets95_paintor_1causals_137regions.txt.gz
# https://github.com/bogdanlab/hichip/blob/main/additional/PrCa_GeneList_Used.csv

# variant list
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  zcat $dataDir/credsets95_paintor_1causals_137regions.txt.gz | 
  # chr23 -> chrX
  sed 's/chr23/chrX/g' | 
  # get coords
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$15,$14}' | 
  sort -k1,1 -k2,2n ) |
cat > $outDir/PrCa.VariantList.tsv
sed 1d $outDir/PrCa.VariantList.tsv > $outDir/PrCa.VariantList.bed

# gene list
cat $dataDir/PrCa_GeneList_Used.csv |
awk -F, 'NR > 1 {print $1}'  \
> $outDir/PrCa.Drivers.txt

# credible set list
module load R/4.0.2
Rscript code/CredibleSetList.R PrCa
