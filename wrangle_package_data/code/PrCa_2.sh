#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/Traits/PrCa/ ; mkdir -p $dataDir/
outDir=output/Traits/ ; mkdir -p $outDir/
. code/SNiPA_to_VarList.sh

# Giambartolomei variant list
# from Giambartolomei et al., 2021
# https://github.com/bogdanlab/hichip/blob/main/additional/credsets95_paintor_1causals_137regions.txt.gz
mkdir -p $outDir/PrCa_Giambartolomei2021/
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  zcat $dataDir/Giambartolomei2021/credsets95_paintor_1causals_137regions.txt.gz | 
  # chr23 -> chrX
  sed 's/chr23/chrX/g' | 
  # get coords
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$15,$14}' | 
  sort -k1,1 -k2,2n ) |
cat > $outDir/PrCa_Giambartolomei2021/VariantList.tsv
sed 1d $outDir/PrCa_Giambartolomei2021/VariantList.tsv > $outDir/PrCa_Giambartolomei2021/VariantList.bed

# Giambartolomei2021_expanded variant list (proxy SNPs of the index SNP)
mkdir -p $outDir/PrCa_Giambartolomei2021_expanded/SNiPA/
# 1.  get index SNPs from CSs
zcat $dataDir/Giambartolomei2021/credsets95_paintor_1causals_137regions.txt.gz |
sed 1d | cut -f17 | sort -u \
> $outDir/PrCa_Giambartolomei2021_expanded/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_VarList PrCa_Giambartolomei2021_expanded

# Beesley LD list
mkdir -p $outDir/PrCa_Beesley/SNiPA/
cat /working/lab_georgiat/jonathB/PROJECTS/fredwards/multicancer_captureHiC/data/assocsignal_dist/SNPs/prostate.bed |
cut -f4 | sort -u \
> $outDir/PrCa_Beesley/index_SNPs.txt
# SNiPA!
SNiPA_to_VarList PrCa_Beesley

# Dadaev variant list
# from Dadaev et al., 2018
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-04109-8/MediaObjects/41467_2018_4109_MOESM4_ESM.xlsx
# Supplementary Data 1
# (a) Full list of variants within the 95% credible set selected by JAM
# (b) List of all variants that had a marginal P-value beyond the threshold
#     specified for genome-wide significance of association with PrCa (P<5x10-8) in the 5 regions for which
#     JAM was unable to resolve candidate variants, alongside detailed variant annotations. 
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  # FM variants
  cat $dataDir/Dadaev2018/Supplementary_Data_1a.tsv | 
  awk -F'\t' -v OFS='\t' 'NR > 1 {print "chr"$5,$6-1,$6,$4,$1,$17}' |
  # chr23 -> chrX
  sed 's/chr23/chrX/g' | 
  # sort
  sort -k1,1 -k2,2n
  ) |
cat > $outDir/PrCa_Dadaev2018/VariantList.tsv
sed 1d $outDir/PrCa_Dadaev2018/VariantList.tsv > $outDir/PrCa_Dadaev2018/VariantList.bed

# GWAS catalog variants
# 1.  get GWAS catalog associated variants list
#     from NHGRI-EBI: https://www.ebi.ac.uk/gwas/efotraits/EFO_0001663 
ls $dataDir/GWAS_catalog/efotraits_EFO_0001663-associations-2022-03-4.csv
# 2.  get index SNPs (p<5e-8)
module load R/4.0.2 ; Rscript code/PrCa_1.R
# 3.  get proxy SNPs (r2>0.8 w/ an index SNP) 
#     from https://snipa.helmholtz-muenchen.de/snipa3/
SNiPA_to_VarList PrCa_GWAS_catalog

# gene list
# https://github.com/bogdanlab/hichip/blob/main/additional/PrCa_GeneList_Used.csv
cat $dataDir/Giambartolomei2021/PrCa_GeneList_Used.csv |
awk -F, 'NR > 1 {print $1}'  |
sort -u |
tee $outDir/PrCa_Giambartolomei2021/Drivers.txt |
tee $outDir/PrCa_Dadaev2018/Drivers.txt |
tee $outDir/PrCa_Giambartolomei2021_expanded/Drivers.txt |
tee $outDir/PrCa_Beesley/Drivers.txt \
> $outDir/PrCa_GWAS_catalog/Drivers.txt


# credible set list
module load R/4.0.2
Rscript code/CredibleSetList.R PrCa_Giambartolomei2021
Rscript code/CredibleSetList.R PrCa_Dadaev2018
# Rscript code/CredibleSetList.R PrCa_GWAS_catalog
