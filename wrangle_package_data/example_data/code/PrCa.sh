#!/bin/bash
module load R/4.0.2
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/traits/PrCa/ ; mkdir -p $dataDir/
outDir=output/traits/PrCa/ ; mkdir -p $outDir/{variants,known_genes}/
. code/example_data/utils.sh

##### KNOWN GENES #####

# Giambartolomei2021 #
mkdir $outDir/known_genes/Giambartolomei2021/
cat $dataDir/known_genes/Giambartolomei2021/PrCa_GeneList_Used.csv |
awk -F, 'NR > 1 {print $1}'  |
sort -u \
> $outDir/known_genes/Giambartolomei2021/known_genes.txt

##### VARIANTS #####

# Giambartolomei2021 #
mkdir $outDir/variants/Giambartolomei2021/
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  zcat $dataDir/variants/Giambartolomei2021/credsets95_paintor_1causals_137regions.txt.gz | 
  # chr23 -> chrX
  sed 's/chr23/chrX/g' | 
  # get coords
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$15,$14}' | 
  sort -k1,1 -k2,2n ) |
cat > $outDir/variants/Giambartolomei2021/variants.tsv
sed 1d $outDir/variants/Giambartolomei2021/variants.tsv > $outDir/variants/Giambartolomei2021/variants.bed

# Giambartolomei2021_proxies #
mkdir -p $outDir/variants/Giambartolomei2021_proxies/SNiPA/
# 1.  get index SNPs from CSs
zcat $dataDir/variants/Giambartolomei2021/credsets95_paintor_1causals_137regions.txt.gz |
sed 1d | cut -f17 | sort -u \
> $outDir/variants/Giambartolomei2021_proxies/SNiPA/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants PrCa Giambartolomei2021_proxies

# Dadaev2018 #
mkdir $outDir/variants/Dadaev2018/
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  # FM variants
  cat $dataDir/variants/Dadaev2018/Supplementary_Data_1a.tsv | 
  awk -F'\t' -v OFS='\t' 'NR > 1 {print "chr"$5,$6-1,$6,$4,$1,$17}' |
  # chr23 -> chrX
  sed 's/chr23/chrX/g' | 
  # sort
  sort -k1,1 -k2,2n ) |
cat > $outDir/variants/Dadaev2018/variants.tsv
sed 1d $outDir/variants/Dadaev2018/variants.tsv > $outDir/variants/Dadaev2018/variants.bed

# Schumacher2018_summstats #

mkdir $outDir/variants/Schumacher2018_summstats.proxies/SNiPA/
summstats_to_VarList PrCa Schumacher2018_summstats
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  zcat $dataDir/variants/Schumacher2018_summstats/29892016-GCST006085-EFO_0001663-build37.f.tsv.gz |
  awk -F'\t' -v OFS='\t' 'NR>1 {print "chr"$2,$3-1,$3,$1,$1}' ) |
cat > $outDir/variants/Schumacher2018_summstats/variants.tsv
sed 1d $outDir/variants/Schumacher2018_summstats/variants.tsv > $outDir/variants/Schumacher2018_summstats/variants.bed

# # Schumacher2018_summstats.proxies #
# mkdir $outDir/variants/Schumacher2018_summstats.proxies/SNiPA/
# cat $outDir/variants/Schumacher2018_summstats/variants.bed | 
# cut -f4

# Schumacher2018_summstats.GCTA-COJO.proxies #
mkdir -p $outDir/variants/Schumacher2018_summstats.GCTA-COJO.proxies/SNiPA/
# 1.  get index SNPs from CSs
cat /working/lab_georgiat/jonathB/PROJECTS/fredwards/multicancer_captureHiC/data/assocsignal_dist/SNPs/prostate.bed |
cut -f4 | sort -u \
> $outDir/variants/Schumacher2018_summstats.GCTA-COJO.proxies/SNiPA/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants PrCa Schumacher2018_summstats.GCTA-COJO.proxies

# # GWASCatalog_associations.proxies #
# mkdir -p $outDir/variants/GWASCatalog_associations.proxies/SNiPA/
# Rscript code/example_data/associations_to_index_SNPs.R PrCa

# GWASCatalog_associations.proxies #
# 1.  get index SNPs (p<5e-8)
mkdir -p $outDir/variants/GWASCatalog_associations.proxies/SNiPA/
awk -F'\t' -vOFS='\t' \
  '$8 == "Prostate cancer" && ($28 + 0) <= 5E-8 {print $22}' \
  data/GWASCatalog/gwas_catalog_v1.0-associations_e105_r2022-03-08.tsv |
sort -u |
sed 's/;\s/\n/g' \
> $outDir/variants/GWASCatalog_associations.proxies/SNiPA/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants PrCa GWASCatalog_associations.proxies
