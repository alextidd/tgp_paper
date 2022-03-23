#!/bin/bash
# # hpcapp01 head node: #
# (
# IBDDir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/traits/IBD/
# mkdir -p $IBDDir/{variants,known_genes}/Huang2017/ ; cd $IBDDir
# FTPPath=ftp://ftp.broadinstitute.org/outgoing/lincRNA/Nasser2020/data/
# wget -P known_genes/Huang2017/ $FTPPath/GeneLists.IBD.txt
# wget -P variants/Huang2017/ $FTPPath/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt
# wget -P variants/Huang2017/ $FTPPath/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt
# echo $FTPPath | tee variants/Huang2017/README.txt > known_genes/Huang2017/README.txt
# )

WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/traits/IBD/
outDir=output/traits/IBD/ ; mkdir -p $outDir/{variants,known_genes}/
. code/example_data/utils.sh

##### KNOWN GENES #####

# Huang2017 #
mkdir $outDir/known_genes/Huang2017/
cat $dataDir/known_genes/Huang2017/GeneLists.IBD.txt |
awk '$0!="" && NR > 1' |
# fix gencode incompatible C1orf106 -> INAVA
sed 's/C1orf106/INAVA/g' \
> $outDir/known_genes/Huang2017/known_genes.txt

##### VARIANTS #####

# Huang2017 #
mkdir $outDir/variants/Huang2017/
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  cat $dataDir/variants/Huang2017/IBDCombined.set1-2.variant.list.txt |
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$4,$7}' |
  sort -k1,1 -k2,2n ) | 
cat > $outDir/variants/Huang2017/variants.tsv
sed 1d $outDir/variants/Huang2017/variants.tsv > $outDir/variants/Huang2017/variants.bed

# Huang2017_proxies #
mkdir -p $outDir/variants/Huang2017_proxies/SNiPA/
# 1.  get index (best) SNPs from CSs
cat $dataDir/variants/Huang2017/IBDCombined.set1-2.variant.list.txt |
sed 1d | cut -f8 | sort -u \
> $outDir/variants/Huang2017_proxies/SNiPA/best_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants IBD Huang2017_proxies
