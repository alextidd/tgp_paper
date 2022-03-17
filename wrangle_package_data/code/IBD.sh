#!/bin/bash

# # wget (hpcapp01 head node)
# (
# IBDDir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/Traits/IBD/ ; mkdir -p $IBDDir ; cd $IBDDir
# FTPPath=ftp://ftp.broadinstitute.org/outgoing/lincRNA/Nasser2020/data/
# wget ${FTPPath}/GeneLists.IBD.txt
# wget ${FTPPath}/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt
# wget ${FTPPath}/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt
# echo $FTPPath > README.txt
# )

WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/Traits/IBD/
outDir=output/Traits/ ; mkdir -p $outDir
. code/SNiPA_to_VarList.sh

# variant list
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  cat $dataDir/IBDCombined.set1-2.variant.list.txt |
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$4}' |
  sort -k1,1 -k2,2n ) | 
cat > $outDir/IBD/VariantList.tsv
sed 1d $outDir/IBD/VariantList.tsv > $outDir/IBD/VariantList.bed

# expanded variant list (proxy SNPs of the index SNP)
mkdir -p $outDir/IBD_index_proxies/SNiPA/
# 1.  get index (best) SNPs from CSs
cat $dataDir/IBDCombined.set1-2.variant.list.txt |
sed 1d | cut -f8 | sort -u \
> $outDir/IBD_index_proxies/best_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_VarList IBD_index_proxies

# gene list
cat $dataDir/GeneLists.IBD.txt |
awk '$0!="" && NR > 1' |
# fix gencode incompatible C1orf106 -> INAVA
sed 's/C1orf106/INAVA/g' |
tee $outDir/IBD_index_proxies/Drivers.txt \
> $outDir/IBD/Drivers.txt


