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

# wrangle
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
dataDir=data/Traits/IBD/
outDir=output/Traits/IBD/ ; mkdir -p $outDir

# variant list
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
  cat $dataDir/IBDCombined.set1-2.variant.list.txt |
  awk -F'\t' -v OFS='\t' '{print $1,$2-1,$2,$3,$4}' |
  sort -k1,1 -k2,2n ) | 
cat > $outDir/IBD.VariantList.tsv
sed 1d $outDir/IBD.VariantList.tsv > $outDir/IBD.VariantList.bed

# gene list
cat $dataDir/GeneLists.IBD.txt |
awk '$0!="" && NR > 1' |
# fix gencode incompatible C1orf106 -> INAVA
sed 's/C1orf106/INAVA/g' \
> $outDir/IBD.Drivers.txt
