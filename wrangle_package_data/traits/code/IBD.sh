#!/bin/bash
# # hpcapp01 head node: #
# (
# IBDDir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/data/Huang2017
# mkdir -p $IBDDir/ ; cd $IBDDir
# FTPPath=ftp://ftp.broadinstitute.org/outgoing/lincRNA/Nasser2020/data/
# wget $FTPPath/GeneLists.IBD.txt
# wget $FTPPath/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt
# wget $FTPPath/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt
# echo $FTPPath > $IBDDir/README.txt
# )

. code/utils.sh

##### VARIANTS #####

# Huang2017 FM #
mkdir output/IBD_Huang2017_FM/
( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
  cat data/Huang2017/IBDCombined.set1-2.variant.list.txt |
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $1,$2-1,$2,$3,$4,$7}' |
  sort -k1,1 -k2,2n ) | 
cat > output/IBD_Huang2017_FM/variants.tsv
sed 1d output/IBD_Huang2017_FM/variants.tsv > output/IBD_Huang2017_FM/variants.bed

# Huang2017 LD #
mkdir output/IBD_Huang2017_LD/
# 1.  get index (best) SNPs from CSs
cat data/Huang2017/IBDCombined.set1-2.variant.list.txt |
sed 1d | cut -f8 | sort -u \
> output/IBD_Huang2017_LD/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants output/IBD_Huang2017_LD/

##### KNOWN GENES #####

# Huang2017 #
get_kgs IBD \
<(  cat data/Huang2017/GeneLists.IBD.txt |
    awk '$0!="" && NR > 1' |
    # fix gencode incompatible C1orf106 -> INAVA
    sed 's/C1orf106/INAVA/g' )