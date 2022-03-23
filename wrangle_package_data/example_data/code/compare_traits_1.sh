#!/bin/bash

# get trait info
VarList=output/Traits/$trait/VariantList.bed
n=$(cat $VarList |
  bedtools slop -i - -g data/hg/hg19.genome -b 2000000 |
  bedtools intersect \
    -a - \
    -b  <( join -o1.1,1.2,1.3,1.4,1.5 -1 5 -2 1 -t$'\t' \
            <(zcat output/GENCODE/proteincoding.gencode.v34lift37.basic.tss.bed.gz | sort -k5,5) \
            <(cat output/Traits/$trait/Drivers.txt | sort -k1,1) ) \
    -wa -wb |
  sort -k5,5 -k11,11 -u |
  wc -l)
nVars=$(cat $VarList | wc -l)
nCSs=$(cat $VarList | sort -k5,5 -u | wc -l)
avgVarsPerCS=$(echo "scale=1; $nVars/$nCSs" | bc -l)
nGenes=$(cat output/Traits/$trait/Drivers.txt | wc -l)

( echo -e "trait\tnPairs\tnVars\tnCSs\tavgVarsPerCS\tnGenes" ;
  echo -e "$trait\t$n\t$nVars\t$nCSs\t$avgVarsPerCS\t$nGenes"
) | cat > output/Traits/$trait/trait_info.tsv

# unite trait info tables
awk 'FNR==1 && NR!=1 { while (/^trait/) getline ; } ; 1 {print}' \
output/Traits/*/trait_info.tsv \
> output/Traits/trait_info.tsv
