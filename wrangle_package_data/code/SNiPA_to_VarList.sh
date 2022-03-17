#!/bin/bash

SNiPA_to_VarList () {
  ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet\tPosteriorProb" ;
    cat $outDir/$1/SNiPA/proxySearch.results.csv |
    awk -F'\t' -vOFS='\t' 'NR > 1 {print "chr"$4,$6-1,$6,$2,$1,"1000"}' ) |
  cat > $outDir/$1/VariantList.tsv
  sed 1d $outDir/$1/VariantList.tsv > $outDir/$1/VariantList.bed
}