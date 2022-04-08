#!/bin/bash
# generate directories and define outfiles
cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/

SNiPA_to_variants () {
  variants_dir=$1 
  ( cd $variants_dir
    ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
      cat proxySearch.results.csv |
      awk -F'\t' -vOFS='\t' 'NR > 1 {print "chr"$4,$6-1,$6,$2,$1}' |
    sort -k1,1 -k2,2n ) |
    cat  > variants.tsv
    sed 1d variants.tsv > variants.bed
  )
}

function sed_ichav () { cat $1 | sed 's/\tichav/\./g' | sed 's/\tCIMBA/\.CIMBA/g' ; }

get_kgs () {
  trait=$1 ; kgs=$2
  cat $kgs > $trait_kgs.tmp
  (for trait_dir in output/${trait}_*/ ; do 
    cat $trait_kgs.tmp > $trait_dir/known_genes.txt ; 
  done)
  rm -f $trait_kgs.tmp
}

