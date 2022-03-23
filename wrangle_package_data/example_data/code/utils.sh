#!/bin/bash
summstats_to_index_SNPs () {
  trait=$1 ; vars=$2
  mkdir -p output/traits/$trait/variants/$vars/SNiPA/
  ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
    zcat data/traits/$trait/variants/$vars/*build37.f.tsv.gz |
    awk -F'\t' -v OFS='\t' 'NR>1 && $9<0.000000005 {print "chr"$2,$3-1,$3,$1,$1}' ) |
  cat > output/traits/$trait/variants/$vars/SNiPa/index_SNPs.txt
}

summstats_to_variants () {
  trait=$1 ; vars=$2
  ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
    zcat data/traits/$trait/variants/$vars/*build37.f.tsv.gz |
    awk -F'\t' -v OFS='\t' 'NR>1 {print "chr"$2,$3-1,$3,$1,$1}' ) |
  cat > output/traits/$trait/variants/$vars/variants.tsv
  sed 1d output/traits/$trait/variants/$vars/variants.tsv \
  > output/traits/$trait/variants/$vars/variants.bed
}

SNiPA_to_variants () {
  trait=$1 ; vars=$2
  ( cd output/traits/$trait/variants/$vars/
    ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
      cat SNiPA/proxySearch.results.csv |
      awk -F'\t' -vOFS='\t' 'NR > 1 {print "chr"$4,$6-1,$6,$2,$1}' ) |
    cat > variants.tsv
    sed 1d variants.tsv > variants.bed
  )
}