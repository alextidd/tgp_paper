#!/bin/bash
# generate directories and define outfiles
base_dir=/working/lab_jonathb/alexandT/
cd $base_dir/tgp_paper/wrangle_package_data/traits/

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
  cat $kgs > ${trait}_kgs.tmp
  (for trait_dir in output/${trait}_*/ ; do 
    cat ${trait}_kgs.tmp > $trait_dir/known_genes.txt ; 
  done)
  rm -f ${trait}_kgs.tmp
}

L2G_to_variants () {
  trait=$1 ; study=$2
  hg19_L2G=$(ls $base_dir/tgp_paper/compare_methods/data/L2G/$study/GCST*-independently-associated-loci.hg19.bed)
  out_dir=output/${trait}_L2Gx${study}_FM/
  
  if [[ ! -f $hg19_L2G ]] ; then 
    echo "must run $base_dir/tgp_paper/compare_methods/data/L2G/liftover.sh first" 
  else 
    ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
      cat $hg19_L2G |
      awk -F'\t' -vOFS='\t' '{print $1,$2,$3,$4,$4}' |
      sort -u ; ) |
    cat > $out_dir/variants.tsv
    sed 1d $out_dir/variants.tsv > $out_dir/variants.bed
  fi
}

cS2G_to_variants () {
  trait=$1 ; study=$2
  hg19_cS2G=$base_dir/tgp_paper/compare_methods/data/cS2G/gwas_catalog_cS2G/gwas_catalog_v1.0-associations_e100_r2021-02-25.annot.hg19.bed
  out_dir=output/${trait}_cS2Gx${study}_assoc/
  PMID=$(cat output/metadata.tsv | awk -F'_' -v study=$study '$2==study {print}' | cut -f5)
  
  if [ $PMID == '' ] ; then
    echo "no PMID found for $study in output/metadata.tsv"
  elif [ ! -f $hg19_cS2G ] ; then 
    echo "$hg19_cS2G does not exist ; must run $base_dir/tgp_paper/compare_methods/data/cS2G/liftover.sh first" 
  elif [ $(grep -c $PMID $hg19_cS2G) -eq 0 ] ; then
    echo "PMID $PMID is not found in $hg19_cS2G"
  else
    ( echo -e "chrom\tstart\tend\tvariant\tCredibleSet" ;
      cat $hg19_cS2G |
      awk -F'\t' -vOFS='\t' -v PMID=$PMID '$4 == PMID {print $1,$2,$3,$8,$8}' |
      sort -u ; ) |
    cat > $out_dir/variants.tsv
    sed 1d $out_dir/variants.tsv > $out_dir/variants.bed
  fi
}