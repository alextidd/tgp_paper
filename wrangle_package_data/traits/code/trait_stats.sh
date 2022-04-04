#!/bin/bash
base_dir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/

(while IFS=$'\t' read -r trait variants variants_file known_genes_file PMID ; do
  echo $trait $variants 
  cat $variants_file |
    bedtools slop -i - -g $base_dir/reference_panels/data/hg/hg19.genome -b 2000000 |
    bedtools intersect \
      -a - \
      -b  <( join -o1.1,1.2,1.3,1.4,1.5 -1 5 -2 1 -t$'\t' \
              <(zcat $base_dir/reference_panels/output/GENCODE/proteincoding.gencode.v34lift37.basic.tss.bed.gz | 
                sort -k5,5) \
              <(cat $known_genes_file | sort -k1,1) |
              sort -k1,1 -k2,2n ) \
      -wa -wb |
    cut -f4,5,10 |
    sort -u \
    > pairs.tmp
  
  # variants #
  n_variants=$(cat $variants_file | wc -l)
  n_CSs=$(cat $variants_file | sort -k5,5 -u | wc -l)
  avg_n_variants_per_CS=$(echo "scale=1; $n_variants/$n_CSs" | bc -l)
  variants_hdr="variants\tn_variants\tn_CSs\tavg_n_variants_per_CS"
  variants_stats="$variants\t$n_variants\t$n_CSs\t$avg_n_variants_per_CS"
  
  # known genes #
  n_known_genes=$(cat $known_genes_file | wc -l)
  
  # pairs #
  n_vxg_pairs=$(          cat pairs.tmp | wc -l)
  n_cxg_pairs=$(          sort -k2,2 -k3,3 -u pairs.tmp | wc -l)
  n_paired_variants=$(    sort -k1,1 -u pairs.tmp | wc -l)
  n_paired_CSs=$(         sort -k2,2 -u pairs.tmp | wc -l)
  n_paired_known_genes=$( sort -k3,3 -u pairs.tmp | wc -l)
  pairs_hdr="n_vxg_pairs\tn_cxg_pairs\tn_paired_variants\tn_paired_CSs\tn_paired_known_genes"
  pairs_stats="$n_vxg_pairs\t$n_cxg_pairs\t$n_paired_variants\t$n_paired_CSs\t$n_paired_known_genes"
  
  hdr="$variants_hdr\tn_known_genes\t$pairs_hdr"
  ( 
    echo -e $hdr ; 
    echo -e "$variants_stats\t$n_known_genes\t$pairs_stats" ; 
  ) | cat > output/$variants/trait_stats.tsv
  
  # unite trait info tables
  ( echo -e $hdr ; 
    awk -v hdr="$hdr" '$0!=hdr' output/*/trait_stats.tsv ) |
  cat  > output/trait_stats.tsv
  
  rm -f pairs.tmp
  
done < <(sed 1d output/metadata.tsv) )


