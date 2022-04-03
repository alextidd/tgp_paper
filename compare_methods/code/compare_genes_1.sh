#!/bin/bash

# Michailidou2017 x BRST.MCF7.CNCR method gene predictions #
trait=BC_Michailidou2017_FM
celltypes=BRST.MCF7.CNCR_celltype # (enriched_tissues BRST.MCF7.CNCR_celltype all_celltypes)

base_dir=/working/lab_jonathb/alexandT/ ; 
WKDIR=$base_dir/tgp_paper/compare_methods/ ; cd $WKDIR
out_dir=output/$trait/$celltypes/ ; mkdir -p $out_dir

variants_bed=$base_dir/tgp_paper/wrangle_package_data/traits/output/$trait/variants.bed
traits_metadata=$base_dir/tgp_paper/wrangle_package_data/traits/output/metadata.tsv
PMID=$(cat $traits_metadata | awk -v trait="$trait" '$2 == trait {print $5}')
EpiMAP_BSS=$(cat data/EpiMAP/metadata/main_metadata_table.tsv | awk '$3=="MCF-7" {print $1}')

(
  echo -e "cs\tsymbol\tscore\tmethod" ;
  
  # ABC #
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz | 
  awk '$24 == "MCF-7-ENCODE"' |
  sort -k1,1 -k2,2n |
  bedtools intersect -a $variants_bed -b - -wa -wb |
  awk -F'\t' -vOFS='\t' '{print $5,$12,$26,"ABC"}' | 
  sort -u ;
  
  # EpiMAP #
  join -1 9 -2 4 -t$'\t' -o 1.5,2.5,1.10 \
  <(  zcat data/EpiMAP/links/links_corr_only/${EpiMAP_BSS}_collated_pred.tsv.gz |
      bedtools intersect -a $variants_bed -b - -wa -wb | sort -k9,9 ) \
  <(  zcat $base_dir/tgp_paper/wrangle_package_data/reference_panels/data/GENCODE/gencode.v34lift37.basic.tss.bed.gz |
      sort -k4,4 ) |
  awk '{print $0"\tEpiMAP"}' | 
  sort -u ;
  
  # tgp #
  cat $base_dir/tgp/out/$trait/$celltypes/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 {print $1,$6,$7,"tgp"}' | 
  sort -u ;
  
  # MAGMA #
  cat data/MAGMA/genes.txt |
  awk -F'\t' -vOFS='\t' '{print "NA",$0,"1","MAGMA"}' | 
  sort -u ;
  
  # cS2G #
  cat data/cS2G/gwas_catalog_cS2G/gwas_catalog_v1.0-associations_e100_r2021-02-25.annot | 
  awk -F'\t' -vOFS='\t' -v PMID="$PMID" '$1 == PMID {print "NA",$6,$7,"cS2G"}' | 
  sort -u ;
  
  # closest #
  bedtools closest \
    -a $variants_bed  \
    -b <(zcat $base_dir/tgp_paper/wrangle_package_data/reference_panels/data/GENCODE/gencode.v34lift37.basic.tss.bed.gz) |
  awk -F'\t' -vOFS='\t' '{print $5,$10,"1","closest"}' | 
  sort -u ;

) | cat > $out_dir/genes.tsv

