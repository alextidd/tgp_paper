#!/bin/bash
module load bedtools/2.27.1
module load R/4.0.2

#trait=BC_Michailidou2017_FM
base_dir=/working/lab_jonathb/alexandT/ ; 
WKDIR=$base_dir/tgp_paper/compare_methods/ ; cd $WKDIR
out_dir=output/$variants/ ; mkdir -p $out_dir

variants_bed=$base_dir/tgp_paper/wrangle_package_data/traits/output/$variants/variants.bed
traits_metadata=$base_dir/tgp_paper/wrangle_package_data/traits/output/metadata.tsv
#PMID=$(cat $traits_metadata | awk -v trait="$trait" '$2 == trait {print $5}')
#EpiMAP_BSS=$(cat data/EpiMAP/metadata/main_metadata_table.tsv | awk '$3=="MCF-7" {print $1}')

(
  echo -e "variant\tcs\tsymbol\tscore\tmethod" ;
    
  # tgp enriched_celltypes #
  cat $base_dir/tgp/out/$variants/enriched_celltypes/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_enriched_celltypes"}' |
  sort -u ;
  
  # tgp enriched tissues #
  cat $base_dir/tgp/out/$variants/enriched_tissues/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_enriched_tissues"}' |
  sort -u ;
  
  # tgp all celltypes #
  cat $base_dir/tgp/out/$variants/all_celltypes/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_all_celltypes"}' |
  sort -u ;
  
  # closest #
  bedtools closest \
    -a $variants_bed  \
    -b <(zcat $base_dir/tgp_paper/wrangle_package_data/reference_panels/data/GENCODE/gencode.v34lift37.basic.tss.bed.gz) |
  cut -f4,5,10 |
  awk '{print $0"\t1\tclosest"}' | 
  sort -u ;
  
  # ABC_all_celltypes 
  # from cS2G: we used ABC links, kept the maximum ABC score across the 167 cell-types 
  # when an enhancer was interacting with a gene in multiple cell-types, and used this 
  # score as raw linking value.
  zcat data/ABC/ABC_AllPredictions_collapsed.bed.gz |
  bedtools intersect -a $variants_bed -b - -wa -wb |
  cut -f4,5,9,10 |
  awk '{print $0"\tABC_all_celltypes"}' | 
  sort -u ;
  
  # # ABC_MCF7 #
  # zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz | 
  # awk '$24 == "MCF-7-ENCODE"' |
  # sort -k1,1 -k2,2n |
  # bedtools intersect -a $variants_bed -b - -wa -wb |
  # awk -F'\t' -vOFS='\t' '{print $1,$2,$3,$5,$12,$26,"ABC"}' | 
  # sort -u ;
  # 
  # # ABC enrichments # (makes predictions for the whole CS, get coords of all input variants)
  # join -1 1 -2 4 -t$'\t' -o 2.1,2.2,2.3,1.1,1.2,1.3,1.4 \
  # <(  cat $base_dir/ABC-GWAS-Paper/ABC-Max/out/ABC/$trait/tgp_settings/GenePredictions.allCredibleSets.tsv |
  #     grep -v '#' | 
  #     awk -F'\t' -vOFS='\t' 'NR > 1 && $6 == "TRUE" && $16 != "NA" {print $2,$9,$16,"ABC_enriched"}' | 
  #     sort -k1,1 ) \
  # <(  cat $variants_bed | cut -f1-3,5 | sort -k4,4 ) |
  # # fix gencode incompatible C1orf106 -> INAVA
  # sed 's/C1orf106/INAVA/g' |
  # sort -u ;
  
  # EpiMAP_all_celltypes #
  # from cS2G: we used linked EpiMap enhancers based on expression-enhancer activity
  # correlation across 833 cell-types, kept the maximum correlation when an enhancer
  # was linked to a gene in multiple tissues, and used this correlation as raw linking
  # value.
  # must run data/EpiMAP/links/EpiMAP_corrlinks_collapsed.sh to generate collapsed links for intersection
  zcat data/EpiMAP/links/EpiMAP_corrlinks_collapsed.bed.gz |
  bedtools intersect -a - -b $variants_bed -wa -wb | 
  awk -F'\t' -vOFS='\t' '{print $9,$10,$4,$5,"EpiMAP_all_celltypes"}' |
  sort -u ;
  
  # # EpiMAP #
  # join -1 9 -2 4 -t$'\t' -o 1.1,1.2,1.3,1.5,2.5,1.10 \
  # <(  zcat data/EpiMAP/links/links_corr_only/${EpiMAP_BSS}_collated_pred.tsv.gz |
  #     bedtools intersect -a $variants_bed -b - -wa -wb | sort -k9,9 ) \
  # <(  zcat $base_dir/tgp_paper/wrangle_package_data/reference_panels/data/GENCODE/gencode.v34lift37.basic.tss.bed.gz |
  #     sort -k4,4 ) |
  # awk '{print $0"\tEpiMAP"}' | 
  # sort -u ;
  # 
  # # EpiMAP GWAS enrichments #
  # zcat data/EpiMAP/gwas_resources/all_gwas_SNP_links.tsv.gz |
  # awk -F'\t' -vOFS='\t' -v uid_PMID="$PMID - " \
  # '$13~uid_PMID && $7 != "NA" && $4 < 2500 {print $1,$2-1,$2,"NA",$6,$7,"EpiMAP_enriched"}' |
  # sort -u ;
  
  # # cS2G #
  # cat data/cS2G/gwas_catalog_cS2G/gwas_catalog_v1.0-associations_e100_r2021-02-25.annot | 
  # awk -F'\t' -vOFS='\t' -v PMID="$PMID" '$1 == PMID {print "chr"$3,$4-1,$4,"NA",$6,$7,"cS2G"}' | 
  # sort -u >> $out_dir/predictions.tsv

  # # MAGMA # - no known genes
  # cat data/MAGMA/genes.txt |
  # awk -F'\t' -vOFS='\t' '{print "NA",$0,"1","MAGMA"}' | 
  # sort -u ;
  
) | cat > $out_dir/predictions.tsv

Rscript code/compare_genes_2.R