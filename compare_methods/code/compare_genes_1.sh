#!/bin/bash
module load bedtools/2.27.1

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
  echo -e "variant\tcs\tsymbol\tscore\tmethod" > $out_dir/predictions.tsv
    
  # tgp #
  cat $base_dir/tgp/out/$trait/$celltypes/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 {print $5,$1,$6,$7,"tgp"}' |
  sort -u >> $out_dir/predictions.tsv
  
  # closest #
  bedtools closest \
    -a $variants_bed  \
    -b <(zcat $base_dir/tgp_paper/wrangle_package_data/reference_panels/data/GENCODE/gencode.v34lift37.basic.tss.bed.gz) |
  cut -f4,5,10 |
  awk '{print $0"\t1\tclosest"}' | 
  sort -u >> $out_dir/predictions.tsv
  
  # ABC_all_celltypes 
  # from cS2G: we used ABC links, kept the maximum ABC score across the 167 cell-types 
  # when an enhancer was interacting with a gene in multiple cell-types, and used this 
  # score as raw linking value.
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz | 
  cut -f1-3,7,21 |
  sort -k5,5nr | 
  sort -u -k1,1 -k2,2n -k4,4 |
  sort -k1,1 -k2,2n |
  bedtools intersect -a $variants_bed -b - -wa -wb |
  cut -f4,5,9,10 |
  awk '{print $0"\tABC_all_celltypes"}' | 
  sort -u >> $out_dir/predictions.tsv
  
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
  (> $out_dir/EpiMAP_corrlinks.tmp.gz ; for file in data/EpiMAP/links/links_corr_only/BSS*_collated_pred.tsv.gz ; do 
    echo $file
    join -1 6 -2 4 -t$'\t' -o 1.1,1.2,1.3,1.4,1.5,2.5,1.7 \
    <(  zcat $file | 
        bedtools intersect -a - -b $variants_bed -wa -wb | 
        awk -F'\t' -vOFS='\t' '{print $1,$2,$3,$10,$11,$4,$5}' | 
        sort -k6,6 ) \
    <(  zcat $base_dir/tgp_paper/wrangle_package_data/reference_panels/data/GENCODE/gencode.v34lift37.basic.tss.bed.gz | 
        sort -k4,4 ) |
    gzip >> $out_dir/EpiMAP_corrlinks_collapsed.bed.gz
  done)
  zcat $out_dir/EpiMAP_corrlinks_collapsed.bed.gz |
  sort -k7,7nr | 
  sort -u -k1,1 -k2,2n -k6,6 | 
  cut -f4- |
  awk '{print $0"\tEpiMAP_all_celltypes"}' | 
  sort -u >> $out_dir/predictions.tsv
  
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
  
) #| cat > $out_dir/predictions.tsv

