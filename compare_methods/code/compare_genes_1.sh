#!/bin/bash

module load bedtools/2.27.1
module load R/4.0.2

base_dir=/working/lab_jonathb/alexandT/ ; 
WKDIR=$base_dir/tgp_paper/compare_methods/ ; cd $WKDIR
traits_metadata=$base_dir/tgp_paper/wrangle_package_data/traits/output/metadata.tsv

(while IFS=$'\t' read -r trait variants variants_file known_genes_file ; do 
  echo $variants
  out_dir=output/$variants/ ; mkdir -p $out_dir
  variants_bed=$base_dir/tgp_paper/wrangle_package_data/traits/output/$variants/variants.bed
  #PMID=$(cat $traits_metadata | awk -v trait="$trait" '$2 == trait {print $5}')
  #EpiMAP_BSS=$(cat data/EpiMAP/metadata/main_metadata_table.tsv | awk '$3=="MCF-7" {print $1}')
  
  (
    echo -e "variant\tcs\tsymbol\tscore\tmethod" ;
      
    # tgp_enriched_celltypes #
    cat $base_dir/tgp/out/$variants/enriched_celltypes/target_gene_predictions_full.tsv |
    awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_enriched_celltypes"}' |
    sort -u ;
    
    # tgp_enriched_tissues #
    cat $base_dir/tgp/out/$variants/enriched_tissues/target_gene_predictions_full.tsv |
    awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_enriched_tissues"}' |
    sort -u ;
    
    # tgp_all_celltypes #
    cat $base_dir/tgp/out/$variants/all_celltypes/target_gene_predictions_full.tsv |
    awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_all_celltypes"}' |
    sort -u ;
    
    # ABC all celltypes #
    # from cS2G: we used ABC links, kept the maximum ABC score across the 167 cell-types 
    # when an enhancer was interacting with a gene in multiple cell-types, and used this 
    # score as raw linking value.
    zcat data/ABC/ABC_AllPredictions_collapsed.bed.gz |
    bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
    cut -f4,5,9,10 |
    awk '{print $0"\tABC_all_celltypes"}' | 
    sort -u ;
    
    # ABC_enriched_celltypes # 
    # (makes predictions for the whole CS)
    cat $base_dir/ABC-GWAS-Paper/ABC-Max/out/ABC/$variants/tgp_settings/GenePredictions.allCredibleSets.tsv |
    grep -v '#' |
    awk -F'\t' -vOFS='\t' 'NR > 1 && $6 == "TRUE" && $16 != "NA" {print "NA",$2,$9,$16,"ABC_enriched_celltypes"}' |
    # fix gencode incompatible C1orf106 -> INAVA
    sed 's/C1orf106/INAVA/g' |
    sort -u ;
    
    # ABC_MCF7 #
    # zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz | 
    # awk '$24 == "MCF-7-ENCODE"' |
    # sort -k1,1 -k2,2n |
    # bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
    # awk -F'\t' -vOFS='\t' '{print $1,$2,$3,$5,$12,$26,"ABC"}' | 
    # sort -u ;
    
    # EpiMAP_all_celltypes #
    # from cS2G: we used linked EpiMap enhancers based on expression-enhancer activity
    # correlation across 833 cell-types, kept the maximum correlation when an enhancer
    # was linked to a gene in multiple tissues, and used this correlation as raw linking
    # value.
    # must run data/EpiMAP/links/EpiMAP_corrlinks_collapsed.sh to generate collapsed links for intersection
    zcat data/EpiMAP/links/EpiMAP_corrlinks_collapsed.bed.gz |
    bedtools intersect -a - -b <(cut -f1-5 $variants_bed) -wa -wb | 
    awk -F'\t' -vOFS='\t' '{print $9,$10,$4,$5,"EpiMAP_all_celltypes"}' |
    sort -u ;
    
    # # EpiMAP #
    # join -1 9 -2 4 -t$'\t' -o 1.1,1.2,1.3,1.5,2.5,1.10 \
    # <(  zcat data/EpiMAP/links/links_corr_only/${EpiMAP_BSS}_collated_pred.tsv.gz |
    #     bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb | sort -k9,9 ) \
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
  
    # MAGMA # maximum positive Z score
    bedtools intersect \
      -a <(cut -f1-5 $variants_bed) \
      -b <(awk 'NR>1 {print "chr"$2"\t"$3"\t"$4"\t"$0}' data/MAGMA/$variants/magma.genes.out | sort -k1,1 -k2,2n) \
      -wa -wb |
    awk -F'\t' -vOFS='\t' '$16 !~ "-" {print $4,$5,$18,$16,"MAGMA"}' |
    sort -u ;
    
    # INQUISIT # (must convert levels (1,2,3) to scores (3,2,1))
    if [ $variants == 'BC_Michailidou2017_FM' ] ; then
      INQ_dir=/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/data/precision_recall/predictions/INQUISIT/
      ( 
        cat $INQ_dir/inq3_predictions.txt | awk -F'\t' -vOFS='\t' 'NR>1{print "NA",$1,$2,"1","INQUISIT"}' ;
        cat $INQ_dir/inq2_predictions.txt | awk -F'\t' -vOFS='\t' 'NR>1{print "NA",$1,$2,"2","INQUISIT"}' ; 
        cat $INQ_dir/inq1_predictions.txt | awk -F'\t' -vOFS='\t' 'NR>1{print "NA",$1,$2,"3","INQUISIT"}' ;
      ) | cat |
      sort -u ;
    fi
    
  ) | cat > $out_dir/predictions.tsv

done < <(sed 1d $traits_metadata) )

# test performance
Rscript code/compare_genes_2.R

# unite performance tables
awk 'FNR==1 && NR!=1 { while (/^variants/) getline ; } ; 1 {print}' \
output/*/performance.tsv \
> output/performance.tsv

# unite performance plot pdfs
pdfunite output/*/performance.pdf output/performance.pdf

# unite performance plot pngs
(
cd output/ ; for png in */performance.png ; do cp $png ${png%/*}.png ; done
zip performance.png.zip *.png ; rm -f *.png
)