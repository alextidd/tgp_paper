#!/bin/bash
# script for loading the predictions made by methods of interest and intersecting with trait variants
# RUNS: 3 traits x 3 celltype groupings ( x 3 methods )
# $trait

module load R/4.0.2
module load bedtools/2.27.1
module load ucsctools/20160223

echo $trait
celltypes=enriched_tissues # (enriched_tissues BRST.MCF7.CNCR_celltype all_celltypes)

base_dir=/working/lab_jonathb/alexandT/ ; 
WKDIR=$base_dir/tgp_paper/compare_methods/ ; cd $WKDIR
out_dir=output/$trait/$celltypes/ ; mkdir -p $out_dir

variants_bed=$base_dir/tgp_paper/wrangle_package_data/traits/output/$trait/variants.bed
traits_metadata=$base_dir/tgp_paper/wrangle_package_data/traits/output/metadata.tsv
PMID=$(cat $traits_metadata | awk -v trait="$trait" '$2 == trait {print $5}')

# Tier 1 = predicted genes for the trait
# Tier 2 = scored predicted genes for the trait
# Tier 3 = scored predicted variant-gene pairs for the trait
echo -e "cs\tsymbol\tscore\tmethod" | gzip > $out_dir/predictions_long.tsv.gz

# tgp predictions #
cat $base_dir/tgp/out/$trait/$celltypes/target_gene_annotations.tsv | 
awk -F'\t' -vOFS='\t' 'NR > 1 {print $1,$7,$8,"tgp"}' |
gzip >> $out_dir/predictions_long.tsv.gz

# cS2G predictions #
# . data/cS2G/liftover.sh
cat data/cS2G/gwas_catalog_cS2G/gwas_catalog_v1.0-associations_e100_r2021-02-25.annot.hg19.bed |
awk -v PMID="$PMID" '$4 == PMID' |
sort -k1,1 -k2,2n |
bedtools intersect -a - -b $variants_bed -wa -wb |
awk -F'\t' -vOFS='\t' '{print $16,$9,$10,"cS2G"}' |
gzip >> $out_dir/predictions_long.tsv.gz

# ABC predictions #
(if [ "$celltypes" = "BRST.MCF7.CNCR_celltype" ] ; then
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz |
  sed 1d | 
  grep MCF-7 | 
  sort -k1,1 -k2,2n | 
  bedtools intersect -a - -b $variants_bed -wa -wb |
  awk '{print $29"\t"$7"\t"$21"\tABC"}' |
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $out_dir/predictions_long.tsv.gz
elif [ "$celltypes" = "enriched_tissues" ] ; then 
  cat $base_dir/ABC-GWAS-Paper/ABC-Max/out/ABC/$trait/tgp_settings/GenePredictions.allCredibleSets.tsv |
  grep -v '#' |
  awk 'NR > 1 && $6 == "TRUE" && $16 != "NA" {print $2"\t"$9"\t"$16"\tABC"}' | 
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $out_dir/predictions_long.tsv.gz
elif [ "$celltypes" = "all_celltypes" ] ; then 
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz |
  sed 1d | 
  sort -k1,1 -k2,2n | 
  bedtools intersect -a - -b $variants_bed -wa -wb |
  awk '{print $29"\t"$7"\t"$21"\tABC"}' |
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $out_dir/predictions_long.tsv.gz
fi)

# EpiMAP predictions #
# to get corrlinks: filter to GWAS paper, scored links only, distance to center < 2500
zcat data/EpiMAP/gwas_resources/all_gwas_SNP_links.tsv.gz |
awk -F'\t' -vOFS='\t' -v uid_PMID="$PMID - " \
  '$13~uid_PMID && $7 != "NA" && $4 < 2500 {print $1,$2-1,$2,$0}' |
bedtools intersect -a - -b $variants_bed -wa -wb |
awk -F'\t' -vOFS='\t' '{print $23,$9,$10,"EpiMAP"}' |
gzip >> $out_dir/predictions_long.tsv.gz

# compare methods
Rscript code/compare_methods_2.R $trait $celltypes

# unite performance tables
awk 'FNR==1 && NR!=1 { while (/^trait/) getline ; } ; 1 {print}' \
output/*/enriched_tissues/performance.tsv \
> output/performance.tsv

# unite performance plots
pdfunite output/*/enriched_tissues/performance.pdf \
  output/performance.pdf
(
cd output/ ; for png in */*/performance.png ; do cp $png ${png%%/*}.png ; done
zip performance.png.zip *.png ; rm -f *.png 
)
  
  
  
  
  
  
##########

# if [ "$trait" == "PrCa_Giambartolomei2021" ] || [ "$trait" == "BC" ] ; then 
#   (if [ "$celltypes" = "BRST.MCF7.CNCR_celltype" ] ; then
#     echo -e "cs\tensg\tscore" | gzip > $out_dir/EpiMAP/intersected_predictions.tsv.gz
#     cat data/EpiMAP/main_metadata_table.tsv | grep MCF-7 | cut -f1 |
#     grep -f - <(ls data/EpiMAP/predictions/*_collated_pred.tsv.gz) |
#     xargs zcat |
#     bedtools intersect -a - -b $variants_bed -wa -wb |
#     awk '{print $11"\t"$4"\t"$5}' |
#     gzip >> $out_dir/EpiMAP/intersected_predictions.tsv.gz
#     # convert to symbols
#     Rscript code/ensg_to_symbol.R $out_dir/EpiMAP/intersected_predictions.tsv.gz FALSE
#     # add to predictions_long
#     cat $out_dir/EpiMAP/intersected_predictions_with_symbols.tsv | 
#     sed 1d |
#     awk '{print $0"\tEpiMAP"}' |
#     gzip >> $out_dir/predictions_long.tsv.gz
#   elif [ "$celltypes" = "enriched_tissues" ] ; then 
#     # run: # code/get_EpiMAP.sh $trait
#     
    # PMID="29059683" ; trait_info="breast cancer"
    # Trait_info=${trait_info^}
    # uid="$PMID - $Trait_info"
    # # to get links: filter to GWAS paper, scored links, distance to center < 2500
    # zcat data/EpiMAP/gwas_resources/all_gwas_SNP_links.tsv.gz |
    # awk -F'\t' -vOFS='\t' -v uid="$uid" \
    # '$13==uid && $7 != "NA" && $4 < 2500 {print $1,$2-1,$2,$0}' |
    # bedtools intersect -a - -b $variants_bed -wa -wb |
    # awk -F'\t' -vOFS='\t' '{print $23,$9,$10,"EpiMAP"}' |
    # gzip >> $out_dir/predictions_long.tsv.gz
#     
#     # cat data/EpiMAP/GWAS_enrichments/$trait/corrlinks_inloci.tsv | 
#     # awk -F'\t' -vOFS='\t' 'NR > 1 && $6 != "NA" {print $1,$2-1,$2,$0}' |  
#     # bedtools intersect -a - -b $base_dir/tgp_paper/wrangle_package_data/traits/output/$trait/CredibleSetList.bed -wa -wb | 
#     # awk -F'\t' -vOFS='\t' '{print $19,$9,$10,"EpiMAP"}' |
#     # gzip >> $out_dir/predictions_long.tsv.gz
#   elif [ "$celltypes" = "all_celltypes" ] ; then 
#     echo -e "cs\tensg\tscore" | gzip > $out_dir/EpiMAP/intersected_predictions.tsv.gz
#     ls data/EpiMAP/predictions/*_collated_pred.tsv.gz |
#     xargs zcat |
#     bedtools intersect -a - -b $variants_bed -wa -wb |
#     awk -F'\t' -vOFS='\t' '{print $11,$4,$5}' |
#     gzip >> $out_dir/EpiMAP/intersected_predictions.tsv.gz
#     # convert to symbols
#     Rscript code/ensg_to_symbol.R $out_dir/EpiMAP/intersected_predictions.tsv.gz FALSE
#     # add to predictions_long
#     cat $out_dir/EpiMAP/intersected_predictions_with_symbols.tsv | 
#     sed 1d |
#     awk '{print $0"\tEpiMAP"}' |
#     gzip >> $out_dir/predictions_long.tsv.gz
#   fi)
# fi
