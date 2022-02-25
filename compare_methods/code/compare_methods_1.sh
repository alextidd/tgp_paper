#!/bin/bash
# script for loading the predictions made by methods of interest and intersecting with trait variants
# RUNS: 3 traits x 3 celltype groupings ( x 3 methods )
module load R/4.0.2
module load bedtools/2.27.1
TRAIT=BC # (BC IBD PrCa) 
CELLTYPES=enriched_tissues # (enriched_tissues BRST.MCF7.CNCR_celltype all_celltypes)

echo $TRAIT $CELLTYPES # CELLTYPES=enriched_tissues 

baseDir=/working/lab_jonathb/alexandT/ ; 
WKDIR=$baseDir/tgp_paper/compare_methods/ ; cd $WKDIR
outDir=output/$TRAIT/$CELLTYPES/ ; mkdir -p $outDir
TraitVarsBed=$baseDir/tgp_paper/wrangle_package_data/output/Traits/$TRAIT/$TRAIT.VariantList.bed

echo -e "cs\tsymbol\tscore\tmethod" | gzip > $outDir/predictions_long.tsv.gz

# tgp predictions #
cat $baseDir/tgp/out/$TRAIT/$CELLTYPES/target_gene_annotations.tsv | 
sed 1d |
cut -f5,7,8 |
awk '{print $0"\ttgp"}' |
gzip >> $outDir/predictions_long.tsv.gz

# ABC predictions #
(if [ "$CELLTYPES" = "BRST.MCF7.CNCR_celltype" ] ; then
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz |
  sed 1d | 
  grep MCF-7 | 
  sort -k1,1 -k2,2n | 
  bedtools intersect -a - -b $TraitVarsBed -wa -wb |
  awk '{print $29"\t"$7"\t"$21"\tABC"}' |
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $outDir/predictions_long.tsv.gz
elif [ "$CELLTYPES" = "enriched_tissues" ] ; then 
  cat $baseDir/ABC-GWAS-Paper/ABC-Max/out/ABC/$TRAIT/tgp_settings/GenePredictions.allCredibleSets.tsv |
  grep -v '#' |
  awk 'NR > 1 && $6 == "TRUE" && $16 != "NA" {print $2"\t"$9"\t"$16"\tABC"}' | 
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $outDir/predictions_long.tsv.gz
elif [ "$CELLTYPES" = "all_celltypes" ] ; then 
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz |
  sed 1d | 
  sort -k1,1 -k2,2n | 
  bedtools intersect -a - -b $TraitVarsBed -wa -wb |
  awk '{print $29"\t"$7"\t"$21"\tABC"}' |
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $outDir/predictions_long.tsv.gz
fi)

# EpiMAP predictions #
mkdir $outDir/EpiMAP/
(if [ "$CELLTYPES" = "BRST.MCF7.CNCR_celltype" ] ; then
  echo -e "cs\tensg\tscore" | gzip > $outDir/EpiMAP/intersected_predictions.tsv.gz
  cat data/EpiMAP/main_metadata_table.tsv | grep MCF-7 | cut -f1 |
  grep -f - <(ls data/EpiMAP/predictions/*_collated_pred.tsv.gz) |
  xargs zcat |
  bedtools intersect -a - -b $TraitVarsBed -wa -wb |
  awk '{print $11"\t"$4"\t"$5}' |
  gzip >> $outDir/EpiMAP/intersected_predictions.tsv.gz
  # convert to symbols
  Rscript code/ensg_to_symbol.R $outDir/EpiMAP/intersected_predictions.tsv.gz FALSE
  # add to predictions_long
  cat $outDir/EpiMAP/intersected_predictions_with_symbols.tsv | 
  sed 1d |
  awk '{print $0"\tEpiMAP"}' |
  gzip >> $outDir/predictions_long.tsv.gz
elif [ "$CELLTYPES" = "enriched_tissues" ] ; then 
  # run: # code/get_EpiMAP.sh $TRAIT
  # zcat data/EpiMAP/GWAS_enrichments/all_gwas_SNP_links.tsv.gz  |     
  # awk -F'\t' 'NR > 1 && $6 != "NA" && $13 == "29059683 - Breast cancer" {print $1"\t"$2-1"\t"$2"\t"$0}' |
  # bedtools intersect -a - -b $baseDir/tgp_paper/wrangle_package_data/output/Traits/BC/BC.CredibleSetList.bed -wa -wb
  
  cat data/EpiMAP/GWAS_enrichments/$TRAIT/corrlinks_inloci.tsv | 
  awk -F'\t' 'NR > 1 && $6 != "NA" {print $1"\t"$2-1"\t"$2"\t"$0}' |  
  bedtools intersect -a - -b $baseDir/tgp_paper/wrangle_package_data/output/Traits/$TRAIT/$TRAIT.CredibleSetList.bed -wa -wb | 
  awk -F'\t' '{print $19"\t"$9"\t"$10"\tEpiMAP"}' |
  gzip >> $outDir/predictions_long.tsv.gz
elif [ "$CELLTYPES" = "all_celltypes" ] ; then 
  echo -e "cs\tensg\tscore" | gzip > $outDir/EpiMAP/intersected_predictions.tsv.gz
  ls data/EpiMAP/predictions/*_collated_pred.tsv.gz |
  xargs zcat |
  bedtools intersect -a - -b $TraitVarsBed -wa -wb |
  awk '{print $11"\t"$4"\t"$5}' |
  gzip >> $outDir/EpiMAP/intersected_predictions.tsv.gz
  # convert to symbols
  Rscript code/ensg_to_symbol.R $outDir/EpiMAP/intersected_predictions.tsv.gz FALSE
  # add to predictions_long
  cat $outDir/EpiMAP/intersected_predictions_with_symbols.tsv | 
  sed 1d |
  awk '{print $0"\tEpiMAP"}' |
  gzip >> $outDir/predictions_long.tsv.gz
fi)
  
# compare methods
Rscript code/compare_methods.R $TRAIT $CELLTYPES




# cat data/EpiMAP/GWAS_enrichments/BC_Michailidou2017/corrlinks_inloci.tsv |
#     sed 1d |
#     awk '{print $1"\t"$2-1"\t"$2+1"\t"$0}' |
#     bedtools intersect -a - \
#     -b /working/lab_georgiat/jonathB/bcac_oncoarray/BCAC.secondarySignalSNPs.9col.bed  \
#     -wa -wb 