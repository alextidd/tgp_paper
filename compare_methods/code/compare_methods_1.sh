#!/bin/bash
# script for loading the predictions made by methods of interest and intersecting with trait variants
# RUNS: 3 traits x 3 celltype groupings ( x 3 methods )
trait=$1

module load R/4.0.2
module load bedtools/2.27.1
module load ucsctools/20160223
baseDir=/working/lab_jonathb/alexandT/ ; 
WKDIR=$baseDir/tgp_paper/compare_methods/ ; cd $WKDIR
celltypes=enriched_tissues # (enriched_tissues BRST.MCF7.CNCR_celltype all_celltypes)

echo $trait $celltypes
outDir=output/$trait/$celltypes/ ; mkdir -p $outDir
TraitVarsBed=$baseDir/tgp_paper/wrangle_package_data/traits/output/$trait/variants.bed

# header
echo -e "cs\tsymbol\tscore\tmethod" | gzip > $outDir/predictions_long.tsv.gz

# tgp predictions #
cat $baseDir/tgp/out/$trait/$celltypes/target_gene_annotations.tsv | 
sed 1d |
cut -f1,7,8 |
awk '{print $0"\ttgp"}' |
gzip >> $outDir/predictions_long.tsv.gz

# cS2G predictions #
(if [ "$trait" == "BC" ] && [ "$celltypes" = "enriched_tissues" ] ; then 
  
fi)

# ABC predictions #
(if [ "$celltypes" = "BRST.MCF7.CNCR_celltype" ] ; then
  zcat data/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz |
  sed 1d | 
  grep MCF-7 | 
  sort -k1,1 -k2,2n | 
  bedtools intersect -a - -b $TraitVarsBed -wa -wb |
  awk '{print $29"\t"$7"\t"$21"\tABC"}' |
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $outDir/predictions_long.tsv.gz
elif [ "$celltypes" = "enriched_tissues" ] ; then 
  cat $baseDir/ABC-GWAS-Paper/ABC-Max/out/ABC/$trait/tgp_settings/GenePredictions.allCredibleSets.tsv |
  grep -v '#' |
  awk 'NR > 1 && $6 == "TRUE" && $16 != "NA" {print $2"\t"$9"\t"$16"\tABC"}' | 
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  gzip >> $outDir/predictions_long.tsv.gz
elif [ "$celltypes" = "all_celltypes" ] ; then 
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
if [ "$trait" == "PrCa_Giambartolomei2021" ] || [ "$trait" == "BC" ] ; then 
  (if [ "$celltypes" = "BRST.MCF7.CNCR_celltype" ] ; then
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
  elif [ "$celltypes" = "enriched_tissues" ] ; then 
    # run: # code/get_EpiMAP.sh $trait
    
    PMID="29059683" ; trait_info="breast cancer"
    Trait_info=${trait_info^}
    uid="$PMID - $Trait_info"
    # to get links: filter to GWAS paper, scored links, distance to center < 2500
    zcat data/EpiMAP/gwas_resources/all_gwas_SNP_links.tsv.gz |
    awk -F'\t' -vOFS='\t' -v uid="$uid" \
    '$13==uid && $7 != "NA" && $4 < 2500'  
    
    
    cat data/EpiMAP/GWAS_enrichments/$trait/corrlinks_inloci.tsv | 
    awk -F'\t' -vOFS='\t' 'NR > 1 && $6 != "NA" {print $1,$2-1,$2,$0}' |  
    bedtools intersect -a - -b $baseDir/tgp_paper/wrangle_package_data/traits/output/$trait/CredibleSetList.bed -wa -wb | 
    awk -F'\t' -vOFS='\t' '{print $19,$9,$10,"EpiMAP"}' |
    gzip >> $outDir/predictions_long.tsv.gz
  elif [ "$celltypes" = "all_celltypes" ] ; then 
    echo -e "cs\tensg\tscore" | gzip > $outDir/EpiMAP/intersected_predictions.tsv.gz
    ls data/EpiMAP/predictions/*_collated_pred.tsv.gz |
    xargs zcat |
    bedtools intersect -a - -b $TraitVarsBed -wa -wb |
    awk -F'\t' -vOFS='\t' '{print $11,$4,$5}' |
    gzip >> $outDir/EpiMAP/intersected_predictions.tsv.gz
    # convert to symbols
    Rscript code/ensg_to_symbol.R $outDir/EpiMAP/intersected_predictions.tsv.gz FALSE
    # add to predictions_long
    cat $outDir/EpiMAP/intersected_predictions_with_symbols.tsv | 
    sed 1d |
    awk '{print $0"\tEpiMAP"}' |
    gzip >> $outDir/predictions_long.tsv.gz
  fi)
fi

# compare methods
Rscript code/compare_methods_2.R $trait $celltypes

# unite performance tables
awk 'FNR==1 && NR!=1 { while (/^trait/) getline ; } ; 1 {print}' \
output/*/enriched_tissues/performance.tsv \
> output/performance.tsv

# unite performance plots
pdfunite output/*/enriched_tissues/performance.pdf \
  output/performance.pdf