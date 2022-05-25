#!/bin/bash
# trait  info passed from run/qsub_run.sh 
# for testing: # trait=BC;variants=BC_cS2GxMichailidou2017_assoc;variants_file=$traits_dir/output/$variants/variants.bed;known_genes_file=$traits_dir/output/$variants/known_genes.txt

module load bedtools/2.27.1
module load R/4.0.2

base_dir=/working/lab_jonathb/alexandT/ 
WKDIR=$base_dir/tgp_paper/compare_methods/ ; cd $WKDIR
traits_dir=$base_dir/tgp_paper/wrangle_package_data/traits/
traits_metadata=$traits_dir/output/metadata.tsv
tissue_matches=output/intermethod_tissue_matches.tsv
ABC_metadata=data/ABC/metadata/CellTypes.Annotated.ABCPaper.txt
EpiMAP_metadata=data/EpiMAP/metadata/main_metadata_table.tsv

# function to get tgp-enriched tissues, for method matching... #
get_tgp_tissues () {
  method=$1
  method_col=$( head -1 $tissue_matches | tr '\t' '\n' | grep -n $method | cut -f1 -d: )
  join -j1 -t$'\t' \
    <(  sort -k1,1 $tissue_matches ) \
    <(  cat $base_dir/tgp/out/$variants/enriched_tissues/tissue_enrichments.tsv |
        awk -F'\t' -vOFS='\t' 'NR > 1 && $10 == "TRUE" {print $2}' |
        sort -k1,1 -u )  |
  cut -f $method_col | sed -z '$ s/\n$//' | tr '\n' '|' |
  grep -E -f - data/$method/metadata/* |
  cut -f1 | sed -z '$ s/\n$//' | tr '\n' '|'
}

echo $variants ; >&2 echo $variants
out_dir=output/$variants/ ; mkdir -p $out_dir
variants_bed=$traits_dir/output/$variants/variants.bed

# predictions.tsv #
(
  echo -e "variant\tcs\tsymbol\tscore\tmethod" ;

  >&2 echo "> tgp_enriched_tissues"
  cat $base_dir/tgp/out/$variants/enriched_tissues/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_enriched_tissues"}' |
  sort -u ;

  >&2 echo "> tgp_enriched_celltypes"
  cat $base_dir/tgp/out/$variants/enriched_celltypes/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_enriched_celltypes"}' |
  sort -u ;

  >&2 echo "> tgp_all_celltypes"
  cat $base_dir/tgp/out/$variants/all_celltypes/target_gene_predictions_full.tsv |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $7 > 0 {print $5,$1,$6,$7,"tgp_all_celltypes"}' |
  sort -u ;

  >&2 echo "> ABC_all_celltypes"
  # from cS2G: we used ABC links, kept the maximum ABC score across the 167 cell-types
  # when an enhancer was interacting with a gene in multiple cell-types, and used this
  # score as raw linking value.
  zcat data/ABC/ABC_AllPredictions_collapsed.bed.gz |
  bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
  awk  -F'\t' -vOFS='\t' '{print $4,$5,$9,$10,"ABC_all_celltypes"}' |
  sort -u ;

  >&2 echo "> ABC_enriched_celltypes"
  # (makes predictions for the whole CS)
  cat $base_dir/ABC-GWAS-Paper/ABC-Max/out/ABC/$variants/tgp_settings/GenePredictions.allCredibleSets.tsv |
  grep -v '#' |
  awk -F'\t' -vOFS='\t' 'NR > 1 && $6 == "TRUE" && $16 != "NA" {print "NA",$2,$9,$16,"ABC_enriched_celltypes"}' |
  # fix gencode incompatible C1orf106 -> INAVA
  sed 's/C1orf106/INAVA/g' |
  sort -u ;

  if [[ $(get_tgp_tissues ABC) != '' ]] ; then
    >&2 echo "> ABC_tgp-enriched_tissues"
    zcat data/ABC/ABC_AllPredictions_all.bed.gz |
    grep -E $(get_tgp_tissues ABC) |
    bedtools intersect  -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
    awk  -F'\t' -vOFS='\t' '{print $4,$5,$9,$10,"ABC_tgp-enriched_tissues"}' |
    sort -u ;
  fi

  if [[ $variants == 'BC_Michailidou2017_FM' ]] ; then
    >&2 echo "> ABC_MCF7"
    zcat data/ABC/ABC_AllPredictions_all.bed.gz |
    awk '$6 == "MCF-7-ENCODE"' |
    sort -k1,1 -k2,2n |
    bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
    awk  -F'\t' -vOFS='\t' '{print $4,$5,$9,$10,"ABC_MCF7"}' |
    sort -u ;
  fi

  >&2 echo "> EpiMAP_all_celltypes"
  # from cS2G: we used linked EpiMap enhancers based on expression-enhancer activity
  # correlation across 833 cell-types, kept the maximum correlation when an enhancer
  # was linked to a gene in multiple tissues, and used this correlation as raw linking
  # value.
  # must run data/EpiMAP/links/EpiMAP_corrlinks_collapsed.sh to generate collapsed links for intersection
  zcat data/EpiMAP/links/EpiMAP_corrlinks_collapsed.bed.gz |
  bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
  awk -F'\t' -vOFS='\t' '{print $4,$5,$9,$10,"EpiMAP_all_celltypes"}' |
  sort -u ;

  if [[ $(get_tgp_tissues EpiMAP) != '' ]] ; then
    >&2 echo "> EpiMAP_tgp-enriched_tissues"
    zcat data/EpiMAP/links/EpiMAP_corrlinks_all.bed.gz |
    grep -E $(get_tgp_tissues EpiMAP)  |
    bedtools intersect  -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
    awk -F'\t' -vOFS='\t' '{print $4,$5,$9,$10,"EpiMAP_tgp-enriched_tissues"}' |
    sort -u ;
  fi

  if [[ $variants == 'BC_Michailidou2017_FM' ]] ; then
    >&2 echo "> EpiMAP_MCF7"
    awk '$3 ~ "MCF-7" {print $1}' data/EpiMAP/metadata/main_metadata_table.tsv |
              sed -z '$ s/\n$//' | tr '\n' '|' |
    zgrep -E -f -  data/EpiMAP/links/EpiMAP_corrlinks_all.bed.gz  |
    sort -k1,1 -k2,2n |
    bedtools intersect -a <(cut -f1-5 $variants_bed) -b - -wa -wb |
    awk  -F'\t' -vOFS='\t' '{print $4,$5,$9,$10,"EpiMAP_MCF7"}' |
    sort -u ;
  fi

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

  >&2 echo "> cS2G_found_variants"
  # liftover
  # . data/cS2G/liftover.sh
  # variant-gene pair scores are the same across traits
  # variants column must be rsIDs
  join -t$'\t' -1 4 -2 1 -o 1.4,1.5,2.2,2.3 \
      <(sort -k4,4 $traits_dir/output/$variants/variants.bed) \
      <(zcat data/cS2G/cS2G_1000GEUR/cS2G.*.SGscore.gz | sort -k1,1) |
  awk -F'\t' -vOFS='\t' '{print $0,"cS2G_found_variants"}' |
  sort -u ;
  
  if [[ $variants =~ '_cS2Gx' ]] ; then
    PMID=$(cat $traits_dir/output/metadata.tsv | awk -v variants=$variants '$2==variants {print $5}')
    if [[ $PMID != '' ]] ; then
      >&2 echo "> cS2G_study_associated_variants"
      cat $base_dir/tgp_paper/compare_methods/data/cS2G/gwas_catalog_cS2G/gwas_catalog_v1.0-associations_e100_r2021-02-25.annot.hg19.bed |
      awk -F'\t' -vOFS='\t' -v PMID=$PMID '$4 == PMID {print $8,$8,$9,$10,"cS2G_study_associated_variants"}' |
      sort -u ;
    fi
  fi
    
  MAGMA=data/MAGMA/$variants/magma.genes.out 
  if [[ -f $MAGMA ]] ; then 
    >&2 echo "> MAGMA"
    # p value cutoff = 0.05 / n_genes
    n_genes=$(cat $MAGMA | sed 1d | cut -f1 | uniq | wc -l)
    max_p_val=$(echo "scale = 10; 0.05 / $n_genes" | bc)
    # maximum positive Z score
    bedtools intersect \
      -a <( cut -f1-5 $variants_bed ) \
      -b <( cat $MAGMA |
            awk -F'\t' -vOFS='\t'  \
              -v max_p_val=$max_p_val \
              'NR>1 && ($9 + 0)<max_p_val {print "chr"$2,$3,$4,$0}' |
              sort -k1,1 -k2,2n ) \
      -wa -wb |
    awk -F'\t' -vOFS='\t' '$16 !~ "-" {print $4,$5,$18,$16,"MAGMA"}' |
    sort -u ;
  fi

  if [[ $variants == 'BC_Michailidou2017_FM' ]] ; then
    >&2 echo "> MAGMA"
    # p value cutoff = 0.05 / n_genes
    n_genes=$(cat data/MAGMA/$variants/magma.genes.out | sed 1d | cut -f1 | uniq | wc -l)
    max_p_val=$(echo "scale = 10; 0.05 / $n_genes" | bc)
    # maximum positive Z score
    bedtools intersect \
      -a <( cut -f1-5 $variants_bed ) \
      -b <( cat data/MAGMA/$variants/magma.genes.out |
            awk -F'\t' -vOFS='\t'  \
              -v max_p_val=$max_p_val \
              'NR>1 && ($9 + 0)<max_p_val {print "chr"$2,$3,$4,$0}' |
              sort -k1,1 -k2,2n ) \
      -wa -wb |
    awk -F'\t' -vOFS='\t' '$16 !~ "-" {print $4,$5,$18,$16,"MAGMA"}' |
    sort -u ;
  fi

  if [[ $variants == 'BC_Michailidou2017_FM' ]] ; then
    >&2 echo "> INQUISIT"
    # (must convert levels (1,2,3) to scores (3,2,1))
    INQ_dir=/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/data/precision_recall/predictions/INQUISIT/
    (
      cat $INQ_dir/inq3_predictions.txt | awk -F'\t' -vOFS='\t' 'NR>1 {print "NA",$1,$2,"1","INQUISIT"}' ;
      cat $INQ_dir/inq2_predictions.txt | awk -F'\t' -vOFS='\t' 'NR>1 {print "NA",$1,$2,"2","INQUISIT"}' ;
      cat $INQ_dir/inq1_predictions.txt | awk -F'\t' -vOFS='\t' 'NR>1 {print "NA",$1,$2,"3","INQUISIT"}' ;
    ) | cat |
    # replace ichav
    sed 's/ichav/\./g' |
    sort -u ;
  fi

) | cat > $out_dir/predictions.tsv

# test performance
Rscript code/compare_genes_2.R $variants

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