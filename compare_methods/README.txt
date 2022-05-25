# output/ structure: #

output/
-$TRAIT/            # (BC, IBD, PrCa) = trait-specific predictions 
--predictions.tsv  # all predictions by all methods (cols: variant, cs, symbol, score, method)
--performance.pdf  # PR, AUPRC
--performance.tsv  # Precision-Recall and summary statistics

# $CELLTYPES
This denotes the group of celltypes from which predictions have been made, which must be consistent across methods for a comparison. 

## enriched_tissues = predictions from tissues in which variants have been enriched in features
tgp = tissues with significant Fisher enrichment (p<0.05, estimate>2) of variants in the top 10% most celltype-specifically H3K27ac-marked DHSs
EpiMAP = tissues with significant Fisher enrichment (p<0.05, estimate>2) of variants in the top 10% most celltype-specifically H3K27ac-marked DHSs (tgp post-processing)
ABC = performed binomial test comparing fraction at which variants overlap ABC enhancers with the fraction at which all common variants overlap ABC enhancers in that celltype (Bonferroni-corrected binomial P value < 0.001)

## ${CELLTYPE}_celltype = predictions from/for one celltype

## all_celltypes = predictions from all celltypes considered by the method