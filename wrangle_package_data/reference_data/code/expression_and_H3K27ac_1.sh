#!/bin/bash
# job submission: # ( cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; qsub -l ncpus=20,mem=196G,walltime=6:00:00 -e log/ -o log/ code/expression_and_H3K27ac_1.sh )

baseDir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $baseDir
WD=$baseDir/data/expression_and_H3K27ac/
OD=$baseDir/output/ 
mkdir -p $WD $OD/{expression,H3K27ac}
module load R/3.6.2

# # copy package #
# cp -R \
# /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/signal_matrix/build_matrix/build_matrix/* \
# $WD

# 1. expression processing #

# generate expression matrices
Rscript code/expression_data_processing.R 

# 2. H3K27ac processing #

# extract bw signal per DHS site using multibigwigsummary
. code/make_H3K27ac_matrix.sh

# reformat raw matrix, normalise, and export
Rscript code/h3k27ac_data_processing.R
