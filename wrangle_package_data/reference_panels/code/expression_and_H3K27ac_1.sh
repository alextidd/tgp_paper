#!/bin/bash
# job submission: # ( cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/ ; qsub -l ncpus=20,mem=196G,walltime=6:00:00 -e log/ -o log/ code/expression_and_H3K27ac_1.sh )
cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/

##### H3K27ac processing #####
module load python/2.7.13

# create output dir
mkdir output/H3K27ac/

# compute signal (only for celltypes in metadata)
accessions=$(cat output/metadata.tsv | awk '$3=="H3K27ac" {print $6}' | tr '\n' '\|') ; accessions=${accessions%??}
multiBigwigSummary BED-file \
  -b $(find data/H3K27ac/*.big*ig | grep -E "$accessions") \
  --BED data/DHS.bed \
  -o output/H3K27ac/H3K27ac_multiBigwigSummary.npz \
  -p 20 \
  --outRawCounts output/H3K27ac/H3K27ac_multiBigwigSummary.tmp

# extract header
head -n1 output/H3K27ac/H3K27ac_multiBigwigSummary.tmp > output/H3K27ac/H3K27ac_multiBigwigSummary.header.tmp

# sort and remove duplicate rows
sed '1d' output/H3K27ac/H3K27ac_multiBigwigSummary.tmp |
sort |
uniq > output/H3K27ac/H3K27ac_multiBigwigSummary.sorted.tmp

# add header
cat output/H3K27ac/H3K27ac_multiBigwigSummary.header.tmp output/H3K27ac/H3K27ac_multiBigwigSummary.sorted.tmp > output/H3K27ac/H3K27ac_input_matrix.tsv

# remove temp
rm -f output/H3K27ac/*.tmp output/H3K27ac/H3K27ac_multiBigwigSummary.npz

#### H3K27ac and expression processing #####
module load R/4.0.2
Rscript code/expression_and_H3K27ac_2.R
