#!/bin/bash
cd /working/lab_jonathb/alexandT/tgp_paper/compare_methods/
(while IFS=$'\t' read -r trait variants variants_file known_genes_file ; do
  echo $variants
  qsub -v trait="$variants" -e log/$variants.ER -o log/$variants.OU -l mem=30GB,walltime=6:00:00 code/compare_methods_1.sh
done < <(sed 1d /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/metadata.tsv) )
