#!/bin/bash
baseDir=/working/lab_jonathb/alexandT/ ; cd $baseDir

# run ABC
( cd ABC-GWAS-Paper/ABC-Max/
  run/run.sh ABC $trait NA tgp_settings
)

# run tgp
module load R/4.0.2
( cd tgp/ 
  Rscript run/run.R $trait
)

# compare methods
( cd tgp_paper/compare_methods/
  code/compare_methods_1.sh $trait
)

# SUBMIT ONE TRAIT: #
# # ( trait="CRC_index_proxies" ; cd /working/lab_jonathb/alexandT/tgp_paper/compare_methods/ ; qsub -v trait=$trait -l mem=50GB,walltime=24:00:00 -e log/ -o log/ code/run.sh )
# SUBMIT ALL AVAILABLE TRAITS: #
# # ( cd /working/lab_jonathb/alexandT/tgp_paper/compare_methods/ ; for path in $(ls ../wrangle_package_data/output/Traits/*/VariantList.bed) ; do trait=$(echo $path | sed 's/.*\/Traits\///g' | sed 's/\/.*//g') ; qsub -v trait=$trait -l mem=50GB,walltime=24:00:00 -N $trait.compare -e log/$trait.ER -o log/$trait.OU code/run.sh ; done)
