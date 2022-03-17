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

# ( trait="IBD_index_proxies" ; cd /working/lab_jonathb/alexandT/tgp_paper/compare_methods/ ; qsub -v trait=$trait -l mem=50GB,walltime=24:00:00 -e log/ -o log/ code/run.sh )
