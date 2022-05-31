#!/bin/bash
# job submission: # ( cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/ ; qsub -l ncpus=20,mem=196G,walltime=6:00:00 -e log/ -o log/ code/expression_and_H3K27ac_1.sh )
cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/

##### expression inputs #####
# # hpcapp01 head node: #
# ( 
#   cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/expression
#   tawk=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/code/t.awk
#   while read accession ; do echo $accession
#     wget https://www.encodeproject.org/files/$accession/@@download/$accession.tsv -O $accession.tmp
#     # get gene_id and TPM, remove trailing tab
#     ( echo -e "gene_id\tTPM" ; 
#       awk -F'\t' -f $tawk -v cols=gene_id,TPM  $accession.tmp | sed 's/\t$//g' ; 
#     ) | cat > $accession.tsv
#     rm -f $accession.tmp
#   done < <(grep ENCFF README.txt) )
#   
#   wget \
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73784/suppl/GSE73784_LNCaP.TPM.expression.data.tsv.gz \
#   -O GSE73784.tsv.tmp.gz
#   gunzip GSE73784.tsv.tmp.gz
#   cut -f1,9 GSE73784.tsv.tmp > GSE73784.tsv
#   rm -f GSE73784.tsv.tmp
#   
#   cat /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/RNA-seq/HGSOC/hgsoc.tpm.txt |
#   cut -f1,4 \
#   > GSE121103.tsv
# )

##### H3K27ac inputs #####
# for file in /working/lab_georgiat/jonathB/genomic_data_downloads/ROADMAP/ChIP-seq/bigwig/H3K27ac/* ; do
#   cp \
#     $file \
#     /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/H3K27ac/$(basename ${file%-*}).bigWig
# done
# cp \
# /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/ChIP-seq/SRX4862819.bw \
# /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/H3K27ac/SRX4862819.bigWig
# cp \
# /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/epi_data_formatting/output/ChIP-seq/signal_tracks/SRR8030178.rpgccoverage.bw \
# /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/H3K27ac/SRR8030178.bigWig
# cp \
# /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/epi_data_formatting/output/ChIP-seq/signal_tracks/PRJNA278653_ishikawa/SRR1917223.rpgccoverage.bw \
# /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/H3K27ac/SRR1917223.bigWig
# cp \
# /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/epi_data_formatting/output/ChIP-seq/ENCFF990CJG.bigWig \
# /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/H3K27ac/ENCFF990CJG.bigWig

##### H3K27ac pre-processing #####
module load python/2.7.13

# create output dir
mkdir output/H3K27ac/

# compute signal (only for celltypes in metadata)
accessions=$(cat output/metadata.tsv | awk '$3=="H3K27ac" {print $5}' | tr '\n' '\|') ; accessions=${accessions%?}
multiBigwigSummary BED-file \
  -b $(find data/H3K27ac/*.bigWig | grep -E "$accessions") \
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
