#!/bin/bash
module load R/4.0.2
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
REVELOutDir=${WKDIR}/output/REVEL/ ; mkdir -p $REVELOutDir

# dbNSFP data
# v4.1a: col79 is the REVEL score, col8/9 is hg19 coordinates, col13 is genename, col14 is ENSGs, col15 is ENSTs
dbNSFPv=4.1a 
ENSG_col=14 
score_col=79
dbNSFPdata=/reference/data/dbNSFP/${dbNSFPv}/

# missense ===
awk -F'\t' -v ENSG_col=$ENSG_col -v score_col=$score_col \
'$0 !~ "#" && $7 ~ /^rs/ && $score_col >= 0.2 && $8!="." { print "chr"$8,$9,$ENSG_col,$score_col }' OFS='\t' \
${dbNSFPdata}/dbNSFP${dbNSFPv}_variant.chr* \
> ${REVELOutDir}/missense_SNVs.tsv

# nonsense ===
awk -F'\t' -v ENSG_col=$ENSG_col \
'$0 !~ "#" && $7 ~ /^rs/ && $6 == "X" && $8!="." { print "chr"$8,$9,$ENSG_col }' OFS='\t' \
${dbNSFPdata}/dbNSFP${dbNSFPv}_variant.chr* \
> ${REVELOutDir}/nonsense_SNVs.tsv

# splicesite ===
# DB = dbscSNV
# http://www.liulab.science/dbscsnv.html
# dbscSNV includes all potential human SNVs within splicing consensus regions (−3 to +8 at the 5’ splice site and −12 to +2 at the 3’ splice site),
# i.e. scSNVs, related functional annotations and two ensemble prediction scores for predicting their potential of altering splicing.
# PROBLEM - does not contain rsIDs - use dbSNP to extract the correct mutation

# data
dbscSNVdata=/reference/data/dbscSNV/1.1/
dbSNP147=/working/lab_georgiat/jonathB/PROJECTS/general/snp_ids/data/link/snp147_byChr_alleles/

# loop per chrom
> ${REVELOutDir}/splicesite_SNVs.tsv
(for i in chr{1..22} chrX ; do
  echo $i

  # recommended adaboost and rf cutoffs = 0.6
  awk -F'\t' '$17 >= 0.6 || $18 >= 0.6 { print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$14"\t"$16 }' ${dbscSNVdata}/dbscSNV1.1.${i} |

  # extract position from dbSNP so we can sort the SNPs from other mutations
  awk -F'\t' 'NR==FNR{ a[$3]=$4FS$5FS$6 ; next }{ print "chr"$0"\t"a[$2] }' <( zcat ${dbSNP147}/snp147.alleles.${i}.gz ) - |

  # match alleles
  awk -F'\t' '$3==$9 && $4==$10 && NR>1 { print $1"\t"$2"\t"$6 }'  |

  # extract ENSGs
  perl -pe 's/\(.*?\)//g' \
  >> ${REVELOutDir}/splicesite_SNVs.tsv

done)

# fix ENSGs column
Rscript code/REVEL_2.R