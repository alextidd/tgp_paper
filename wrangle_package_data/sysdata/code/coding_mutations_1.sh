#!/bin/bash
module load R/4.0.2
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/sysdata/ ; cd $WKDIR
out_dir=${WKDIR}/output/coding_mutations/ ; mkdir -p $out_dir

# dbNSFP data
# v4.1a: col79 is the REVEL score, col8/9 is hg19 coordinates, col13 is genename, col14 is ENSGs, col15 is ENSTs
dbNSFPv=4.1a 
ENSG_col=14 
score_col=79
dbNSFPdata=/reference/data/dbNSFP/$dbNSFPv/

# missense ===
( echo -e "chrom\tposition\tensgs\tscore" ; 
  awk -F'\t' -vOFS='\t' -v ENSG_col=$ENSG_col -v score_col=$score_col \
  '$0 !~ "#" && $7 ~ /^rs/ && $score_col >= 0.2 && $8!="." { print "chr"$8,$9,$ENSG_col,$score_col }'  \
  ${dbNSFPdata}/dbNSFP${dbNSFPv}_variant.chr* ; ) |
cat > $out_dir/missense_SNVs.tsv

# nonsense ===
( echo -e "chrom\tposition\tensgs" ; 
  awk -F'\t' -vOFS='\t' -v ENSG_col=$ENSG_col \
  '$0 !~ "#" && $7 ~ /^rs/ && $6 == "X" && $8!="." { print "chr"$8,$9,$ENSG_col }' \
  ${dbNSFPdata}/dbNSFP${dbNSFPv}_variant.chr* ; ) |
cat > $out_dir/nonsense_SNVs.tsv

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
echo -e "chrom\tposition\tensgs" > $out_dir/splicesite_SNVs.tsv
(for i in chr{1..22} chrX ; do
  echo $i

  # recommended adaboost and rf cutoffs = 0.6
  awk -F'\t' -vOFS='\t' '$17 >= 0.6 || $18 >= 0.6 { print $1,$2,$3,$4,$10,$14,$16 }' ${dbscSNVdata}/dbscSNV1.1.${i} |

  # extract position from dbSNP so we can sort the SNPs from other mutations
  awk -F'\t' -vOFS='\t' 'NR==FNR{ a[$3]=$4FS$5FS$6 ; next }{ print "chr"$0,a[$2] }' <( zcat ${dbSNP147}/snp147.alleles.${i}.gz ) - |

  # match alleles
  awk -F'\t' -vOFS='\t' '$3==$9 && $4==$10 && NR>1 { print $1,$2,$6 }'  |

  # extract ENSGs
  perl -pe 's/\(.*?\)//g' \
  >> $out_dir/splicesite_SNVs.tsv
done)

# fix ENSGs column
Rscript code/coding_mutations_2.R