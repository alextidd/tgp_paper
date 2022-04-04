#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/ ; cd $WKDIR
inDir=data/GENCODE/ ; mkdir $inDir
outDir=output/GENCODE/ ; mkdir $outDir
pcENSGs=$outDir/proteincoding.gencode.v34lift37.basic.ENSGs.txt

# copy files
sourceDir=/working/lab_georgiat/jonathB/DATA/gene_annotations/GENCODE/GENCODE_v34lift37.basic/
cp $sourceDir/gencode.v34lift37.basic.{{tss,promoter,intron,exon}.bed.gz,biotypes.txt} \
  $inDir

# get protein-coding ENSGs
cat $inDir/gencode.v34lift37.basic.biotypes.txt |
grep protein_coding | cut -f3 | sort -k1,1 -u \
> $pcENSGs

# fix up GENCODE annotations files for compatibility
(for ann in tss promoter intron exon ; do
  echo $ann
  
  # covert 1-based coords to 0-based, get ENSTs without version extension
  zcat $inDir/gencode.v34lift37.basic.${ann}.bed.gz | 
  awk -F'\t' '{gsub(/\.[0-9].*$/,"",$6);print $1,$2-1,$3,$4,$5,$6}' OFS='\t' |
  
  # remove symbol/ENSG columns except for tss file
  if [ $ann = 'tss' ] ; then cat ; else cut -f1-3,6 ; fi | 
  
  # sort
  sort -k1,1 -k2,2n |
  
  # save
  gzip > $outDir/gencode.v34lift37.basic.${ann}.bed.gz
  
done)
