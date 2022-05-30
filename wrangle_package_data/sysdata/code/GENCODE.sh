#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/sysdata/ ; cd $WKDIR
data_dir=data/GENCODE/ ; mkdir $data_dir
out_dir=output/GENCODE/ ; mkdir $out_dir
filterENSGs=$out_dir/filter.gencode.v34lift37.basic.ENSGs.txt
pcENSGs=$out_dir/proteincoding.gencode.v34lift37.basic.ENSGs.txt

# copy files
sourceDir=/working/lab_georgiat/jonathB/DATA/gene_annotations/GENCODE/GENCODE_v34lift37.basic/
cp $sourceDir/gencode.v34lift37.basic.{{tss,promoter,intron,exon}.bed.gz,biotypes.txt} \
  $data_dir

# filter biotypes
biotypes=(antisense lincRNA lncRNA miRNA nonsense_mediated_decay protein_coding pseudogene)
awk -v biotypes="${biotypes[*]}" -F'\t' '
  BEGIN { n=split(biotypes,a," "); for (i=1;i<=n;++i) print i,biotype[a[i]] }
  $4 in biotype' $data_dir/gencode.v34lift37.basic.biotypes.txt |
cut -f3 | sort -k1,1 -u \
> $filterENSGs

# get protein-coding ENSGs
cat $data_dir/gencode.v34lift37.basic.biotypes.txt |
grep protein_coding | cut -f3 | sort -k1,1 -u \
> $pcENSGs

# fix up GENCODE annotations files for compatibility
(for ann in tss promoter intron exon ; do
  echo $ann
  
  zcat $data_dir/gencode.v34lift37.basic.${ann}.bed.gz | 
  
  # filtered biotypes only
  grep -F -f $filterENSGs |
  
  # remove older versions of duplicated ENST
  grep -v ENST00000442069.1 |
  grep -v ENST00000429980.1 |
  
  # covert 1-based coords to 0-based, get ENSTs without version extension
  awk -F'\t' '{gsub(/\.[0-9].*$/,"",$6);print $1,$2-1,$3,$4,$5,$6}' OFS='\t' |
  
  # remove symbol/ENSG columns except for tss file
  if [ $ann = 'tss' ] ; then cat ; else cut -f1-3,6 ; fi | 
  
  # sort
  sort -k1,1 -k2,2n |
  
  # save
  gzip > $out_dir/filter.gencode.v34lift37.basic.${ann}.bed.gz
  
done)
