#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/
data_dir=${WKDIR}/data/TADs/ ; mkdir $data_dir
out_dir=${WKDIR}/output/TADs/ ; mkdir $out_dir


## Trench ## (breast)
# https://docs.google.com/spreadsheets/d/1LAdGdI5F2ZWYWD-GPLV3fNd3RxUXdokJ7TADrQngGNA/edit?usp=sharing
mkdir $data_dir/Trench/
cat $data_dir/Trench/T47D.TADs.nonTADS.bed |
cut -f1-3 \
> $out_dir/T47D.bed


## Taberlay 2016 ## (prostate)
mkdir $data_dir/Taberlay2016/
touch $data_dir/Taberlay2016/LNCAP.bed
# manually copied from https://genome.cshlp.org/content/suppl/2016/05/06/gr.201517.115.DC1/Supplemental_Table_S4.txt
cp $data_dir/Taberlay2016/LNCAP.bed $out_dir/


## Rao 2014 ## (blood, cervix, breast, vascular, lung, skin) 
mkdir $data_dir/Rao2014/
# # hpcapp01 head node:
# wget -r -nH -nc -np -e robots=off -A domainlist.txt.gz --cut-dirs=100 \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/ \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/data/TADs/Rao2014/

# domainlist contains non-exclusive intervals, so must merge first (and exclude mouse samples)
module load bedtools/2.27.1
for file in $(ls $data_dir/Rao2014/*domainlist.txt.gz | grep -v "CH12-LX") ; do
  filename=$(basename $file)
  celltype=${filename#*_} ; celltype=${celltype%%_*} ; celltype=${celltype^^}
  outfile=$out_dir/${celltype}.bed
  echo $celltype
  echo $outfile
  zcat ${file} |
  sed '1d' |
  awk -F'\t' '{print "chr"$1"\t"$2"\t"$3}' |
  sort -k1,1 -k2,2n |
  bedtools merge  \
  > $outfile
done


## Barutcu 2015 ##
mkdir $data_dir/Barutcu2015/
# # hpcapp01 head node:
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66733/suppl/GSE66733_Hi-C_MCF7_MCF10A_processed_HiCfiles.tar.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/data/TADs/Barutcu2015/

tar -xvzf \
  $data_dir/Barutcu2015/GSE66733_Hi-C_MCF7_MCF10A_processed_HiCfiles.tar.gz \
  -C $data_dir/Barutcu2015/
cat $data_dir/Barutcu2015/Hi-C_MCF7_MCF10A_processed_HiCfiles/TAD_boundaries/HiCStein-MCF7-WT__hg19__* | 
grep -v "header" |
sed -e 's/\(:\||\)/\t/g' |
cut -f3,5,6 |
# some non-exclusive TAD bounaries? -> merge
bedtools merge \
> $out_dir/MCF7.bed


## Smith 2021 ##
mkdir $data_dir/Smith2021/
# # hpcapp01 head node:
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158007/suppl/GSE158007_HCT116_40kb_TADs.bed.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/data/TADs/Smith2021/

zcat $data_dir/Smith2021/GSE158007_HCT116_40kb_TADs.bed.gz \
> $out_dir/HCT116.bed


## Javierre 2016 ##
mkdir $data_dir/Javierre2016/
# https://osf.io/u8tzp/ > $data_dir/Javierre2016/TAD_definitions.tar.gz
tar xvzf $data_dir/Javierre2016/TAD_definitions.tar.gz -C $data_dir/Javierre2016/
cat $data_dir/Javierre2016/TADs_nCD4_mean_merged.bed |
awk -F'\t' -vOFS='\'t 'NR>1{print "chr"$1,$2,$3}' \
> $out_dir/CD4.bed


# ## Meir 2020 ## ! cannot find TADs file
# mkdir $data_dir/Meir2020/
# # # hpcapp01 head node:
# # wget \
# # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144357/suppl/GSE144357_RAW.tar \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/data/TADs/Meir2020/
# 
# tar -xvf \
#   $data_dir/Meir2020/GSE144357_RAW.tar \
#   -C $data_dir/Meir2020/