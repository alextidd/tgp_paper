#!/bin/bash
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/ ; cd $WKDIR
data_dir=$WKDIR/data/TADs/ ; mkdir $data_dir
out_dir=$WKDIR/output/TADs/ ; mkdir $out_dir
module load bedtools/2.27.1
module load ucsctools/20160223

echo "Trench ====================================================================="
mkdir $data_dir/Trench/
# manually coppied from https://docs.google.com/spreadsheets/d/1LAdGdI5F2ZWYWD-GPLV3fNd3RxUXdokJ7TADrQngGNA/edit?usp=sharing
cat $data_dir/Trench/T47D.TADs.nonTADS.bed |
cut -f1-3 \
> $out_dir/BRST.T47D.CNCR.bed

echo "Taberlay2016 ====================================================================="
mkdir $data_dir/Taberlay2016/
# manually copied from https://genome.cshlp.org/content/suppl/2016/05/06/gr.201517.115.DC1/Supplemental_Table_S4.txt
cp $data_dir/Taberlay2016/LNCAP.bed $out_dir/PR.LNCAP.CNCR.bed

echo "Rao2014 ====================================================================="
mkdir $data_dir/Rao2014/
# # hpcapp01 head node:
# wget -r -nH -nc -np -e robots=off -A domainlist.txt.gz --cut-dirs=100 \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/ \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/data/TADs/Rao2014/

# domainlist contains non-exclusive intervals, so must merge first (and exclude mouse samples)
for file in $(ls $data_dir/Rao2014/*domainlist.txt.gz | grep -v "CH12-LX") ; do
  filename=$(basename $file)
  cellline=${filename#*_} ; cellline=${cellline%%_*} ; cellline=${cellline^^}
  celltype=$(cat output/metadata.tsv | grep TADs | grep $cellline | cut -f1)
  echo $celltype
  zcat $file |
  sed '1d' |
  awk -F'\t' '{print "chr"$1"\t"$2"\t"$3}' |
  sort -k1,1 -k2,2n |
  bedtools merge  \
  > $out_dir/$celltype.bed
done

echo "Barutcu2015 ====================================================================="
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
> $out_dir/BRST.MCF7.CNCR.bed

echo "Du2021 ====================================================================="
mkdir $data_dir/Du2021/
# # hpcapp01 head node:
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158007/suppl/GSE158007_HCT116_40kb_TADs.bed.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/data/TADs/Smith2021/
zcat $data_dir/Du2021/GSE158007_HCT116_40kb_TADs.bed.gz |
bedtools merge \
> $out_dir/CLN.HCT116.CNCR.bed

echo "Javierre2016 ====================================================================="
mkdir $data_dir/Javierre2016/
# https://osf.io/u8tzp/ > $data_dir/Javierre2016/TAD_definitions.tar.gz
# tar xvzf $data_dir/Javierre2016/TAD_definitions.tar.gz -C $data_dir/Javierre2016/
cat $data_dir/Javierre2016/TADs_nCD4_mean_merged.bed |
awk -F'\t' -vOFS='\'t 'NR>1 {print "chr"$1,$2,$3}' \
> $out_dir/BLD.CD4T.bed

echo "Schmitt2016 ====================================================================="
mkdir $data_dir/Schmitt2016/
# manually copied from https://ars.els-cdn.com/content/image/1-s2.0-S2211124716314814-mmc4.xlsx, Sheet "OV"
cat $data_dir/Schmitt2016/1-s2.0-S2211124716314814-mmc4_OV.tsv \
> $out_dir/OVRY.bed

echo "LaGreca2022 ====================================================================="
mkdir $data_dir/LaGreca2022/
# manually copied from https://cdn.elifesciences.org/articles/66034/elife-66034-fig6-data1-v2.txt
liftOver \
  $data_dir/LaGreca2022/hg38.elife-66034-fig6-data1-v2.txt \
  data/hg/hg38ToHg19.over.chain \
  $data_dir/LaGreca2022/hg19.elife-66034-fig6-data1-v2.txt \
  $data_dir/LaGreca2022/unmapped.elife-66034-fig6-data1-v2.txt 
cat $data_dir/LaGreca2022/hg19.elife-66034-fig6-data1-v2.txt \
> $out_dir/ENDM.ISHIKAWA.CNCR.bed

echo "ENCODE ====================================================================="
mkdir $data_dir/ENCODE/
# # hpcapp01 head node #
# ( cd /working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/
# while read ENCFF ; do
#   echo $ENCFF
#   wget https://www.encodeproject.org/files/$ENCFF/@@download/$ENCFF.bed.gz \
#   -P data/TADs/ENCODE/
# done < <(cat output/metadata.tsv | awk '$3=="TADs" && $5 ~ "ENCFF" {print $5}')
# )
for file in $(ls $data_dir/ENCODE/*bed.gz) ; do
  accession=$(basename ${file%.bed.gz}) 
  celltype=$(cat output/metadata.tsv | grep $accession | cut -f1)
  zcat $file | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > $out_dir/$celltype.bed
done

echo "Generate RDS ====================================================================="
## Save to RDS ##
module load R/4.0.2
Rscript code/TADs_2.R


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