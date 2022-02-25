#!/bin/bash
# Reformat all contact data to BEDPE with normalised scores
# chrA | startA | endA | chrB | startB | endB | normalised_score
# File names: study_sample_assay.bedpe

# reformat HiChIP_interactions.txt.gz: FitHiChIP interaction files (with header)
# -> (chr1 s1 e1 chr2 s2 e2 cc Coverage1 isPeak1 Bias1 Mapp1 GCContent1 RESites1 Coverage2 isPeak2 Bias2 Mapp2 GCContent2
#     RESites2 p exp_cc_Bias p_Bias dbinom_Bias P-Value_Bias Q-Value_Bias)
# -> Score = Q-Value_Bias

module load R/4.0.2
WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrange_package_data/ ; cd $WKDIR
dataDir=${WKDIR}/data/contact/ ; mkdir $dataDir
outputDir=${WKDIR}/output/contact/ ; mkdir -p $outputDir/pre_QC

echo "Trench ====================================================================="

# reformat FitHiChIP
Rscript code/wrangle_package_data/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP /working/lab_georgiat/jonathB/PROJECTS/trench_lab/DoD_screen_GenomicAssays/integrate_CCVs_HiChIP/data/collect_data/HiChIP/FitHiChIP/Tumcells/MCF7.FitHiChIP.Peak2ALL.Q0.01.bed \
  --out.bedpe $outputDir/pre_QC/Trench_MCF7_HiChIP.bedpe

echo "Shi 2021 ====================================================================="
mkdir $dataDir/Shi2021

# # wget on hpcapp01
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4568nnn/GSM4568375/suppl/GSM4568375_HaCaT_unst_27ac_helen_combined_stringent_merged_Q0.01_WashU.bed.gz \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4568nnn/GSM4568376/suppl/GSM4568376_MyLa_27ac_helen_combined_stringent_merged_Q0.01_WashU.bed.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Shi2021/

# reformat WashU longrange
for file in $dataDir/Shi2021/GSM*WashU.bed.gz ; do
  celltype=$(basename $file) ; celltype=${celltype#*_} ; celltype=${celltype%%_*} ; celltype=${celltype^^} ; echo $celltype
  zcat $file |
  cut -f1-4 |
  sed 's/:\|-\|,/\t/g' \
  > $outputDir/pre_QC/Shi2021_${celltype}_HiChIP.bedpe
done

echo "Chen 2021 ====================================================================="
mkdir $dataDir/Chen2021

# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE173699&format=file

# reformat - duplicates
for file in ${dataDir}/Chen2021/GSM*.gz ; do
  filename=$(basename ${file%.gz}) ; echo $filename
  zcat $file |
  sed 's/\s/\t/g' |
  cut -f1-6,8 \
  > ${outputDir}/pre_QC/Chen2021_${filename}
done

# merge duplicates and sum reads
mergefile=${outputDir}/pre_QC/Chen2021_HCT116_HiChIP.bedpe
join \
    -j 1 -t $'\t' -e 0 -a 1 -a 2 \
    -o 0,1.8,2.8 \
    <( cat ${outputDir}/pre_QC/Chen2021_*1.intra.loop_counts.bedpe |
       awk '{print $1"-"$2"-"$3"-"$4"-"$5"-"$6"\t"$0}' |
       sort -k1,1 ) \
    <( cat ${outputDir}/pre_QC/Chen2021_*2.intra.loop_counts.bedpe |
       awk '{print $1"-"$2"-"$3"-"$4"-"$5"-"$6"\t"$0}' |
       sort -k1,1 ) |
sed 's/-/\t/g' |
awk '{ for(i=7;i<=NF;i++) t+=$i ; print $0"\t"t ; t=0}' |
cut -f1-6,9 \
> ${mergefile}
rm -f ${outputDir}/pre_QC/Chen2021_GSM*

echo "Giambartolomei 2021 ====================================================================="
mkdir $dataDir/Giambartolomei2021

# # wget on hpcapp01
# wget https://wangftp.wustl.edu/~dli/Claudia/HiChIP_LNCaP.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Giambartolomei2021/

# (expand anchors - bins are 5000bp)
zcat $dataDir/Giambartolomei2021/HiChIP_LNCaP.gz |
sed 's/,\|:\|-/\t/g' |
awk '{FS=OFS="\t"}{print $1,$2-2499,$3+2499,$4,$5-2499,$6+2499,$7}' \
> $outputDir/pre_QC/Giambartolomei2021_LNCAP_HiChIP.bedpe

echo "Liu 2021 ====================================================================="
mkdir $dataDir/Liu2021

# # wget on hpcapp01
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5066nnn/GSM5066590/suppl/GSM5066590_LK2_HiChIP_H3K27ac.interactions.all.mango.txt.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Liu2021/

zcat $dataDir/Liu2021/GSM5066590_LK2_HiChIP_H3K27ac.interactions.all.mango.txt.gz | 
cut -f1-6,8 \
> $outputDir/pre_QC/Liu2021_LK2_HiChIP.bedpe

echo "Ma 2021 ====================================================================="
mkdir $dataDir/Ma2021

# # wget on hpcapp01
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178598/suppl/GSE178598_FitHiChIP.interactions_Q0.01_AoSMC_merge.bed.gz \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178598/suppl/GSE178598_FitHiChIP.interactions_Q0.01_HAEC_merge.bed.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Ma2021/

# reformat FitHiChIP
for file in $dataDir/Ma2021/*bed.gz ; do
  celltype=${file%_merge.bed.gz} ; celltype=${celltype##*_} ; CELLTYPE=${celltype^^}
  outfile=${outputDir}/pre_QC/Ma2021_${CELLTYPE}_HiChIP.bedpe

  Rscript code/wrangle_package_data/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP $file \
  --out.bedpe $outfile
done

echo "Corces 2020 ====================================================================="
mkdir $dataDir/Corces2020

# # wget on hpcapp01
# wget -r -nH -nc --no-parent -R "index.html*" -e robots=off --cut-dirs=100 \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4441nnn/GSM44418{30..41}/ \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Corces2020/

# reformat FitHiChIP
mkdir $outputDir/replicates
(for file in $dataDir/Corces2020/*bed.gz ; do
  celltype=$(basename ${file##*CTRL-}) ; celltype=${celltype##*RCLN-} ; celltype=${celltype/-*/}
  rep=$(ls $dataDir/Corces2020/*.bed.gz | grep $celltype | nl -nln | grep $(basename $file) | cut -f1 |tr -d ' ')
  (for chr in chr{1..22} ; do
    echo $celltype $rep $chr
    infile=${outputDir}/replicates/Corces2020_${celltype}_${chr}_FitHiChIP_replicate${rep}.bed
    outfile=${outputDir}/replicates/Corces2020_${celltype}_${chr}_HiChIP_replicate${rep}.bedpe
    zcat $file | awk -v chr=$chr 'NR==1 || $1==chr' > $infile    
    Rscript code/wrangle_package_data/FitHiChIP_to_bedpe.R \
    --in.FitHiChIP $infile \
    --out.bedpe $outfile
    rm -f $infile
  done)
done)


echo "Quality control ====================================================================="
Rscript code/wrangle_package_data/contact_2.R


# rejected files: #

# echo "Schmidt 2020 =====================================================================" # YY1, not H3K27ac
# mkdir $dataDir/Schmidt2020/
# 
# # # wget on hpcapp01
# # wget \
# # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_{GM12878_primary,HUVEC,HeLa}_HiCCUPS_looplist.txt.gz \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Schmidt2020/
# 
# for file in $dataDir/Rao2014/GSE63525_*_HiCCUPS_looplist.txt.gz ; do
#   celltype=$(basename $file) ; celltype=${celltype#*_} ; celltype=${celltype%%_*} ; celltype=${celltype^^} ; echo $celltype
#   zcat $file | 
#   sed 1d |
#   cut -f1-6 |
#   awk '{print "chr"$1,$2,$3,"chr"$4,$5,$6,"1"}' OFS='\t' \
#   > $outputDir/pre_QC/Schmidt2020_${celltype}_HiChIP.bedpe
# done

# echo "Zirkel 2018 =====================================================================" # CTCF, not H3K27ac
# mkdir $dataDir/Zirkel2018
# 
# # # wget on hpcapp01
# # wget \
# # https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2936nnn/GSM2936365/suppl/GSM2936365_prolif_huv.bedpe.gz \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Zirkel2018/
# 
# # add score 
# zcat $dataDir/Zirkel2018/GSM2936365_prolif_huv.bedpe.gz |
# awk '{print $0"\t1"}' \
# > $outputDir/pre_QC/Zirkel2018_HUVEC_HiChIP.bedpe

# echo "Shi 2021 =====================================================================" # old files
# mkdir $dataDir/Shi2021
# 
# # # wget on hpcapp01
# # accessions=( GSM4093604 GSM4093605 GSM4093606 )
# # samples=( HaCaT_stimulated_HiChIP HaCaT_unstimulated_HiChiP MyLa_HiChIP )
# # samples_out=( HaCaT_stimulated_HiChIP HaCaT_unstimulated_HiChIP MyLa_HiChIP ) # fix HiChiP typo
# # for (( i=0; i<${#accessions[@]}; i++ )) ; do
# # wget \
# # -O /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Shi2021/${samples_out[i]}_interactions.txt.gz \
# # https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4093nnn/${accessions[i]}/suppl/${accessions[i]}_${samples[i]}_interactions.txt.gz
# # done
# 
# # reformat FitHiChIP
# for file in ${dataDir}/Shi2021/*HiChIP_interactions.txt.gz ; do
#   celltype=$( basename $file | sed 's/_HiChIP.*//g' | sed 's/_/\./g' ) ; CELLTYPE=${celltype^^} ; echo $CELLTYPE
#   outfile=${outputDir}/pre_QC/Shi2021_${CELLTYPE}_HiChIP.bedpe
# 
#   Rscript code/wrangle_package_data/FitHiChIP_to_bedpe.R \
#   --in.FitHiChIP $file \
#   --out.bedpe $outfile
# done

# echo "Huang 2020 =====================================================================" # CTCF, not H3K27ac
# mkdir $dataDir/Huang2020
# 
# # # wget on hpcapp01
# # wget  \
# # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137849/suppl/GSE137849_CTCF_loops_mango_fdr0.05.tsv.gz \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/data/contact/Huang2020/
# 
# # get control reads
# zcat $dataDir/Huang2020/GSE137849_CTCF_loops_mango_fdr0.05.tsv.gz |
# sed '1d' |
# awk '$7>0' |
# cut -f1-7 \
# > $outputDir/pre_QC/Huang2020_HELAS3_HiChIP.bedpe