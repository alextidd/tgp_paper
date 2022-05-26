#!/bin/bash
# Reformat all HiChIP data to BEDPE with normalised scores
# chrA | startA | endA | chrB | startB | endB | normalised_score
# File names: study_sample_assay.bedpe

# reformat HiChIP_interactions.txt.gz: FitHiChIP interaction files (with header)
# -> (chr1 s1 e1 chr2 s2 e2 cc Coverage1 isPeak1 Bias1 Mapp1 GCContent1 RESites1 Coverage2 isPeak2 Bias2 Mapp2 GCContent2
#     RESites2 p exp_cc_Bias p_Bias dbinom_Bias P-Value_Bias Q-Value_Bias)
# -> Score = -log10(Q-Value_Bias)

module load R/4.0.2
get_cellline() { cellline=${celltype#*.} ; cellline=${cellline%.*} ; }

WKDIR=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/ ; cd $WKDIR
data_dir=$WKDIR/data/HiChIP/ ; mkdir $data_dir
out_dir=$WKDIR/output/HiChIP/ ; mkdir -p $out_dir/pre_QC

echo "Trench ====================================================================="

# FitHiChIP to bedpe
Trench_FitHiChIP_dir=/working/lab_georgiat/jonathB/PROJECTS/trench_lab/DoD_screen_GenomicAssays/integrate_CCVs_HiChIP/data/collect_data/HiChIP/FitHiChIP/
(for celltype in BRST.HMEC BRST.MCF7.CNCR BRST.T47D.CNCR ; do echo $celltype
  get_cellline
  if [[ "$cellline" == "HMEC" ]] ; then cellline=B80T5 ; fi
  file=$(ls $Trench_FitHiChIP_dir/*cells/$cellline.FitHiChIP.Peak2ALL.Q0.01.bed)
  Rscript code/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP $file \
  --out.bedpe $out_dir/pre_QC/$celltype.bedpe
done)

echo "Shi2021 ====================================================================="
mkdir $data_dir/Shi2021

# # wget on hpcapp01
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4568nnn/GSM4568375/suppl/GSM4568375_HaCaT_unst_27ac_helen_combined_stringent_merged_Q0.01_WashU.bed.gz \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4568nnn/GSM4568376/suppl/GSM4568376_MyLa_27ac_helen_combined_stringent_merged_Q0.01_WashU.bed.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Shi2021/

# reformat WashU longrange
for celltype in SKIN.HACAT BLD.MYLA.CNCR ; do echo $celltype
  get_cellline
  file=$(find $data_dir/Shi2021/ -iname GSM*$cellline*WashU.bed.gz)
  zcat $file |
  cut -f1-4 |
  sed 's/:\|-\|,/\t/g' \
  > $out_dir/pre_QC/$celltype.bedpe
done

echo "Chen2021 ====================================================================="
mkdir $data_dir/Chen2021

# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE173699&format=file

# reformat - duplicates
for file in $data_dir/Chen2021/GSM*.gz ; do
  filename=$(basename ${file%.gz}) ; echo $filename
  zcat $file |
  sed 's/\s/\t/g' |
  cut -f1-6,8 \
  > $out_dir/pre_QC/Chen2021_$filename
done

# merge duplicates and sum reads
join \
    -j 1 -t $'\t' -e 0 -a 1 -a 2 \
    -o 0,1.8,2.8 \
    <( cat ${out_dir}/pre_QC/Chen2021_*1.intra.loop_counts.bedpe |
       awk '{print $1"-"$2"-"$3"-"$4"-"$5"-"$6"\t"$0}' |
       sort -k1,1 ) \
    <( cat ${out_dir}/pre_QC/Chen2021_*2.intra.loop_counts.bedpe |
       awk '{print $1"-"$2"-"$3"-"$4"-"$5"-"$6"\t"$0}' |
       sort -k1,1 ) |
sed 's/-/\t/g' |
awk '{ for(i=7;i<=NF;i++) t+=$i ; print $0"\t"t ; t=0}' |
cut -f1-6,9 \
> $out_dir/pre_QC/CLN.HCT116.CNCR.bedpe
rm -f $out_dir/pre_QC/Chen2021_GSM*

echo "Giambartolomei2021 ====================================================================="
mkdir $data_dir/Giambartolomei2021

# # wget on hpcapp01
# wget https://wangftp.wustl.edu/~dli/Claudia/HiChIP_LNCaP.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Giambartolomei2021/

# (expand anchors - bins are 5000bp)
zcat $data_dir/Giambartolomei2021/HiChIP_LNCaP.gz |
sed 's/,\|:\|-/\t/g' |
awk '{FS=OFS="\t"}{print $1,$2-2499,$3+2499,$4,$5-2499,$6+2499,$7}' \
> $out_dir/pre_QC/PR.LNCAP.CNCR.bedpe

echo "Liu2021 ====================================================================="
mkdir $data_dir/Liu2021

# # wget on hpcapp01
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5066nnn/GSM5066590/suppl/GSM5066590_LK2_HiChIP_H3K27ac.interactions.all.mango.txt.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Liu2021/

# mango cols: (c1 s1 e1 c2 s2 e2 PETs FDR)
# filter to significant interactions (FDR < 0.05) and -log10(FDR) to get loop score
zcat $data_dir/Liu2021/GSM5066590_LK2_HiChIP_H3K27ac.interactions.all.mango.txt.gz |
awk -F'\t' -vOFS='\t' '$8 + 0 < 0.05 {neglog10p= -log($8+0)/log(10) ; print $1,$2,$3,$4,$5,$6,neglog10p}' \
> $out_dir/pre_QC/LNG.LK2.CNCR.bedpe

echo "Ma2021 ====================================================================="
mkdir $data_dir/Ma2021

# # wget on hpcapp01
# wget \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178598/suppl/GSE178598_FitHiChIP.interactions_Q0.01_AoSMC_merge.bed.gz \
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178598/suppl/GSE178598_FitHiChIP.interactions_Q0.01_HAEC_merge.bed.gz \
# -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Ma2021/

# reformat FitHiChIP
for celltype in VAS.AOSMC VAS.HAEC ; do echo $celltype
  get_cellline
  file=$(find $data_dir/Ma2021/ -iname *$cellline*.bed.gz)
  Rscript code/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP $file \
  --out.bedpe $out_dir/pre_QC/$celltype.bedpe
done

echo "Bhattacharyya2019 ====================================================================="
mkdir $data_dir/Bhattacharyya2019/
# WashU epigenome browser sessions: (http://epigenomegateway.wustl.edu/browser/)
# GM12878 H3K27ac loop callsSession ID: b491c3d0-65f7-11e9-b334-5ff263937318
# K562 H3K27ac loop callsSession ID: 019e06b0-65f4-11e9-921c-577d3df57445

# https://zenodo.org/record/3255048/files/FitHiChIP_Source_Data_June2019.zip?download=1
# (cd $data_dir/Bhattacharyya2019/ ; unzip FitHiChIP_Source_Data_June2019.zip)

(for celltype in BLD.GM12878 BLD.K562.CNCR BLD.CD4.TCELL ; do echo $celltype
  get_cellline

  infile=$(ls $data_dir/Bhattacharyya2019/FitHiChIP_Source_Data_June2019/$cellline*/H3K27*ac/Combined_Replicates_HIChIP_Loops/Table_*.xlsx)
  outfile=$out_dir/pre_QC/$celltype
  if [[ celltype == 'BLD.CD4.TCELL' ]] ; then outfile=$out_dir/pre_QC/BLD.CD4T ; fi

  # xlsx to tsv
  libreoffice --headless --convert-to csv --outdir $out_dir/pre_QC/ $infile

  # tsv to FitHiChIP
  ( echo -e "chr1\ts1\te1\tchr2\ts2\te2\tQ-Value_Bias" ;
    tac $out_dir/pre_QC/Table_*.csv |
    awk '/Chromosome/ {exit} 1' |
    tac |
    awk -F, -vOFS='\t' '{print $1,$2,$3,$1,$4,$5,$14}'
  ) | cat > $outfile.FitHiChIP
  rm -f $out_dir/pre_QC/Table_*.csv

  # FitHiChIP to BEDPE
  Rscript code/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP $outfile.FitHiChIP \
  --out.bedpe $outfile.bedpe
  rm -f $outfile.FitHiChIP
done)

echo "Quinn2021 ====================================================================="
mkdir $data_dir/Quinn2021/

cp \
  /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/HiChIP/TOV112D/FitHiChIP/output/FitHiChIP_Peak2ALL_b5000_L5000_U2000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q0.01.bed \
  $data_dir/Quinn2021/

# reformat FitHiChIP
Rscript code/FitHiChIP_to_bedpe.R \
--in.FitHiChIP $data_dir/Quinn2021/FitHiChIP.interactions_FitHiC_Q0.01.bed \
--out.bedpe $out_dir/pre_QC/OVRY.TOV112D.CNCR.bedpe

echo "Omara2019 ====================================================================="
mkdir $data_dir/Omara2019/

cp \
  /working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/epigenomic_data/epi_data_formatting/output/HiChIP/Ishikawa/FitHiChIP/FitHiChIP_Peak2ALL_b5000_L5000_U2000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q0.01.bed \
  $data_dir/Omara2019/

# reformat FitHiChIP
Rscript code/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP $data_dir/Omara2019/FitHiChIP.interactions_FitHiC_Q0.01.bed \
  --out.bedpe $out_dir/pre_QC/ENDM.ISHIKAWA.CNCR.bedpe

echo "Chandra2021 ====================================================================="
mkdir $data_dir/Chandra2021/

# > $data_dir/Chandra2021/41588_2020_745_MOESM3_ESM_ST3.tsv
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-00745-3/MediaObjects/41588_2020_745_MOESM3_ESM.xlsx
# manually copied from Supplementary Table 3 - Summary of promoter interactions (H3K27ac HiChIP) (F4:K,N)

( echo -e "chr1\ts1\te1\tchr2\ts2\te2\tQ.Value_Bias" ;
  sed 1d $data_dir/Chandra2021/41588_2020_745_MOESM3_ESM_ST3.tsv ; ) |
cat > $data_dir/Chandra2021/41588_2020_745_MOESM3_ESM_ST3.FitHiChIP

# reformat FitHiChIP
Rscript code/FitHiChIP_to_bedpe.R \
  --in.FitHiChIP $data_dir/Chandra2021/41588_2020_745_MOESM3_ESM_ST3.FitHiChIP \
  --out.bedpe $out_dir/pre_QC/BLD.MONO.bedpe

echo "Quality control ====================================================================="
Rscript code/HiChIP_2.R

###################
# rejected files: #

# echo "Corces 2020 =====================================================================" # not using brain panel anymore
# mkdir $data_dir/Corces2020
# 
# # # wget on hpcapp01
# # wget -r -nH -nc --no-parent -R "index.html*" -e robots=off --cut-dirs=100 \
# # https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4441nnn/GSM44418{30..41}/ \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Corces2020/
# 
# # reformat FitHiChIP
# mkdir $out_dir/replicates
# (for file in $data_dir/Corces2020/*bed.gz ; do
#   celltype=$(basename ${file##*CTRL-}) ; celltype=${celltype##*RCLN-} ; celltype=${celltype/-*/}
#   rep=$(ls $data_dir/Corces2020/*.bed.gz | grep $celltype | nl -nln | grep $(basename $file) | cut -f1 |tr -d ' ')
#   (for chr in chr{1..22} ; do
#     echo $celltype $rep $chr
#     infile=${out_dir}/replicates/Corces2020_${celltype}_${chr}_FitHiChIP_replicate${rep}.bed
#     outfile=${out_dir}/replicates/Corces2020_${celltype}_${chr}_HiChIP_replicate${rep}.bedpe
#     zcat $file | awk -v chr=$chr 'NR==1 || $1==chr' > $infile    
#     Rscript code/FitHiChIP_to_bedpe.R \
#     --in.FitHiChIP $infile \
#     --out.bedpe $outfile
#     rm -f $infile
#   done)
# done)

# echo "Schmidt 2020 =====================================================================" # YY1, not H3K27ac
# mkdir $data_dir/Schmidt2020/
# 
# # # wget on hpcapp01
# # wget \
# # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_{GM12878_primary,HUVEC,HeLa}_HiCCUPS_looplist.txt.gz \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Schmidt2020/
# 
# for file in $data_dir/Rao2014/GSE63525_*_HiCCUPS_looplist.txt.gz ; do
#   celltype=$(basename $file) ; celltype=${celltype#*_} ; celltype=${celltype%%_*} ; celltype=${celltype^^} ; echo $celltype
#   zcat $file | 
#   sed 1d |
#   cut -f1-6 |
#   awk '{print "chr"$1,$2,$3,"chr"$4,$5,$6,"1"}' OFS='\t' \
#   > $out_dir/pre_QC/Schmidt2020_${celltype}_HiChIP.bedpe
# done

# echo "Zirkel 2018 =====================================================================" # CTCF, not H3K27ac
# mkdir $data_dir/Zirkel2018
# 
# # # wget on hpcapp01
# # wget \
# # https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2936nnn/GSM2936365/suppl/GSM2936365_prolif_huv.bedpe.gz \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Zirkel2018/
# 
# # add score 
# zcat $data_dir/Zirkel2018/GSM2936365_prolif_huv.bedpe.gz |
# awk '{print $0"\t1"}' \
# > $out_dir/pre_QC/Zirkel2018_HUVEC_HiChIP.bedpe

# echo "Shi 2021 =====================================================================" # old files
# mkdir $data_dir/Shi2021
# 
# # # wget on hpcapp01
# # accessions=( GSM4093604 GSM4093605 GSM4093606 )
# # samples=( HaCaT_stimulated_HiChIP HaCaT_unstimulated_HiChiP MyLa_HiChIP )
# # samples_out=( HaCaT_stimulated_HiChIP HaCaT_unstimulated_HiChIP MyLa_HiChIP ) # fix HiChiP typo
# # for (( i=0; i<${#accessions[@]}; i++ )) ; do
# # wget \
# # -O /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Shi2021/${samples_out[i]}_interactions.txt.gz \
# # https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4093nnn/${accessions[i]}/suppl/${accessions[i]}_${samples[i]}_interactions.txt.gz
# # done
# 
# # reformat FitHiChIP
# for file in ${data_dir}/Shi2021/*HiChIP_interactions.txt.gz ; do
#   celltype=$( basename $file | sed 's/_HiChIP.*//g' | sed 's/_/\./g' ) ; CELLTYPE=${celltype^^} ; echo $CELLTYPE
#   outfile=${out_dir}/pre_QC/Shi2021_${CELLTYPE}_HiChIP.bedpe
# 
#   Rscript code/FitHiChIP_to_bedpe.R \
#   --in.FitHiChIP $file \
#   --out.bedpe $outfile
# done

# echo "Huang 2020 =====================================================================" # CTCF, not H3K27ac
# mkdir $data_dir/Huang2020
# 
# # # wget on hpcapp01
# # wget  \
# # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137849/suppl/GSE137849_CTCF_loops_mango_fdr0.05.tsv.gz \
# # -P /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/HiChIP/Huang2020/
# 
# # get control reads
# zcat $data_dir/Huang2020/GSE137849_CTCF_loops_mango_fdr0.05.tsv.gz |
# sed '1d' |
# awk '$7>0' |
# cut -f1-7 \
# > $out_dir/pre_QC/Huang2020_HELAS3_HiChIP.bedpe