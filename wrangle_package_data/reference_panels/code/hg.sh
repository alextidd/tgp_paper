#!/bin/bash 

# hg19.genome from https://drive.google.com/file/d/12BxeHUiJIo-jkuzLvp0pOnF3v4RFi2uu/view?usp=sharing

# ##### on hpcapp01 head node:
# dir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/data/hg
# mkdir -p $dir
# cd $dir
# wget \
#   'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' \
#   -O hg38ToHg19.over.chain.gz
# gunzip hg38ToHg19.over.chain.gz
# #####