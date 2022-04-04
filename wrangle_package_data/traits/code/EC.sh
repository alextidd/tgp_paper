#!/bin/bash
. code/utils.sh

##### VARIANTS #####

# Wang2022 LD #
mkdir output/EC_Wang2022_LD/
# 1.  get index SNPs from CSs
cat data/Wang2022/Table2.tsv | sed 1d | cut -f2 | tr ',' '\n' |tr -d '^ ' \
> output/EC_Wang2022_LD/index_SNPs.txt
# 2.  get proxy SNPs (r2>0.8 w/ an index SNP) from 
#     https://snipa.helmholtz-muenchen.de/snipa3/
# 3.  generate variant list
SNiPA_to_variants output/EC_Wang2022_LD/

##### KNOWN GENES #####

# Dietlein2020 #
get_kgs EC \
<(  cat data/Dietlein2020/NIHMS1546846-supplement-2_Supplementary_Table_3.tsv |
    awk -F'\t' -vOFS='\t' '$2 == "Endometrium" && $4 != "level D" {print $1}' )


