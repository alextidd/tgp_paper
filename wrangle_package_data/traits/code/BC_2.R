setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/")
library(tidyverse)

st14 <- read.delim("data/Michailidou2017/Supplementary_Table_14.tsv") %>%
  as_tibble()

index_SNPs <- st14 %>% 
  filter(var_name == top_snp) %>%
  distinct(top_snp, rs_number)

missing_rsIDs <- tibble(top_snp = paste0(setdiff(st14$top_snp, index_SNPs$top_snp))) %>%
  separate(top_snp, c("chrom", "position", "ref", "alt")) %>%
  mutate(chrom = paste0("chr", chrom))

write_delim(missing_rsIDs, "output/BC_Michailidou2017_LD/missing_index_rsIDs.tsv",
            delim = '\t')

# upload to UCSC Table Browser... > output/BC_Michailidou2017_LD/found_index_rsIDs.tsv