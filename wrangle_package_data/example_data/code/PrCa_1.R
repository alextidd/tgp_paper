library(tidyverse)

wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/" ; setwd(wkdir)
outDir <- paste0(wkdir, "output/Traits/PrCa_GWASCatalog_index_proxies/")

assoc_SNPs <- read.csv(
  file = "data/Traits/PrCa/GWASCatalog/efotraits_EFO_0001663-associations-2022-03-4.csv",
  quote = "\"") %>% as_tibble %>%
  filter(Location != "Mapping not available") %>%
  separate_rows(Variant.and.risk.allele, sep = ", ") %>%
  mutate(variant = Variant.and.risk.allele %>% gsub("\\-.*", "", .),
         p =  as.numeric(gsub(" x 10", "e", P.value))) %>%
  distinct(variant, p)

index_SNPs <- assoc_SNPs %>%
  filter(p < 5e-8) %>%
  group_by(variant) %>%
  filter(p == min(p)) %>%
  distinct

write.table(index_SNPs %>% dplyr::select(variant),
            paste0(outDir, "index_SNPs.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")
####################################################################################
