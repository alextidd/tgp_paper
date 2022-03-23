library(tidyverse)
trait <- commandArgs(trailingOnly = T)[1]

wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/" ; setwd(wkdir)

# define in and out
assoc_file <- paste("data/traits/", trait, "variants/GWASCatalog_associations/", sep = "/") %>%
  list.files(pattern = "efotraits", full.names = T)
if(length(assoc_file) != 1){stop("Unique GWASCatalog association file not found in ", paste("data/traits/", trait, "variants/GWASCatalog_associations/", sep = "/"))}
index_file <- paste(wkdir, "output/traits/", trait, "variants/GWASCatalog_associations.proxies/SNiPA/index_SNPs.txt", sep = "/")

# get assoc SNPs
assoc_SNPs <- read.csv(
  file = assoc_file,
  quote = "\"") %>% as_tibble %>%
  filter(Location != "Mapping not available") %>%
  separate_rows(Variant.and.risk.allele, sep = ", ") %>%
  mutate(variant = Variant.and.risk.allele %>% gsub("\\-.*", "", .),
         p =  as.numeric(gsub(" x 10", "e", P.value))) %>%
  distinct(variant, p)

# get index SNPs
index_SNPs <- assoc_SNPs %>%
  filter(p < 5e-8) %>%
  group_by(variant) %>%
  filter(p == min(p)) %>%
  distinct

write.table(index_SNPs %>% dplyr::select(variant),
            index_file,
            col.names = F, row.names = F, quote = F, sep = "\t")
####################################################################################
