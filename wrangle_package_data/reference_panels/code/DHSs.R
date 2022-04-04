library(tidyverse)

# only mapped regions
DHSs <- read.delim("data/DHS.bed", header = F) %>% 
  transmute(chrom = V1, start = V2, end = V3,
            DHS = paste0(chrom, ":", start, "-", end)) %>%
  filter(chrom %in% paste0("chr", c(1:22, "X")))
  
dir.create("output/DHSs/")
saveRDS(DHSs, "output/DHSs/DHSs.rds")
