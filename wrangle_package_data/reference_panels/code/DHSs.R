# only mapped regions
DHSs <- read.delim("data/DHS.bed", header = F) %>% 
  dplyr::transmute(chrom = V1, start = V2, end = V3,
            DHS = paste0(chrom, ":", start, "-", end)) %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X")))
  
saveRDS(DHSs, "output/DHSs.rds")
