setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/") 

TADs_metadata <- read.delim("output/metadata.tsv", header = T) %>%
  dplyr::filter(object == "TADs")
TADs <- list() ; for(ct in TADs_metadata$celltype){
  print(ct)
  file <- paste0("output/TADs/", ct, ".bed")
  TADs[[ct]] <- read.delim(file, header = F) %>%
    dplyr::rename_with(., ~ c("chrom", "start", "end"), 1:3) %>%
    dplyr::as_tibble()
}

# save
saveRDS(TADs, file = paste0("output/TADs.rds"))
