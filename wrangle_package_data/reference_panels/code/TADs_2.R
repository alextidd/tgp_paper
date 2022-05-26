setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/") 
devtools::load_all("/working/lab_jonathb/alexandT/tgp")
TADs_metadata <- read_tibble("output/metadata.tsv", header = T) %>%
  dplyr::filter(object == "TADs")

TADs <- list() ; for(ct in TADs_metadata$celltype){
  print(ct)
  file <- paste0("output/TADs/", ct, ".bed")
  TADs[[ct]] <- import_BED(file)
}

# save
saveRDS(TADs, file = paste0("output/TADs/TADs.rds"))
