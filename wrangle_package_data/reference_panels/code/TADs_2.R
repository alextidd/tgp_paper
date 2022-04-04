setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/") ; devtools::load_all("/working/lab_jonathb/alexandT/tgp")
metadata <- read_tibble("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/output/metadata.tsv", header = T)

TADs <- list() ; for(file in list.files("output/TADs/", pattern = "bed", full.names = T)){
  celltype <- file %>% 
    basename %>%
    {dplyr::filter(metadata, file == ., object == "TADs")} %>%
    dplyr::pull(celltype)
  print(celltype)
  TADs[[celltype]] <- import_BED(file)
}

# save
saveRDS(TADs, file = paste0("output/TADs/TADs.rds"))
