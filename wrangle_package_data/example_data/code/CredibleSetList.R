# generate ABC-compatible cs lists for BC and PrCa
trait <- commandArgs(trailingOnly = T)[1] # trait = "PrCa"
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/" ; setwd(wkdir)
inDir <- paste0(wkdir, "output/Traits/", trait, "/")
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")
library(magrittr)

# read variants
variants <- read_tibble(
  paste0(inDir, "VariantList.tsv"), 
  header = T
) %>%
  dplyr::mutate(CredibleSet = as.character(CredibleSet)) 

# get best SNPs
if(trait %in% c("PrCa_Giambartolomei2021", "PrCa_Dadaev2018")){
  best <- variants %>%
    dplyr::group_by(CredibleSet) %>% 
    # best SNP
    dplyr::filter(PosteriorProb == max(PosteriorProb)) %>%
    # if more than one bestSNP per CS, just randomly sample one (this only affects the distance predictions, not the ABC predictions)
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::transmute(chrom,
                     BestSNPPos = end,
                     BestSNP = variant,
                     CredibleSet)
}
if(trait %in% c("BC")){
  best <- read_tibble(paste0(inDir, "best_SNP_per_CS.tsv"), header = T)
}

variants_annotated <- 
  # add Coding, SpliceSite, Promoter
  list(
    promoters %>% dplyr::transmute(chrom, start, end, Promoter = T),
    splicesite %>% dplyr::transmute(
      chrom,
      start = position - 1,
      end = position,
      SpliceSite = T
    ),
    dplyr::bind_rows(missense, nonsense) %>% dplyr::transmute(
      chrom,
      start = position - 1,
      end = position,
      Coding = T
    )
  ) %>% 
  purrr::map(~ bed_intersect_left(variants, ., keepBcoords = F)) %>%
  purrr::reduce(dplyr::full_join) %>%
  dplyr::full_join(variants) %>%
  dplyr::mutate(dplyr::across(c(Promoter, SpliceSite, Coding), ~ tidyr::replace_na(.x, FALSE)))

# generate set summaries (for ABC input, so use chr,start,end)
sets <- variants_annotated %>% 
  dplyr::group_by(CredibleSet) %>%
  # summarise
  dplyr::summarise(chr = unique(chrom),
                   start = min(start),
                   end = max(end),
                   nSNP = dplyr::n_distinct(variant),
                   AnyCoding = any(Coding),
                   AnyPromoter = any(Promoter),
                   AnySpliceSite = any(SpliceSite),
                   Disease = trait)  %>%
  dplyr::ungroup() %>%  
  dplyr::select(chr, start, end, CredibleSet, everything()) %>%
  # add BestSNP
  dplyr::left_join(best %>% dplyr::rename(chr = chrom)) 

# write
write.table(sets,
            paste0(inDir, "CredibleSetList.tsv"),
            row.names = F,
            quote = F,
            sep = "\t")
write.table(sets,
            paste0(inDir, "CredibleSetList.bed"),
            col.names = F, row.names = F,
            quote = F,
            sep = "\t")

