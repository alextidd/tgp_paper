opts <- commandArgs(trailingOnly = T)
trait <- opts[1] # trait = "BC" #
celltypes <- opts[2] # celltypes = "enriched_tissues" #

baseDir<- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/compare_methods/" ; setwd(wkdir)
variant_to_gene_max_distance = 2e6
max_n_drivers_per_CS = Inf

devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

cat(trait, celltypes, "#############\n")
outDir <- paste("output", trait, celltypes, "", sep = "/")
      
# get trait drivers / variants ====
driversFile <- paste0(baseDir, "wrangle_package_data/traits/output/", trait, "/known_genes.txt")
drivers <- read_tibble(driversFile)$V1 %>%
  check_driver_symbols(driversFile)
variantsFile <- paste0(baseDir, "wrangle_package_data/traits/output/", trait, "/variants/bed")
variants <- import_BED(variantsFile,
                       metadata_cols = c("variant", "cs"))

# ##
# # ABC settings #
# tested_ABC <- read_tibble("../../ABC-GWAS-Paper/ABC-Max/out/ABC/IBD/ABC_settings/PR_tested_CredibleSets.tsv", header = T) 
# variants <- variants %>% dplyr::filter(cs %in% tested_ABC$CredibleSet)
# variant_to_gene_max_distance = 1e6
# max_n_drivers_per_CS = 1
# ##

# txv masterlist ====
txv_master <- variants %>%
  # expand variant coords to search range
  valr::bed_slop(both = variant_to_gene_max_distance,
                 genome = ChrSizes,
                 trim = T) %>%
  # get all protein-coding TSSs within search range of each variant
  bed_intersect_left(., TSSs %>% dplyr::filter(ensg %in% pcENSGs),
                     suffix = c(".variant", ".TSS")) %>%
  # restore variant coords
  dplyr::select(-c(start, end)) %>%
  dplyr::left_join(variants, by = c("chrom", "variant", "cs")) %>%
  # add distance
  dplyr::mutate(distance = abs(end.TSS - end))

# gxc masterlist ====
gxc_master <- txv_master %>% 
  # convert to cs-x-gene level 
  dplyr::distinct(cs, symbol)

# generate distance-based predictions ====
distance_predictions <- txv_master %>%
  dplyr::group_by(variant) %>%
  # Calculate inverse of the absolute bp distance for each variant-transcript pair
  dplyr::mutate(invDistanceToTSS = 1/distance,
                # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                invDistanceToTSSRank = 1/rank(distance, ties.method = "min")) %>%
  dplyr::select(variant, cs, symbol, dplyr::starts_with("inv"))

# get method predictions ====
method_predictions <- read_tibble(
  gzfile(paste("output", trait, celltypes, "predictions_long.tsv.gz", sep = "/")),
  header = T
  ) %>% 
  # get max score per gxc per method
  dplyr::group_by(cs, symbol, method) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::distinct() %>%
  # widen
  tidyr::pivot_wider(id_cols = c(cs, symbol),
                     names_from = method,
                     values_from = score) 

# gather all predictions ===
predictions <- gxc_master %>%
  dplyr::left_join(method_predictions) %>%
  dplyr::left_join(distance_predictions)
predictions[is.na(predictions)] <- 0
# # TODO: add measure of inter-method concurrence (common SNP-gene links)
# predictions %>%
#   tidyr::pivot_longer(where(is.numeric),
#                       values_to = "score",
#                       names_to = "method") %>%
#   dplyr::group_by(cs, symbol, method) %>%
#   dplyr::filter(score == max(score), score > 0) %>%
#   dplyr::group_by(cs, symbol) %>%
#   dplyr::summarise(predictions = dplyr::n_distinct(method)) %>%
#   dplyr::group_by(predictions) %>%
#   dplyr::count()

# performance ====
performance <- predictions %>%
  get_PR(.,
         txv_master,
         drivers,
         pcENSGs,
         max_n_drivers_per_CS)

# write performance
write_tibble(performance$summary %>% 
               dplyr::mutate(trait = trait) %>%
               dplyr::select(trait, everything()) %>%
               dplyr::select(-prediction_type), 
             paste0(outDir, "performance.tsv"))

# plot performance ===

# plot labels
performance_to_plot <- performance %>%
  purrr::map(
    ~ .x %>%
      dplyr::group_by(prediction_method) %>%
      dplyr::mutate(prediction_method = prediction_method %>%
                      paste0(" (", PR_AUC %>% max %>% round(2), ")"))
  )

# subtitle
plot_subtitle <- paste0(
  "\nTrait = ", trait,
  ", Celltypes = " , celltypes,
  "\nmax n drivers per CS = ", max_n_drivers_per_CS, 
  ", max distance = ", variant_to_gene_max_distance)

PR <- performance_to_plot %>%
  plot_PR +
  ggplot2::ggtitle(paste("Precision-Recall\n", plot_subtitle)) +
  ggrepel::geom_text_repel(
    data = . %>% dplyr::filter(prediction_type == "max"))
AUPRC <- performance_to_plot %>%
  plot_AUPRC +
  ggplot2::ggtitle(paste("Area under Precision-Recall curve\n", plot_subtitle))

# save plots
{
  pdf(paste0(outDir, "performance.pdf"), onefile = T, width = 15)
  gridExtra::grid.arrange(
    PR + ggplot2::theme(legend.position = "none"),
    AUPRC + ggplot2::theme(legend.position = "none"),
    ncol = 2)
  dev.off()
}


