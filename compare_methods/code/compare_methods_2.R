opts <- commandArgs(trailingOnly = T)
trait <- opts[1] # trait = "BC" #
celltypes <- opts[2] # celltypes = "enriched_tissues" #

baseDir <- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/compare_methods/" ; setwd(wkdir)
variant_to_gene_max_distance = 1e6
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")
library(tidyverse)

cat(trait, celltypes, "#############\n")
outDir <- paste("output", trait, celltypes, "", sep = "/")
      
# get trait drivers / variants ====
drivers <- read_tibble(paste0(baseDir, "wrangle_package_data/output/Traits/", trait, "/", trait, ".Drivers.txt"))$V1 %>%
  # add info 
  { dplyr::filter(TSSs, symbol %in% .)}
variants <- read_tibble(paste0(baseDir, "wrangle_package_data/output/Traits/", trait, "/", trait, ".VariantList.bed")) %>%
  dplyr::select(chrom = V1, start = V2, end = V3, variant = V4, cs = V5)

# txv masterlist ====
txv_master <- variants %>%
  # expand variant coords to search range
  valr::bed_slop(both = variant_to_gene_max_distance,
                 genome = ChrSizes,
                 trim = T) %>%
  # get all TSSs within search range of each variant
  bed_intersect_left(., TSSs,
                     suffix = c(".variant", ".TSS")) %>%
  # restore variant coords
  dplyr::select(-c(start, end)) %>%
  dplyr::inner_join(variants) %>%
  # add distance
  dplyr::mutate(distance = abs(end.TSS - end))
  

# cs <- read_tibble("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data//data/Traits/IBD//IBDCombined.set1-2.cs.txt", header = T)
# cs_to_consider <- cs %>% dplyr::filter(!AnyCoding, !AnySpliceSite)
# varscores <- read_tibble("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data//data/Traits/IBD//IBDCombined.set1-2.variant.list.txt", header = T)
# vars_to_consider <- varscores %>% dplyr::filter(PosteriorProb > 0.1)
# variants <- read_tibble(paste0(baseDir, "wrangle_package_data/output/Traits/", trait, "/", trait, ".VariantList.bed")) %>%
#   dplyr::select(chrom = V1, start = V2, end = V3, variant = V4, cs = V5)
# txv_master <- variants %>%
#   # expand variant coords to search range
#   valr::bed_slop(both = variant_to_gene_max_distance,
#                  genome = ChrSizes,
#                  trim = T) %>%
#   # get all TSSs within search range of each variant
#   bed_intersect_left(., TSSs,
#                      suffix = c(".variant", ".TSS")) %>%
#   # restore variant coords
#   dplyr::select(-c(start, end)) %>%
#   dplyr::inner_join(variants) %>%
#   # add distance
#   dplyr::mutate(distance = abs(end.TSS - end)) %>% 
#   dplyr::filter(cs %in% cs_to_consider$CredibleSet, variant %in% vars_to_consider$variant) 
# txv_master %>%
#   dplyr::distinct(cs, symbol) %>%
#   dplyr::filter(symbol %in% drivers$symbol) %>%
#   dplyr::group_by(cs) %>%
#   dplyr::filter(dplyr::n() == 1)

# generate distance-based prediction methods ====
distance_predictions <- txv_master %>%
  dplyr::group_by(variant) %>%
  # Calculate inverse of the absolute bp distance for each variant-transcript pair
  dplyr::mutate(invDistanceToTSS = 1/distance,
                # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                invDistanceToTSSRank = 1/rank(distance, ties.method = "min")) %>%
  # pivot
  tidyr::pivot_longer(c(invDistanceToTSS, invDistanceToTSSRank),
                      names_to = "prediction_method",
                      values_to = "score") %>%
  dplyr::group_by(cs, symbol, prediction_method) %>%
  dplyr::summarise(score = max(score))

# gxc masterlist ====
# The gene-x-cs universe (masterlist of all possible gene x credible set pairs < variant_to_gene_max_distance apart)
# only for CSs within variant_to_gene_max_distance of a trait driver TSS
gxc_master <- txv_master %>% 
  # convert to cs-x-gene level 
  dplyr::distinct(cs, symbol)

# get method predictions ====
predictions_long <- read_tibble(
  gzfile(paste("output", trait, celltypes, "predictions_long.tsv.gz", sep = "/")),
  header = T
  ) %>% 
  # rename method -> prediction_method
  dplyr::rename(prediction_method = method) %>%
  # add distance predictions
  dplyr::bind_rows(distance_predictions)

predictions <- gxc_master %>%
  # get full gxc universe per method
  tidyr::crossing(dplyr::tibble(prediction_method = unique(predictions_long$prediction_method))) %>%
  dplyr::left_join(predictions_long) %>%
  dplyr::mutate(score = replace_na(score, 0)) %>%
  # max score per cs x symbol x prediction_method
  dplyr::group_by(cs, symbol, prediction_method) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::distinct() %>%
  # max score per cs x prediction_method
  dplyr::group_by(cs, prediction_method) %>%
  dplyr::mutate(max = as.numeric(score == max(score) & score > 0)) %>%
  dplyr::ungroup() %>%
  # gather prediction types
  tidyr::pivot_longer(c(score, max),
                      names_to = "prediction_type",
                      values_to = "prediction") %>%
  # add driver annotation
  dplyr::mutate(driver = symbol %in% drivers$symbol) %>%
  # get testable
  get_testable 

# performance ====
performance <- list()
performance$summary <- predictions %>%
  # score is not binary, cannot be summarised, filter out
  dplyr::filter(prediction_type == "max") %>%
  dplyr::mutate(prediction = as.logical(prediction)) %>%
  dplyr::group_by(prediction_method, prediction_type) %>%
  dplyr::group_modify(
    ~ data.frame(True = .x %>% condition_n_gene_x_cs_pairs(driver),
                 Positive = .x %>% condition_n_gene_x_cs_pairs(prediction),
                 TP = .x %>% condition_n_gene_x_cs_pairs(prediction & driver),
                 FP = .x %>% condition_n_gene_x_cs_pairs(prediction & !driver),
                 TN = .x %>% condition_n_gene_x_cs_pairs(!prediction & !driver),
                 FN = .x %>% condition_n_gene_x_cs_pairs(!prediction & driver))
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(p = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                Precision = TP / (TP + FP),
                Recall = TP / (TP + FN),
                Sensitivity = TP / (TP + FN),
                Specificity = TN / (TN + FP),
                Fscore = (Precision * Recall) / (Precision + Recall)) %>%
  dplyr::ungroup()
write_tibble(performance$summary, 
             filename = paste0(outDir, "performance.tsv"))

# format for PR function input
PR_in <- predictions %>%
  dplyr::select(prediction_type, prediction_method, prediction, driver) %>%
  dplyr::group_by(prediction_method, prediction_type) %>%
  # refactor driver predictions for PR function
  dplyr::mutate(driver = ifelse(driver, "positive", "negative") %>% factor(c("positive", "negative")))

performance$PR <- PR_in %>%
  # calculate PR curve
  yardstick::pr_curve(driver, prediction) %>%
  dplyr::left_join(PR_in %>%
                     # calculate AUPRC
                     yardstick::pr_auc(driver, prediction) %>%
                     dplyr::select(prediction_method,
                                   prediction_type,
                                   PR_AUC = .estimate),
                   by = c("prediction_method", "prediction_type")) %>%
  dplyr::ungroup()

# Add area under curve metric to summary
performance$summary <- performance$summary %>%
  dplyr::left_join(performance$PR %>%
                     dplyr::distinct(prediction_method,
                                     prediction_type,
                                     PR_AUC),
                   by = c("prediction_method", "prediction_type"))

# generate PR plot ====
PR <- performance %>%
  purrr::map(
    ~ .x %>%
      dplyr::group_by(prediction_method) %>%
      dplyr::mutate(prediction_method = prediction_method %>%
                      paste0(" (",
                             PR_AUC %>% max %>% round(2),
                             ")"))
  ) %>%
  plot_PR +
  ggplot2::ggtitle(
    paste0(
      "Precision-Recall\n",
      "Trait = ", trait, "\n",
      "Celltypes = " , celltypes)
  ) +
  ggrepel::geom_text_repel(data = . %>% filter(prediction_type == "max"))

# generate AUPRC plot ====
AUPRC <- performance %>%
  plot_AUPRC +
  ggplot2::ggtitle(
    paste0(
      "Area Under Precision-Recall Curve\n",
      "Trait = ", trait, "\n",
      "Celltypes = " , celltypes)
  )

# save plots
pdf(paste0(outDir, "performance.pdf"), onefile = T)
print(PR)
print(AUPRC)
dev.off()

# # save plots to list
# gridExtra::grid.arrange(
#   PR, 
#   AUPRC + ggplot2::theme(legend.position = "none"), 
#   ncol = 2)
# PR_plots[[trait]][[celltypes]] <- recordPlot()


