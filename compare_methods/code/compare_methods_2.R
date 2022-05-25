run_variants <- commandArgs(trailingOnly = T)[1]
# for testing: # run_variants="BC_Michailidou2017_FM"
cat(run_variants, "#############\n") ; message(run_variants, " #############\n")

library(tidyverse) ; library(ggplot2) ; library(purrr) ; library(gridExtra) ; library(grid)
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

# variables ====
variant_to_gene_max_distance = 2e6
max_n_known_genes_per_CS = Inf

# directories ====
base_dir <- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- paste0(base_dir, "compare_methods/") ; setwd(wkdir)
out_dir <- paste0("output/", run_variants, "/")

# metadata ====
methods_metadata <- read_tibble("output/metadata.tsv", header = T)
traits_metadata <- read_tibble(paste0(base_dir, "/wrangle_package_data/traits/output/metadata.tsv"), header = T)

# functions ====
get_performance <- function(confusion_mat, 
                            run_variants){
  confusion_mat %>%
    rowwise() %>%
    dplyr::mutate(p_value = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FP),
                  Recall = TP / (TP + FN),
                  Sensitivity = TP / (TP + FN),
                  Specificity = TN / (TN + FP),
                  F_score = (Precision * Recall) / (Precision + Recall)) %>%
    ungroup() %>%
    mutate(variants = run_variants,
           method_name = method %>% gsub("\\_.*", "", .)) %>%
    select(variants, method_name, method, everything())
}

# get variants files ====
run_traits_metadata <- traits_metadata %>%
  filter(variants == run_variants)
known_genes <- read_tibble(run_traits_metadata$known_genes_file)$V1 %>%
  check_known_genes(run_traits_metadata$known_genes_file)
variants <- import_BED(run_traits_metadata$variants_file, metadata_cols = c("variant", "cs"))
vxt_master <- get_vxt_master(variants, TSSs, variant_to_gene_max_distance)

# get distance-based predictions ====
closest <- vxt_master %>%
  group_by(cs, symbol) %>%
  filter(distance == min(distance)) %>%
  transmute(
    variant, cs, symbol,
    score = case_when(distance == 0 ~ 1, TRUE ~ 1 / distance),
    method = "closest")

# get gene predictions ====
annotations_in <- paste(out_dir, "predictions.tsv", sep = "/") %>%
  read_tibble(header = T) %>%
  bind_rows(closest) %>%
  select(-variant)
annotations <- annotations_in %>%
  right_join(
    vxt_master %>%
      distinct(cs, symbol) %>%
      crossing(tibble(method = unique(annotations_in$method)))
  ) %>%
  mutate(score = replace_na(score, 0)) %>%
  pivot_wider(id_cols = c(cs, symbol),
              names_from = method,
              values_from = score,
              values_fn = max)

# get performance ====
performance <- annotations %>%
  get_PR(vxt_master, known_genes, pcENSGs, max_n_known_genes_per_CS) %>%
  map(~ .x %>% mutate(method = prediction_method %>% gsub("\\_.*", "", .)))

no_preds <- performance$summary %>% filter(Positive == 0)
if(nrow(no_preds) > 0){message(paste(no_preds$prediction_method,collapse=", "), "\n... has no positive predictions - check predictions file for errors!")}

# colours ====
colours_df <- methods_metadata %>%
  right_join(performance$summary %>% select(method, prediction_method))
if(nrow(filter(colours_df, is.na(colour)))>0){stop(filter(colours_df, is.na(colour))$method, " does not have an assigned colour in the methods metadata!")}
colours <- colours_df$colour
names(colours) <- colours_df$prediction_method

# plot performance ====
plot_settings <- list(
  scale_colour_manual(values = colours, limits = force),
  scale_fill_manual(values = colours, limits = force),
  theme_classic(),
  theme(
    title = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.position = "none"))

PR_plot <- performance %>%
  map(~ filter(.x, prediction_type == "max")) %>%
  plot_PR(colour = prediction_method) +
  ggrepel::geom_text_repel(data = . %>%
                             dplyr::filter(prediction_type == "max")) +
  geom_line(data = performance$PR %>% filter(prediction_type == "score",
                                             .threshold %ni% c(0, Inf, -Inf),
                                             recall != 1,
                                             precision != 0),
            alpha = 0.3) +
  labs(x = paste0("recall (True = ",
                unique(performance$summary$True),
                ")")) +
  plot_settings

performance_plot <- performance$summary %>%
  mutate(fsc = replace_na(F_score, 0)) %>%
  pivot_longer(cols = c(F_score, Precision, Recall),
                      names_to = "metric",
                      values_to = "performance") %>%
  ggplot(aes(x = reorder(prediction_method, fsc),
             y = performance,
             fill = prediction_method)) +
  geom_col() +
  facet_grid(~ metric,
                      scales = "free", space = "free_y") +
  coord_flip() +
  theme(axis.title = element_blank()) +
  plot_settings

out_plot <- grid.arrange(
  performance_plot,
  PR_plot,
  ncol = 2, nrow = 1,
  top = textGrob(paste0("Method comparison (max score per CS)",
                        "\nVariants = ", run_variants,
                        "\nmax n known genes per CS = ", max_n_known_genes_per_CS,
                        ", max distance = ", variant_to_gene_max_distance),
                 gp = gpar(fontsize = 20)))

# write output ====
{
  write_tibble(performance$summary, paste0(out_dir, "performance.tsv"))
  ggplot2::ggsave(paste0(out_dir, "performance.pdf"), out_plot, width = 15)
  ggplot2::ggsave(paste0(out_dir, "performance.png"), out_plot, width = 15)
}


# triplets ====
#
# # save all
# all_vxt_master[[run_variants]] <- vxt_master
# all_known_genes[[run_variants]] <- known_genes
# all_annotations[[run_variants]] <- annotations
# 
# pred <- all_annotations %>%
#   bind_rows(.id = "variants") %>%
#   # replace missin scores with 0
#   mutate(across(where(is.numeric), replace_na, 0)) %>%
#   # annotate known genes
#   left_join(all_known_genes %>%
#               bind_rows(.id = "variants") %>%
#               transmute(variants, symbol, known_gene = T)) %>%
#   mutate(known_gene = known_gene %>% replace_na(FALSE)) %>%
#   # get testable
#   filter(symbol %in% filter(TSSs, ensg %in% pcENSGs)$symbol) %>% 
#   group_by(variants, cs) %>%
#   group_modify(function(x,y){
#     bind_cols(x, tibble(n_known_genes = x %>% condition_n_genes(known_gene)))
#   }) %>%
#   filter(n_known_genes > 0,
#          n_known_genes <= max_n_known_genes_per_CS) %>%
#   # get maximum score per CS-x-gene-x-method
#   group_by(variants, cs, symbol, known_gene) %>%
#   summarise(across(where(is.numeric), max)) %>%
#   # gather prediction methods
#   pivot_longer(
#     where(is.numeric),
#     names_to = "prediction_method",
#     values_to = "score"
#   ) %>%
#   # max prediction
#   group_by(variants, prediction_method, cs) %>%
#   mutate(max = as.numeric(score == max(score) & score > 0)) %>%
#   # gather prediction types
#   pivot_longer(
#     c(score, max), 
#     names_to = "prediction_type",
#     values_to = "prediction"
#   ) %>%
#   ungroup()
# PR_in <- pred %>%
#   # format for PR function input
#   select(prediction_type, prediction_method, prediction, known_gene) %>%
#   group_by(prediction_method, prediction_type) %>%
#   # refactor known gene predictions for PR function
#   mutate(known_gene = ifelse(known_gene, "positive", "negative") %>% factor(c("positive", "negative")))
# 
# performance$PR <- PR_in %>%
#   # calculate PR curve
#   yardstick::pr_curve(known_gene, prediction) %>%
#   left_join(PR_in %>%
#                      # calculate AUPRC
#                      yardstick::pr_auc(known_gene, prediction) %>%
#                      select(prediction_method,
#                                    prediction_type,
#                                    PR_AUC = .estimate),
#                    by = c("prediction_method", "prediction_type")) %>%
#   ungroup()
# 
# # get summary statistics (various performance metrics)
# performance$summary <- pred %>%
#   # score is not binary, cannot be summarised, filter out
#   filter(!grepl("score", prediction_type)) %>%
#   mutate(prediction = as.logical(prediction)) %>%
#   group_by(prediction_method, prediction_type) %>%
#   group_modify(
#     ~ data.frame(True = .x %>% condition_n_gene_x_cs_pairs(known_gene),
#                  Positive = .x %>% condition_n_gene_x_cs_pairs(prediction),
#                  TP = .x %>% condition_n_gene_x_cs_pairs(prediction & known_gene),
#                  FP = .x %>% condition_n_gene_x_cs_pairs(prediction & !known_gene),
#                  TN = .x %>% condition_n_gene_x_cs_pairs(!prediction & !known_gene),
#                  FN = .x %>% condition_n_gene_x_cs_pairs(!prediction & known_gene),
#                  n_known_genes = dplyr::filter(.x, known_gene)$symbol %>% dplyr::n_distinct())
#   ) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(p_value = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
#                 OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
#                 Precision = TP / (TP + FP),
#                 Recall = TP / (TP + FN),
#                 Sensitivity = TP / (TP + FN),
#                 Specificity = TN / (TN + FP),
#                 F_score = (Precision * Recall) / (Precision + Recall)) %>%
#   dplyr::ungroup() %>%
#   # add area under curve metric to summary
#   dplyr::left_join(performance$PR %>%
#                      dplyr::distinct(prediction_method,
#                                      prediction_type,
#                                      PR_AUC),
#                    by = c("prediction_method", "prediction_type"))
# 
# 
# 
# all_annots <- all_annotations %>%
#   bind_rows(.id = "variants") %>%
#   mutate(across(where(is.numeric), replace_na, 0),
#          cs = paste0(variants, "|", cs))
# 
# all_kgs <- rbind(all_known_genes)
# 
# performance <- annotations %>%
#   get_PR(vxt_master, known_genes, pcENSGs, max_n_known_genes_per_CS) %>%
#   map(~ .x %>% mutate(method = prediction_method %>% gsub("\\_.*", "", .)))
# 
# all_vars <- all_variants %>%
#   bind_rows(.id = "variants") 
# # make sure lack of prediction for a trait by a method is filled out
# all_trues <- all_vars %>%
#   distinct(variants, True) %>%
#   crossing(tibble(prediction_method = unique(all_performance$prediction_method)))
# all_performance <- all_vars %>%
#   right_join(all_trues) %>%
#   mutate(across(where(is.numeric), replace_na, 0)) %>%
#   # sum performance
#   group_by(prediction_method) %>%
#   summarise(across(where(is.numeric), sum))  %>%
#   rowwise() %>%
#   dplyr::mutate(p_value = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
#                 OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
#                 Precision = TP / (TP + FP),
#                 Recall = TP / (TP + FN),
#                 Sensitivity = TP / (TP + FN),
#                 Specificity = TN / (TN + FP),
#                 F_score = (Precision * Recall) / (Precision + Recall)) %>%
#   ungroup() %>%
#   mutate(method = prediction_method %>% gsub("\\_.*", "", .)) 
# 
# # colours ====
# colours_df <- methods_metadata %>% 
#   left_join(all_performance %>% select(method, prediction_method)) %>%
#   filter(!is.na(prediction_method))
# colours <- colours_df$colour
# names(colours) <- colours_df$prediction_method
# 
# # plot ====
# performance_plot <- all_performance %>%
#   mutate(fsc = F_score) %>%
#   pivot_longer(cols = c(F_score, Precision, Recall),
#                names_to = "metric",
#                values_to = "performance") %>%
#   ggplot(aes(x = reorder(prediction_method, fsc),
#              y = performance,
#              fill = prediction_method)) +
#   geom_col() +
#   facet_grid(~ metric,
#              scales = "free", space = "free_y") +
#   coord_flip() +
#   theme(axis.title = element_blank()) +
#   scale_fill_manual(values = colours, limits = force) +
#   theme_classic()
#   
# # write output ====
# write_tibble(all_performance, 
#              paste0("output/all_variants/performance.tsv"))
# {
#   ggplot2::ggsave(paste0("output/all_variants/performance.pdf"), performance_plot, width = 15)
#   ggplot2::ggsave(paste0("output/all_variants/performance.png"), performance_plot, width = 15)
# }
