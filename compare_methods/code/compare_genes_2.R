library(tidyverse) ; library(ggplot2) ; library(purrr) ; library(gridExtra) ; library(grid)
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

# variables ====
variant_to_gene_max_distance = 2e6
max_n_known_genes_per_CS = Inf
base_dir <- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- paste0(base_dir, "compare_methods/") ; setwd(wkdir)

# functions ====
get_performance <- function(confusion_mat, 
                            variants_name){
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
    mutate(variants = variants_name,
           method_name = method %>% gsub("\\_.*", "", .)) %>%
    select(variants, method_name, method, everything())
}

# loop ====
all_variants <- list() ; for(variants_name in 
    unique(read_tibble(paste0(base_dir, "/wrangle_package_data/traits/output/metadata.tsv"), header = T)$variants)
    ){ # for testing: # all_variants<-list();variants_name="BC_Michailidou2017_FM"

  cat(variants_name, "#############\n")
  out_dir <- paste0("output/", variants_name, "/")
    
  # get methods metadata ====
  methods_metadata <- read_tibble("output/metadata.tsv", header = T)
  
  # get variants files ====
  variants_metadata <- paste0(base_dir, "/wrangle_package_data/traits/output/metadata.tsv") %>%
    read_tibble(header = T) %>%
    filter(variants == variants_name)
  known_genes <- read_tibble(variants_metadata$known_genes_file)$V1 %>%
    check_known_genes(variants_metadata$known_genes_file)
  variants <- import_BED(variants_metadata$variants_file, metadata_cols = c("variant", "cs"))
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
  all_variants[[variants_name]] <- performance$summary
  
  no_preds <- performance$summary %>% filter(Positive == 0)
  if(nrow(no_preds) > 0){message(paste(no_preds$prediction_method,collapse=", "), "\n... has no positive predictions - check predictions file for errors!")}

  # plot performance ====
  colours_df <- methods_metadata %>% 
    left_join(performance$summary %>% select(method, prediction_method)) %>%
    filter(!is.na(prediction_method))
  colours <- colours_df$colour
  names(colours) <- colours_df$prediction_method
  plot_settings <- list(
    scale_colour_manual(values = colours, limits = force),
    scale_fill_manual(values = colours, limits = force),
    labs(x = paste0("recall (True = ", 
                    unique(performance$summary$True), 
                    ")")),
    theme_classic(),
    theme(
      title = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      legend.position = "none"))
  plot_subtitle <- paste0(
    "Variants = ", variants_name,
    "\nmax n known genes per CS = ", max_n_known_genes_per_CS, 
    ", max distance = ", variant_to_gene_max_distance)
  
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
    plot_settings
  performance_plot <- performance$summary %>%
    mutate(fsc = F_score) %>%
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
    top = textGrob(paste0("Method comparison (max score per CS)\n",
                          plot_subtitle),
                   gp = gpar(fontsize = 20)))
  
  # write output ====
  { write_tibble(performance$summary, paste0(out_dir, "performance.tsv"))
    ggplot2::ggsave(paste0(out_dir, "performance.pdf"), out_plot, width = 15)
    ggplot2::ggsave(paste0(out_dir, "performance.png"), out_plot, width = 15) }
  
}

# triplets ====
all_performance <- all_variants %>%
  bind_rows(.id = "variants") %>%
  group_by(method) %>%
  summarise(across(where(is.numeric), sum)) %>%
  get_performance(variants_name = "all_variants")
all_performance_plot <- plot_performance(methods_metadata, all_performance, variants_name = "all_variants",
                 max_n_known_genes_per_CS, variant_to_gene_max_distance)
# write output
write_tibble(all_performance, 
             paste0("output/all_variants/performance.tsv"))
{
  ggplot2::ggsave(paste0("output/all_variants/performance.pdf"), performance_plot, width = 15)
  ggplot2::ggsave(paste0("output/all_variants/performance.png"), performance_plot, width = 15)
}