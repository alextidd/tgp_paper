library(tidyverse) ; library(ggplot2) ; library(purrr) ; library(gridExtra) ; library(grid)
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

# variables ====
variant_to_gene_max_distance = 2e6
max_n_known_genes_per_CS = Inf
base_dir <- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- paste0(base_dir, "compare_methods/") ; setwd(wkdir)

# functions ====
get_group_conusion_mat <- function(cxg_predictions){
  cxg_predictions %>%
  group_modify(
    ~ data.frame(True = .x %>% condition_n_gene_x_cs_pairs(known_gene),
                 Positive = .x %>% condition_n_gene_x_cs_pairs(max),
                 TP = .x %>% condition_n_gene_x_cs_pairs(max & known_gene),
                 FP = .x %>% condition_n_gene_x_cs_pairs(max & !known_gene),
                 FN = .x %>% condition_n_gene_x_cs_pairs(!max & known_gene),
                 TN = .x %>% condition_n_gene_x_cs_pairs(!max & !known_gene))
  ) 
}
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
plot_performance <- function(metadata,
                             performance,
                             variants_name,
                             max_n_known_genes_per_CS,
                             variant_to_gene_max_distance){
  colours <- metadata$colour ; names(colours) <- metadata$method
  plot_settings <- list(
    scale_colour_manual(values = colours, limits = force),
    scale_fill_manual(values = colours, limits = force),
    labs(x = paste0("recall (True = ", 
                    unique(performance$True), 
                    ")")),
    theme_classic(),
    theme(
      title = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15)))
  plot_subtitle <- paste0(
    "Variants = ", variants_name,
    "\nmax n known genes per CS = ", max_n_known_genes_per_CS, 
    ", max distance = ", variant_to_gene_max_distance)
  
  # plots
  performance_plot <- performance %>%
    mutate(fsc = F_score) %>%
    pivot_longer(cols = c(F_score, Precision, Recall),
                 names_to = "metric",
                 values_to = "performance") %>%
    ggplot(aes(x = reorder(method, fsc),
               y = performance,
               fill = method_name)) +
    geom_col() +
    facet_grid(~ metric) +
    coord_flip() +
    plot_settings +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    ggtitle("Performance metrics")
  
  PR_plot <- performance %>% 
    ggplot(aes(x = Recall, y = Precision, colour = method_name, label = method)) +
    geom_point() +
    plot_settings +
    ggrepel::geom_text_repel(max.overlaps = n_distinct(performance$method)) +
    xlim(c(0,1)) + ylim(c(0,1)) + coord_equal() + 
    ggtitle("Precision-Recall curve")
  
  out <- grid.arrange(
    performance_plot,
    PR_plot,
    ncol = 2, nrow = 1,
    top = textGrob(paste0("Method comparison (max score per CS)\n",
                          plot_subtitle),
                   gp = gpar(fontsize = 20)))
  return(out)
}

all_variants <- list() ; for(variants_name in 
    unique(read_tibble(paste0(base_dir, "/wrangle_package_data/traits/output/metadata.tsv"), header = T)$variants)
    ){

  cat(variants_name, "#############\n")
  out_dir <- paste0("output/", variants_name, "/")
    
  # get methods metadata ====
  metadata <- read_tibble("output/metadata.tsv", header = T)
  
  # get variants files ====
  variants_metadata <- paste0(base_dir, "/wrangle_package_data/traits/output/metadata.tsv") %>%
    read_tibble(header = T) %>%
    filter(variants == variants_name)
  known_genes <- read_tibble(variants_metadata$known_genes_file)$V1 %>%
    check_known_genes(variants_metadata$known_genes_file)
  variants <- import_BED(variants_metadata$variants_file, metadata_cols = c("variant", "cs"))
  
  # get gene predictions ====
  vxg <- read_tibble(
    paste(out_dir, "predictions.tsv", sep = "/"),
    header = T) 
  
  # gene master list - protein coding genes annotated, nearby known genes annotated ====
  cxg_universe <- variants %>%
    get_vxt_master(TSSs, variant_to_gene_max_distance) %>%
    mutate(known_gene = symbol %in% known_genes$symbol) %>%
    # protein coding, only variants with kg within 2Mb
    get_testable %>%
    distinct(cs, symbol, known_gene) %>%
    crossing(tibble(method = unique(vxg$method)))
  
  # cxg predictions in the cxg universe of the analysis
  cxg_predictions <- vxg %>%
    select(cs, symbol, score, method) %>%
    group_by(cs, method) %>%
    filter(score == max(score) & score > 0) %>%
    mutate(max = T) %>%
    # add all genes near CSs near known genes per method, GENCODE pc genes only
    right_join(cxg_universe, by = c("cs", "symbol", "method")) %>%
    mutate(across(where(is.logical), ~ replace_na(.x, FALSE)))
  
  # performance ====
  confusion_mat <- cxg_predictions %>%
    group_by(method) %>%
    get_group_conusion_mat
  performance <- get_performance(confusion_mat, variants_name)
  
  all_variants[[variants_name]] <- confusion_mat
  
  # plot_performance
  performance_plot <- plot_performance(metadata, 
                                       performance, 
                                       variants_name, 
                                       max_n_known_genes_per_CS,
                                       variant_to_gene_max_distance)
  
  # write output
  write_tibble(performance, 
               paste0(out_dir, "performance.tsv"))
  {
    ggplot2::ggsave(paste0(out_dir, "performance.pdf"), performance_plot, width = 15)
    ggplot2::ggsave(paste0(out_dir, "performance.png"), performance_plot, width = 15)
  }
  
}

# triplets ====
all_performance <- all_variants %>%
  bind_rows(.id = "variants") %>%
  group_by(method) %>%
  summarise(across(where(is.numeric), sum)) %>%
  get_performance(variants_name = "all_variants")
all_performance_plot <- plot_performance(metadata, all_performance, variants_name = "all_variants",
                 max_n_known_genes_per_CS, variant_to_gene_max_distance)
# write output
write_tibble(all_performance, 
             paste0("output/all_variants/performance.tsv"))
{
  ggplot2::ggsave(paste0("output/all_variants/performance.pdf"), performance_plot, width = 15)
  ggplot2::ggsave(paste0("output/all_variants/performance.png"), performance_plot, width = 15)
}