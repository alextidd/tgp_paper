library(tidyverse) ; library(ggplot2) ; library(purrr)
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

trait <- "BC_Michailidou2017_FM" 
celltypes <- "BRST.MCF7.CNCR_celltype"
cat(trait, celltypes, "#############\n")
variant_to_gene_max_distance = 2e6

base_dir<- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/compare_methods/" ; setwd(wkdir)
out_dir <- paste("output", trait, celltypes, "", sep = "/")

# get trait files ====
known_genes_file <- paste0(base_dir, "wrangle_package_data/traits/output/", trait, "/known_genes.txt")
known_genes <- read_tibble(known_genes_file)$V1 %>%
  check_known_genes(known_genes_file)
variants_file <- paste0(base_dir, "wrangle_package_data/traits/output/", trait, "/variants.tsv")
variants <- read_tibble(variants_file, header = T) %>% rename(cs = CredibleSet)

# get gene predictions ====
vxg <- read_tibble(
  gzfile(paste(out_dir, "predictions.tsv", sep = "/")),
  header = T)

# gene master list - protein coding genes annotated, known genes annotated ====
g_universe <- variants %>%
  get_vxt_master(TSSs, variant_to_gene_max_distance) %>%
  mutate(known_gene = symbol %in% known_genes$symbol) %>%
  # protein coding, only variants with kg within 2Mb
  get_testable %>%
  distinct(symbol, known_gene) %>%
  mutate(g_universe = paste0("Michailidou2017_FM_variants_", variant_to_gene_max_distance/1000000, "Mb")) %>% 
  crossing(tibble(method = unique(vxg$method)))
method_g_universe <- vxg %>%
  # all genes within 2Mb
  valr::bed_slop(both = variant_to_gene_max_distance,
                 genome = ChrSizes,
                 trim = T) %>%
  mutate(known_gene = symbol %in% known_genes$symbol) %>%
  get_testable %>%
  dplyr::distinct(method, symbol, known_gene) %>%
  mutate(g_universe = paste0("method_FM_variants_", variant_to_gene_max_distance/1000000, "Mb"))
g_master <- bind_rows(g_universe, method_g_universe)

vxg_predictions <- vxg %>%
  # generate CS names for methods that don't provide them, as a placeholder for grouping
  group_by(method) %>%
  mutate(cs = case_when(is.na(cs) ~ paste0("cs.", row_number()),
                        TRUE ~ cs)) %>%
  # add all genes near CSs near known genes per method, GENCODE pc genes only
  right_join(g_master) %>%
  mutate(known_gene = known_gene %>% replace_na(FALSE),
         score = score %>% replace_na(0)) %>%
  # max protein-coding gene predicted per CS
  group_by(method, cs, g_universe) %>%
  mutate(max = ( score == max(score) & score > 0 )) %>%
  group_by(method)

# stats
vxg_predictions %>%
  group_by(method, g_universe) %>%
  group_modify(
    ~ data.frame(True = .x %>% condition_n_genes(known_gene),
                 Positive = .x %>% condition_n_genes(max),
                 TP = .x %>% condition_n_genes(max & known_gene),
                 FP = .x %>% condition_n_genes(max & !known_gene),
                 FN = .x %>% condition_n_genes(!max & known_gene),
                 TN = .x %>% condition_n_genes(!max & !known_gene))
  ) %>%
  rowwise() %>%
  dplyr::mutate(p_value = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                Precision = TP / (TP + FP),
                Recall = TP / (TP + FN),
                Sensitivity = TP / (TP + FN),
                Specificity = TN / (TN + FP),
                Fscore = (Precision * Recall) / (Precision + Recall)) %>%
  ungroup()

# permutations === 
kgs <- unique(known_genes$symbol)
permutations <- list()
trials <- 1000000
gene_predictions %>% 
  group_by(method) %>%
  group_modify(~ tibble(P = .x %>% condition_n_genes(max),
                        TP = .x %>% condition_n_genes(max & known_gene),
                        T = .x %>% condition_n_genes(known_gene))
    ) %>%
  group_split %>%
  map(function(i.df){
    successful_perms <- map(1:trials, function(trial){
      tibble(trial = trial,
             TP = length(intersect(sample(preds$symbol, i.df$P), unique(kgs))))  
    }) %>% 
      reduce(bind_rows) %>%
      filter(TP > i.TP) %>%
      nrow
    out <- i.df %>% mutate(successful_perms = successful_perms,
                           trials = trials,
                           p_value = successful_perms / trials)
    return(out)
  }) %>% reduce(bind_rows) -> permutations[[as.character(trials)]]
permutations
write_tibble(permutations[["1000"]], paste0(out_dir, "permutations_1000.tsv"))
write_tibble(permutations[["1000000"]], paste0(out_dir, "permutations_1000000.tsv"))
permutations %>%
  names %>%
  lapply(function(n){
    write_tibble(permutations[[n]],
                 paste0(out_dir, "permutations_", names(permutations), ".tsv"))
  })


# plot settings
metadata <- read_tibble("output/metadata.tsv", header = T)
colours <- metadata$colour ; names(colours) <- metadata$method
theme_set(
  theme_classic() +
    theme(
      title = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15)))
recall_lab <- paste0("recall (True = ", n_distinct(known_genes$symbol), ")")
pr_settings <- list(
  coord_equal(),
  scale_colour_manual(values = colours, limits = force),
  scale_fill_manual(values = colours, limits = force),
  labs(x = recall_lab),
  xlim(c(0,1)),
  ylim(c(0,1)))

# # stats
# method_stats <- max_predictions %>%
#   group_by(method) %>%
#   group_modify(
#     ~ data.frame(True = .x %>% condition_n_genes(known_gene),
#                  Positive = .x %>% condition_n_genes(score > 0),
#                  TP = .x %>% condition_n_genes(score > 0 & known_gene),
#                  FP = .x %>% condition_n_genes(score > 0 & !known_gene),
#                  FN = .x %>% condition_n_genes(score == 0 & known_gene))
#   ) %>%
#   rowwise() %>%
#   mutate(Precision = TP / (TP + FP),
#          Recall = TP / (TP + FN),
#          F_score = (Precision * Recall) / (Precision + Recall)) %>%
#   ungroup() %>%
#   mutate(F_score = replace_na(F_score, 0)) -> method_stats
# write_tibble(method_stats, paste0(out_dir, "method_stats.tsv"))
# 
# tier1_PR <- method_stats %>%
#   ggplot(aes(x = Recall, y = Precision, colour = method, label = method)) +
#   geom_point() + 
#   pr_settings +
#   ggrepel::geom_text_repel()
# tier1_Fscore <- method_stats  %>%
#   arrange(F_score) %>%
#   ggplot(aes(x = reorder(method, F_score), 
#              y = F_score, 
#              fill = method)) +
#   geom_col() +
#   scale_fill_manual(values = colours,
#                     limits = force) +
#   coord_flip() +
#   theme(axis.title.y = element_blank())
# 
# # STEP 2 ====
# pr_in <- gene_predictions %>%
#   # add known genes
#   full_join(crossing(kg_master, tibble(method = unique(gene_predictions$method)))) %>%
#   mutate(
#     known_gene = ifelse(is.na(known_gene), "negative", "positive") %>% factor(c("positive", "negative")),
#     score = score %>% replace_na(0)
#   ) 
# pr <- pr_in %>%
#   yardstick::pr_curve(known_gene, score) 
# auprc <- pr_in %>% 
#   yardstick::pr_auc(known_gene, score)
# 
# tier2_PR <- pr %>%
#   filter(.threshold %ni% c(0, Inf, -Inf), recall != 1, precision != 0) %>%
#   ggplot(aes(x = recall, y = precision, colour = method)) +
#   geom_line() +
#   pr_settings
# tier2_AUPRC <- auprc %>%
#   mutate(AUPRC = .estimate) %>%
#   ggplot(aes(x = reorder(method, AUPRC),
#              y = AUPRC,
#              fill = method)) +
#   geom_col() +
#   scale_fill_manual(values = colours[names(colours) %in% method_stats$method]) +
#   coord_flip() +
#   theme(axis.title.y = element_blank())
# 
# # save plots
# out <- gridExtra::grid.arrange(
#   tier1_PR + 
#     theme(legend.position = "none") +
#     ggtitle("Tier 1 PR\nmax gene lists"),
#   tier1_Fscore +
#     ggtitle("Tier 1 F scores\nmax gene lists"),
#   tier2_PR +
#     theme(legend.position = "none") +
#     ggtitle("Tier 2 PR curves\nscored gene lists"),
#   tier2_AUPRC +
#     ggtitle("Tier 2 AUPRC\nscored gene lists"),
#   ncol = 2, nrow = 2)
# 
# {
#   ggsave(paste0(out_dir, "performance.pdf"), out)
#   ggsave(paste0(out_dir, "performance.png"), out)
# }
