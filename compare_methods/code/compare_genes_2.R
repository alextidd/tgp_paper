# get gene predictions ====
library(tidyverse) ; library(ggplot2)

trait <- "BC_Michailidou2017_FM" 
celltypes <- "BRST.MCF7.CNCR_celltype"

base_dir<- "/working/lab_jonathb/alexandT/tgp_paper/"
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/compare_methods/" ; setwd(wkdir)
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

cat(trait, celltypes, "#############\n")
out_dir <- paste("output", trait, celltypes, "", sep = "/")

# get methods metadata ====
metadata <- read_tibble("output/metadata.tsv", header = T)
colours <- metadata$colour ; names(colours) <- metadata$method

# get trait known_genes / variants ====
known_genes_file <- paste0(base_dir, "wrangle_package_data/traits/output/", trait, "/known_genes.txt")
known_genes <- read_tibble(known_genes_file)$V1 %>%
  check_known_genes(known_genes_file)

# master list - protein coding genes, known genes annotated ====
g_master <- TSSs %>% 
  distinct(symbol, ensg) %>%
  filter(ensg %in% pcENSGs) %>%
  transmute(
    symbol,
    known_gene = symbol %in% known_genes$symbol
    ) %>%
  distinct()
kg_master <- g_master %>% 
  filter(known_gene)

# get gene predictions ====
gene_predictions <- read_tibble(
  gzfile(paste(out_dir, "genes.tsv", sep = "/")),
  header = T
) %>% 
  # protein-coding
  filter(symbol %in% g_master$symbol) %>%
  # remove CIMBA
  filter(!grepl("CIMBA", cs)) %>%
  # generate CS names for methods that don't provide them, as a placeholder for grouping
  group_by(method) %>%
  mutate(cs = case_when(is.na(cs) ~ paste0("cs.", row_number()),
                        TRUE ~ cs))

# STEP 1 ====
gene_predictions %>%  
  # max protein-coding gene predicted per CS
  group_by(method, cs) %>%
  filter(score == max(score)) %>% 
  # add known genes
  full_join(crossing(kg_master, tibble(method = unique(gene_predictions$method)))) %>%
  mutate(known_gene = known_gene %>% replace_na(FALSE),
         score = score %>% replace_na(0)) %>%
  # stats
  group_by(method) %>%
  group_modify(
    ~ data.frame(True = .x %>% condition_n_genes(known_gene),
                 Positive = .x %>% condition_n_genes(score > 0),
                 TP = .x %>% condition_n_genes(score > 0 & known_gene),
                 FP = .x %>% condition_n_genes(score > 0 & !known_gene),
                 FN = .x %>% condition_n_genes(score == 0 & known_gene))
  ) %>%
  rowwise() %>%
  mutate(Precision = TP / (TP + FP),
         Recall = TP / (TP + FN),
         F_score = (Precision * Recall) / (Precision + Recall)) %>%
  ungroup() -> method_stats

method_stats %>%
  ggplot(aes(x = Recall, y = Precision, colour = method)) +
  geom_point() + 
  scale_colour_manual(values = colours) +
  coord_equal() 

# STEP 2 ====
pr_in <- gene_predictions %>%
  # add known genes
  full_join(crossing(kg_master, tibble(method = unique(gene_predictions$method)))) %>%
  mutate(
    known_gene = ifelse(is.na(known_gene), "negative", "positive") %>% factor(c("positive", "negative")),
    score = score %>% replace_na(0)
  ) 
pr <- pr_in %>%
  yardstick::pr_curve(known_gene, score) 
auprc <- pr_in %>% 
  yardstick::pr_auc(known_gene, score)

pr %>%
  filter(.threshold %ni% c(0, Inf, -Inf), recall != 1, precision != 0) %>%
  ggplot(aes(x = recall, y = precision, colour = method)) +
  geom_line() +
  scale_colour_manual(values = colours) +
  coord_equal() +
  theme_bw()
