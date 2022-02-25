# L2G
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/compare_methods/" ; setwd(wkdir)
allVarsDir <- "/working/lab_georgiat/jonathB/DATA/GWAS/BCAC/GRCH38/chr_pos_hg19_GRCH38_ID_Pval_beta_SE/overall/"
L2GPredDir <- "output/BC_expanded/methods/L2G/"
window = 5e5

# packages
library(tidyverse)
library(ggplot2)
library(ggrepel)

# load linked w and wo CSs
w_CSs <- read.delim(paste0(L2GPredDir, "L2G_w_CS.tsv")) %>% tibble
wo_CSs <- read.delim(paste0(L2GPredDir, "L2G_wo_CS.tsv")) %>% tibble %>%
  separate(hg38_variant, "chrom", extra = "drop", remove = F) %>%
  mutate(chrom = paste0("chr", chrom))

# get CSs (join by rsID)
allBCVars <- read.delim(paste0(L2GPredDir, "BCACFM.CCRV.STRONG.AND.SECONDARY.variants.with_hg38_var.bed"), header =F) %>% tibble %>%
  rename(CredibleSet = V1, SignalStrength = V2, hg38_variant = V3) %>%
  # remove CIMBAs
  filter(!grepl("CIMBA", CredibleSet)) %>%
  mutate(BCAC_FM = CredibleSet %>% gsub("ichav.*", "", .),
         ichav = CredibleSet)

#### PLOTTING ####

# Plot window
plots <- list()
for(chrom in unique(wo_CSs$chrom)){
  print(chrom)
  chrom.allVars <- allVarsDir %>%
    paste0(chrom, ".txt.gz") %>%
    read.delim(header = F) %>% tibble %>%
    separate(V2, c("chr", "position"), extra = "drop", remove = F) %>%
    transmute(chrom = V1, start = as.numeric(position), end = as.numeric(position),
              hg38_variant = V3, hg19_variant = V2,
              rsID = V4, p = V7)
  # get lead BCAC SNPs
  chrom.BCVars <- allBCVars %>%
    inner_join(chrom.allVars) %>%
    group_by(BCAC_FM, ichav) %>%
    mutate(BCAC_FM_lead_variant = p == min(p),
           group = ichav) %>%
    ungroup()
  for(var in wo_CSs$hg38_variant[wo_CSs$chrom == chrom]){
    print(var)
    i.var <- wo_CSs %>%
      filter(hg38_variant == var) %>%
      left_join(chrom.allVars) %>%
      mutate(group = "L2G lead variant")
    if(nrow(i.var) > 1){
      i.var <- i.var %>%
        filter(grepl(wo_CSs %>%
                       filter(hg38_variant == var) %>%
                       pull(hg38_variant_og) %>%
                       gsub(".*[0-9]\\_.*[0-9]\\_", "", .),
                     hg19_variant))
    }
    i.allVars <- chrom.allVars %>%
      filter(start > i.var$start - window, end < i.var$end + window)
    i.vars_w_CS <- w_CSs %>%
      inner_join(i.allVars) %>%
      mutate(group = "L2G lead variant\n(signal already found)")
    i.BCVars <- chrom.BCVars %>%
      inner_join(i.allVars)
    
    # all
    i.all <- i.var %>%
      full_join(i.BCVars) %>%
      full_join(i.vars_w_CS) %>%
      full_join(i.allVars) %>%
      mutate(group = replace_na(group, "background variant"))
    
    # get colours
    fixed_colours <- c("L2G lead variant" = "cornflowerblue",
                       "L2G lead variant\n(signal already found)" = "darkblue",
                       "background variant" = "grey")
    n_colours <- length(setdiff(unique(i.all$group), names(fixed_colours)))
    manual_colours <- c("lightsalmon", "lightgoldenrod1", "lightpink", "lightgreen", "lavender", "burlywood", "slategray1", "orange", "cyan") %>%
      head(n_colours) %>%
      setNames(setdiff(unique(i.all$group), names(fixed_colours))) %>%
      c(fixed_colours, .)
    manual_colours <- manual_colours[names(manual_colours) %in% unique(i.all$group)]
    
    # get CSs involved
    i.CSs <- filter(allBCVars, BCAC_FM %in% i.BCVars$BCAC_FM) %>%
      inner_join(chrom.allVars) %>%
      group_by(BCAC_FM) %>%
      transmute(chrom = chrom, start = min(start), end = max(end),
                BCAC_FM,
                title = paste0(BCAC_FM, " (", chrom,":", start, "-", end, ")\n")) %>%
      distinct
    
    # plot
    plots[[var]] <- i.all %>%
      ggplot(aes(x = end, y = -log10(p), colour = group, label = rsID)) + #, shape = BCAC_FM_lead_variant)) +
      scale_colour_manual(values = manual_colours) +
      geom_point(data = . %>% filter(group == "background variant"), alpha = 0.5) +
      geom_point(data = . %>% filter(grepl("ichav", group)), size = 3) +
      geom_point(data = . %>% filter(group == "L2G lead variant"), size = 3) +
      geom_point(data = . %>% filter(BCAC_FM_lead_variant), size = 3, shape = 1, colour = "black") +
      geom_point(data = . %>% filter(group == "L2G lead variant\n(signal already found)"), shape = 17) +
      geom_text_repel(data = . %>% arrange(p) %>% head(5)) +
      labs(title = paste0("L2G: ", var, " (", paste0(i.var$chrom, ":", i.var$end), ")\n",
                          "L2G Target Gene: ", i.var$TargetGene, "\n",
                          "BCAC: ", paste0(i.CSs$title, collapse = "BCAC: ")),
           x = "position (bp)") +
      theme_bw()
  }
}

# save plots
pdf(paste0(L2GPredDir, "/nomatch_variants_Manhattan.pdf"))
walk(plots, print)
dev.off()

# LDlinkR for leftover lead SNPs (still no clear ichav assigned from manual inspection)
leftover <- read.delim(paste0(L2GPredDir, "L2G_to_LDlinkR.tsv")) %>% tibble %>%
  separate(variant_fixed, "chrom", remove = F, extra = "drop") %>%
  mutate(chrom = paste0("chr", chrom), hg38_variant = variant_fixed)

# get all CCVs within 500kb of each leftover L2G lead SNP -> for LDlinkR
to_LDlinkR <- tibble()
for(chrom in unique(leftover$chrom)){
  print(chrom)
  chrom.allVars <- allVarsDir %>%
    paste0(chrom, ".txt.gz") %>%
    read.delim(header = F) %>% tibble %>%
    separate(V2, c("chr", "position"), extra = "drop", remove = F) %>%
    transmute(chrom = V1, start = as.numeric(position), end = as.numeric(position),
              hg38_variant = V3, hg19_variant = V2,
              rsID = V4, p = V7)
  # get lead BC SNPs
  chrom.BCVars <- allBCVars %>%
    inner_join(chrom.allVars) %>%
    group_by(BCAC_FM, ichav) %>%
    mutate(BCAC_FM_lead_variant = p == min(p),
           group = ichav) %>%
    ungroup()
  for(var in leftover$hg38_variant[leftover$chrom == chrom]){
    i.var <- leftover %>%
      filter(hg38_variant == var) %>%
      left_join(chrom.allVars) %>%
      mutate(group = "L2G lead variant")
    if(nrow(i.var) > 1){
      i.var <- i.var %>%
        filter(grepl(wo_CSs %>%
                       filter(hg38_variant == var) %>%
                       pull(hg38_variant_og) %>%
                       gsub(".*[0-9]\\_.*[0-9]\\_", "", .),
                     hg19_variant))
    }
    # get all CCVs within window
    i.BCVars <- chrom.BCVars %>%
      filter(start > i.var$start - window, end < i.var$end + window)
    
    to_LDlinkR <- bind_rows(i.var %>% transmute(rsID, hg38_variant, SNP = "L2G"),
                            i.BCVars %>% transmute(rsID, hg38_variant, SNP = "BCAC")) %>%
      mutate(L2G_lead = i.var$hg38_variant) %>%
      bind_rows(., to_LDlinkR)
  }
}
to_LDlinkR %>%
  write.table(paste0(L2GPredDir, "L2G_to_LDlinkR_rsIDs.tsv"),
              quote = F, sep = "\t", row.names = F)

# run on hpcapp01 head node! ###########################################################
# # LDlinkR #
# wkdir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/" ; setwd(wkdir)
# L2GPredDir <- "output/Predictions/L2G/"
# library(LDlinkR)
# library(dplyr)
# to_LDlinkR <- read.delim(paste0(L2GPredDir, "L2G_to_LDlinkR_rsIDs.tsv"))
# proxies <- data.frame()
# for(i.L2G_lead in unique(to_LDlinkR$L2G_lead)){
#   print(i.L2G_lead)
#   i.to_LDlinkR <- to_LDlinkR %>% filter(L2G_lead == i.L2G_lead)
#   i.lead_rsID <- filter(i.to_LDlinkR, SNP == "L2G")$rsID
#   if(nrow(i.to_LDlinkR) > 1){
#     mat <- LDmatrix(i.to_LDlinkR$rsID,
#              pop = "CEU",
#              r2d = "r2",
#              token = "c34f6899a9d3")
#     if(i.lead_rsID %in% colnames(mat)){
#       out <- mat %>%
#         select(BCAC_proxy_rsID = RS_number,
#                r2 = filter(i.to_LDlinkR, SNP == "L2G")$rsID) %>%
#         left_join(i.to_LDlinkR %>% select(BCAC_proxy_rsID = rsID,
#                                           BCAC_proxy = hg38_variant)) %>%
#         transmute(L2G_lead = i.L2G_lead,
#                   BCAC_proxy,
#                   r2) %>%
#         filter(L2G_lead != BCAC_proxy) %>%
#         filter(r2 == max(r2, na.rm = T))
#       proxies <- bind_rows(proxies, out)
#     } else { cat("rsID not found\n") }
#   } else { cat("no CSs nearby\n") }
# }
# proxies %>%
#   write.table(paste0(L2GPredDir, "L2G_to_LDlinkR_proxies.tsv"),
#               quote = F, sep = "\t", row.names = F)
#####################################################################################

# get proxies
proxies <- read.delim(paste0(L2GPredDir, "L2G_to_LDlinkR_proxies.tsv")) %>% tibble
proxies %>%
  left_join(allBCVars, by = c("BCAC_proxy" = "hg38_variant")) %>%
  mutate(ichav = gsub(".*\\_", "", ichav)) %>%
  distinct(L2G_lead, BCAC_FM, ichav, r2) %>%
  write.table(paste0(L2GPredDir, "L2G_to_LDlinkR_best_proxies.tsv"),
              quote = F, sep = "\t", row.names = F)

# final intersected predictions (from Google Sheets L2G_cleanup)
preds <- read.delim(paste0(L2GPredDir, "L2G_AllPredictions.tsv")) %>% tibble
preds[preds==""] <- NA
preds <- preds %>%
  filter(!is.na(BCAC_FM),
         !is.na(ichav),
         !is.na(L2G)) %>%
  transmute(TargetGene = L2G, CredibleSet = paste0(BCAC_FM, ichav))
write.table(preds, "output/BC_expanded/methods/L2G/IntersectedPredictions.tsv",
            row.names = F, quote = F, sep = "\t")
system("gzip -f output/BC_expanded/methods/L2G/IntersectedPredictions.tsv")


# preds <- preds %>%
#   filter(!is.na(BCAC_FM),
#          !is.na(ichav),
#          !is.na(L2G)) %>%
#   separate(hg19_variant, c("chrom", "position"), extra = "drop", remove = F) %>%
#   transmute(chrom = paste0("chr", chrom), start = as.numeric(position), end = as.numeric(position),
#             TargetGene = L2G, CredibleSet = paste0(BCAC_FM, ichav))