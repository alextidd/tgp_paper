wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/" ; setwd(wkdir) 
library(idr2d)
devtools::load_all("/working/lab_jonathb/alexandT/tgp/")
HiChIP_dir <- "output/HiChIP/"
replicates_dir <- paste0(HiChIP_dir, "replicates/")
byChr_dir <- paste0(HiChIP_dir, "by_chr/")
preQC_dir <- paste0(HiChIP_dir, "pre_QC/")
postQC_dir <- paste0(HiChIP_dir, "post_QC/") ; dir.create(postQC_dir, showWarnings = F)
interaction_max_distance = 2e6
HiChIP_metadata <- read_tibble("output/metadata.tsv", header = T) %>% 
  dplyr::filter(object == "HiChIP")

# Dealing with replicates ====
replicate_prefixes <- list.files(replicates_dir, pattern = "replicate.*bedpe") %>% sub("replicate.*", "", .) %>% unique

# # install (run on hpcapp01 head node)
# if (!requireNamespace("remotes", quietly = TRUE)) { install.packages("remotes") }
# remotes::install_github("kkrismer/idr2d")

# load data and run idr2d
if(length(replicate_prefixes)>0){
  for(rp in replicate_prefixes){
    cat(rp, "\n")
    files <- list.files(replicates_dir, pattern = rp, full.names = T)
    rep1_df <- read.delim(files[1], header = F)
    rep2_df <- read.delim(files[2], header = F)
    idr_results <- estimate_idr2d(rep1_df, rep2_df, value_transformation = "identity")
    summary(idr_results)
    # results
    rep1_idr_df <- idr_results$rep1_df
    # filter for low IDR
    out <- rep1_idr_df %>%
      dplyr::filter(idr < 0.05) %>%
      # value = -log10(idr) / max(-log10(idr)) (normalised to 1)
      dplyr::transmute(dplyr::across(chr_a:end_b),
                       value = (-log10(idr)/max(-log10(idr))))
    # write
    write.table(out,
                paste0(byChr_dir,
                       rp %>% sub("\\_$", "", .),
                       ".bedpe"),
                quote = F, sep = "\t", row.names = F, col.names = F)
  }
  # merge per-chromosome Corces files back together
  system(paste(paste0("HiChIP_dir=", wkdir, "output/HiChIP/"),
               "celltypes=$(ls $HiChIP_dir/by_chr/Corces*.bedpe | sed 's/.*Corces2020\\_//g' | sed 's/\\_.*//g' | sort -u)",
               "(for celltype in ${celltypes[@]} ; do echo $celltype",
               "outfile=$HiChIP_dir/pre_QC/Corces2020_${celltype}_HiChIP.bedpe ; > $outfile",
               "(for chrfile in $HiChIP_dir/by_chr/*${celltype}* ; do cat $chrfile >> $outfile ; done)",
               "done) ; ", sep = " ; "))
}

# QC

# files in metadata only
HiChIP <- list() 
for(ct in HiChIP_metadata$celltype){ 
  file <- paste0(preQC_dir, ct, ".bedpe")
  df <- import_BEDPE_to_List(file, metadata_cols = "score") %>%
    purrr::map(~ .x %>%
                 # shift starts +1 so that bins are all mutually exclusive
                 dplyr::mutate(start = start + 1) %>%
                 # for infinite score values, set equal to the maximum non-infinite score
                 dplyr::mutate(score = dplyr::case_when(is.infinite(score) ~ max(score[!is.infinite(score)]),
                                                        TRUE ~ score)) %>%
                 # scale to 1
                 dplyr::mutate(score = score/max(score)))
  
  # # find loops that both connect the same two regions (both ends intersect)
  # partial_dups <- df %>%
  #   purrr::map(~ bed_intersect_left(.x, .x, suffix = c("A", "B"), suffixMetadataCols = T) %>%
  #                dplyr::filter(InteractionIDA != InteractionIDB) %>%
  #                dplyr::rowwise() %>%
  #                dplyr::mutate(overlapping = paste(sort(c(InteractionIDA, InteractionIDB)), collapse = "_x_")) %>%
  #                dplyr::ungroup()) %>%
  #   purrr::reduce(inner_join, by = "overlapping") %>%
  #   dplyr::distinct(overlapping) %>%
  #   tidyr::separate(overlapping, c("InteractionID_A", "InteractionID_B"), sep = "\\_x\\_", remove = F) %>%
  #   tidyr::pivot_longer(c(InteractionID_A, InteractionID_B), values_to = "InteractionID") %>%
  #   dplyr::select(-name)
  
  # find duplicate, mirrored loop entries within each dataset and filter to exclude the non-maximum scores
  # find loops with > max distance
  test <- dplyr::full_join(df$first, df$last, by = c("InteractionID","score")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(loop_ends = paste(min(start.x, start.y), min(end.x, end.y), max(start.x, start.y), max(end.x, end.y), sep = ".")) %>%
    dplyr::group_by(loop_ends)
  dup_IDs <-  test %>%
    dplyr::filter(
      # if mirrored loops have the same score, exclude all but first 
      # if they have different scores, exclude all but maximum
      rank(score, ties.method = "first") != 1) %>%
    dplyr::pull(InteractionID)
  far_IDs <- test %>%
    dplyr::filter(
      # find loops with > max distance
      min(abs(end.x - start.y), abs(end.y - start.x)) > interaction_max_distance) %>%
    dplyr::pull(InteractionID)
  
  # duplicate / far loops message
  if(length(dup_IDs) > 0){
    message(length(dup_IDs), " / ", dplyr::n_distinct(test$InteractionID),
            " loops are duplicates. Each group of loops with identical ends is filtered to only include the maximum-scoring loop.")
  }
  if(length(far_IDs) > 0){
    message(length(far_IDs), " / ", dplyr::n_distinct(test$InteractionID),
            " loops are > ", interaction_max_distance, "bp apart and will be excluded.")
  }
  
  # filter
  df <- df %>%
    purrr::map(~ .x %>%
                 # exclude lower-scoring duplicate loops and far loops
                 dplyr::filter(InteractionID %ni% c(dup_IDs, far_IDs)) %>%
                 # filter score > median(score), unless median(score) == max(score) (ie binary scoring)
                 { if (median(.$score) < max(.$score)) dplyr::filter(., score > median(score)) else . }
               )

  # quantile normalisation and decile binning
  df <- df %>%
    lapply(function(x){
      x %>%
        dplyr::select(score) %>%
        data.matrix %>%
        preprocessCore::normalize.quantiles() %>%
        tibble::as_tibble() %>%
        dplyr::rename(score = V1) %>%
        dplyr::bind_cols(x %>% dplyr::select(-score), .) # %>%
        # # split into deciles, unless only one score (binary loops)
        # dplyr::mutate(score = dplyr::case_when(
        #   dplyr::n_distinct(x$score) > 1 ~ dplyr::ntile(score, 10)/10,
        #   TRUE ~ 1)
        #   )
    })
  
  # write post-QC file
  df %>%
    purrr::reduce(dplyr::inner_join, by = c("InteractionID", "score")) %>%
    dplyr::select(chrom.x:end.x, chrom.y:end.y, score) %>%
    write.table(file = paste0(postQC_dir, ct, ".bedpe"),
                quote = F,
                row.names = F,
                col.names = F,
                sep = "\t")
  
  # save to HiChIP object
  HiChIP[[ct]] <- df
}

# save
saveRDS(HiChIP, file = "output/HiChIP/HiChIP.rds")

