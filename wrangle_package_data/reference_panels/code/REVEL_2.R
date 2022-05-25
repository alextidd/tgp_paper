# REVEL - deleterious coding variants - only those found in GENCODE Basic protein-coding # TODO: fix this in upstream scripts
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/"
setwd(wkdir) 
devtools::load_all("/working/lab_jonathb/alexandT/tgp")
REVELDir <- paste0(wkdir, "output/REVEL/")

c("missense", "nonsense", "splicesite") %>%
  sapply(function(mutation_type){
    print(mutation_type)
    df <- paste0(REVELDir, mutation_type, "_SNVs.tsv") %>%
      read_tibble(header = T) 
    
    df %>%
      tidyr::separate_rows(ensgs) %>%
      dplyr::group_by(chrom, position) %>%
      {if("score" %in% colnames(df)) 
        # get an average REVEL score per position
        dplyr::summarise(.,
                         score = mean(score),
                         ensgs = ensgs %>% unique %>% paste(collapse = ";")) 
        else 
          dplyr::summarise(.,
                           ensgs = ensgs %>% unique %>% paste(collapse = ";"))
      } %>%
      write_tibble(paste0(REVELDir, mutation_type, ".tsv"))
  })


# wrangle_REVEL <- function(df){
#   df %>% 
#     tidyr::separate_rows(ensgs) %>%
#     dplyr::group_by(chrom, position) %>%
#     {if("score" %in% colnames(df)) 
#       # get an average REVEL score per position
#       dplyr::summarise(.,
#                        score = mean(score),
#                        ensgs = ensgs %>% unique %>% paste(collapse = ";")) 
#       else 
#         dplyr::summarise(.,
#                          ensgs = ensgs %>% unique %>% paste(collapse = ";"))
#     }
# }
# # fix and rewrite
# paste0(REVELDir, "missense_SNVs.tsv") %>% 
#   read_tibble %>% 
#   dplyr::select(chrom = V1, position = V2, ensgs = V3) %>%
#   wrangle_REVEL %>%
#   write_tibble(paste0(REVELDir, "missense_SNVs.tsv"))
# paste0(REVELDir, "nonsense_SNVs.tsv") %>% 
#   read_tibble %>% 
#   dplyr::select(chrom = V1, position = V2, ensgs = V3) %>%
#   wrangle_REVEL %>%
#   write_tibble(paste0(REVELDir, "nonsense_SNVs.tsv"))
# paste0(REVELDir, "splicesite_SNVs.tsv") %>% 
#   read_tibble %>%
#   dplyr::select(chrom = V1, position = V2, ensgs = V3) %>%
#   wrangle_REVEL %>%
#   write_tibble(paste0(REVELDir, "splicesite_SNVs.tsv"))
