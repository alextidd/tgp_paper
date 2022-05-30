# deleterious coding mutations - only those found in GENCODE Basic protein-coding # TODO: fix this in upstream scripts
wkdir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/sysdata/"
setwd(wkdir) 
devtools::load_all("/working/lab_jonathb/alexandT/EG2")
coding_mutations_dir <- paste0(wkdir, "output/coding_mutations/")

c("missense", "nonsense", "splicesite") %>%
  sapply(function(mutation_type){
    print(mutation_type)
    df <- paste0(coding_mutations_dir, mutation_type, "_SNVs.tsv") %>%
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
      write_tibble(paste0(coding_mutations_dir, mutation_type, ".tsv"))
  })


