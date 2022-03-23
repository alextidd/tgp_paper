setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/") 
inDir <- "data/expression_and_H3K27ac/output/"
outDir <- "output/"
packageDir <- "/working/lab_jonathb/alexandT/tgp/" ; devtools::load_all(packageDir)

# load 
all_metadata <- read_tibble(paste0(packageDir, "reference_data/data/all_metadata.tsv"), header = T)

# read matrices, fix names #
mats <- list() ; for(object in c("expression", "H3K27ac")){
  for(method in c("signal", "spec")){
    # get filepath
    file <- paste0(inDir, object, "/", object, "_qn_matrix_", method, ".rds")
    method_name <- method ; if(method == "spec"){ method_name <- "specificity" } ; cat(object, method_name, "\n")
    # read mat
    mat <- readRDS(file)
    # fix colnames
    colnames(mat) <- all_metadata %>%
      dplyr::filter(object == object) %>%
      {dplyr::left_join(dplyr::tibble(colname = colnames(mat)), ., by = "colname")} %>%
      dplyr::pull(celltype)
    # convert coords to bed (H3K27ac)
    mat <- mat %>% 
      { if(object == "H3K27ac") 
          dplyr::as_tibble(., rownames = "DHS") %>%
          tidyr::separate(DHS, into = c("chrom", "start", "end"), remove = F, convert = T) %>%
          dplyr::select(chrom, start, end, DHS, everything())
        else if(object == "expression") 
          dplyr::as_tibble(., rownames = "ensg")
      } 
    # write to list
    mats[[object]][[method_name]] <- mat
  }
} ; expression <- mats$expression ; H3K27ac <- mats$H3K27ac

# save objects
saveRDS(expression, file = "output/expression/expression.rds")
saveRDS(H3K27ac, file = "output/H3K27ac/H3K27ac.rds")

# expressed object
expressed <- readRDS(paste0(inDir, "expression/expression_raw_matrix.rds")) 
# replace codes with names
expr_metadata <- all_metadata %>% dplyr::filter(object == "expression")
colnames(expressed) <- expr_metadata[match(colnames(expressed),
                                          expr_metadata$file),]$name
# convert to binary
expressed[expressed!=0]<-1
# to tibble
expressed <- expressed %>%
  dplyr::as_tibble(rownames = "ensg")
# save
saveRDS(expressed, "output/expression/expressed.rds")
