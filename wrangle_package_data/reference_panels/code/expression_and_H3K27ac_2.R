#!/usr/bin/Rscript
setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/")
dir.create("output/expression/")
dir.create("output/H3K27ac/")

# load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))

# functions
specificity_score <- function(x){
  norm_vec <- function(x) sqrt(sum(x^2))
  # function to calculate cell type specificity scores for a vector of values x
  # Inputs: x is vector of scores for cell types
  # Output: A vector of length cell types, giving the specificity score for each cell type
  
  n_cells <- length(x)
  spec_score <- vector(length=n_cells)
  
  # compute specificity score for each cell type
  for (i in 1:n_cells) {
    
    if(is.na(x[i])) {
      spec_score[i] <- NA
    } else {
      spec_score[i] <- x[i]/norm_vec(x)
    }
  }
  return(spec_score)
}
bin <- function(x, x_dims){
  out <- x %>% as_tibble %>% mutate(across(where(is.numeric), ~ ntile(.x, 10)/10)) %>% as.matrix
  dimnames(out) <- x_dims
  return(out)
  }

##### EXPRESSION #####
# import TSSs
TSSs <- read.delim(gzfile("data/GENCODE/gencode.v34lift37.basic.tss.bed.gz"),
                   header = F) %>% as_tibble
names(TSSs) <- c("chrom", "start", "end", "ensg", "symbol", "enst")

# import expression data, combine and tidy
expression_raw <- dir("data/expression/", full.names = T) %>%
  map(~ .x %>% read.delim %>% as_tibble %>% 
        rename(gene_id = 1, TPM = 2) )
names(expression_raw) <- dir("data/expression/") %>% tools::file_path_sans_ext()
expression_raw <- expression_raw %>% 
  # bind rows
  bind_rows(.id = "acc") %>%
  # get only ensgs, remove ensg extensions  
  filter(grepl("ENSG", gene_id)) %>%
  transmute(acc, ensg = gene_id %>% gsub("\\.[0-9]*", "", .), TPM) %>%
  distinct %>%
  # spread each gene, sum duplicate gene entries
  spread(acc, TPM) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  group_by(ensg) %>%
  summarise(across(everything(), sum)) %>%
  # remove gene rows where all values are 0
  filter(rowSums(across(where(is.numeric))) > 0) %>%
  column_to_rownames("ensg") %>%
  # expression = log2(TPM + 1) 
  mutate(across(everything(), ~ log2(.x + 1))) %>%
  as.matrix

##### H3K27ac #####
# import H3K27ac matrix
H3K27ac_raw <- read.delim("output/H3K27ac/H3K27ac_input_matrix.tsv") %>%
  mutate(region = paste0(X..chr.,":",X.start.,"-",X.end.)) %>%
  select(region, everything(), -X..chr., -X.start., -X.end.) %>% 
  column_to_rownames("region") %>% as.matrix
colnames(H3K27ac_raw) <- colnames(H3K27ac_raw) %>%
  gsub("X\\.", "", .) %>% gsub("\\..*", "", .)

##### POST-PROCESSING #####

# create list, replace codes with names
metadata <- read.delim("output/metadata.tsv", header = T) %>% as_tibble
raw <- list(H3K27ac = H3K27ac_raw, expression = expression_raw)
dims <- list() ; for(x in c("H3K27ac", "expression")){
  dims[[x]] <- dimnames(raw[[x]])
  dims[[x]][[2]] <- metadata %>%
      filter(object == x) %>%
      {dplyr::left_join(dplyr::tibble(colname = colnames(raw[[x]])), ., by = "colname")} %>%
      dplyr::pull(celltype)
}

# quantile normalise
qn <- list() ; c("H3K27ac", "expression") %>%
  sapply(function(x){
    # qn
    x_qn <- raw[[x]] %>% preprocessCore::normalize.quantiles(copy = T) 
    x_qn[is.na(x_qn)] <- 0    
    dimnames(x_qn) <- dims[[x]]
    return(x_qn)
  }, simplify = F, USE.NAMES = T) -> qn

# bin
binned <- list() ; qn %>% names %>%
  sapply(function(x){
    x_qn <- qn[[x]]
    # calculate signal scores, bin
    signal <- x_qn %>% bin(dims[[x]])
    # calculate specificity scores, bin
    specificity <- x_qn %>% apply(1, specificity_score) %>% t %>% bin(dims[[x]])
    # return list
    list(specificity = specificity, signal = signal) %>% 
      map(~.x %>% {
          if(x == "H3K27ac")          dplyr::as_tibble(., rownames = "DHS") 
          else if(x == "expression")  dplyr::as_tibble(., rownames = "ensg")
      })
  }, simplify = F, USE.NAMES = T
  ) -> binned
expression <- binned$expression 
H3K27ac <- binned$H3K27ac 

# expressed = expression binary
expressed <- expression_raw 
expressed[expressed!=0] <- 1 
dimnames(expressed) <- dims$expression
expressed <- expressed %>%
  as_tibble(rownames = "ensg")

##### H3K27ac SPECIFICITY RANKINGS FOR CELLTYPE ENRICHMENT ####
# Re-implementation of method in B. Soskic et al., Chromatin activity at GWAS loci identifies T cell states driving complex immune diseases. Nat. Genet. (2019).
# For each regulatory site, quantify H3K27ac per region for each cell type in the reference panel.
# Compute specificity of each site across reference panel of cell types, and rank sites by specicicity.
# Intersect SNPs with DHS sites and compute the mean rank of the intersected DHS.
# Assuming SNPs overlap sites randomly, compute significance of mean rank of sites where SNPs overlap based on deivation from uniform distribution.
H3K27ac_specificity <- qn$H3K27ac %>% apply(1, specificity_score) %>% t 
dimnames(H3K27ac_specificity) <- dims$H3K27ac
H3K27ac_specificity_ranked <- H3K27ac_specificity %>%
  as_tibble(rownames = "DHS") %>%
  mutate(across(where(is.numeric), rank))

##### SAVE ####
saveRDS(expression, "output/expression/expression.rds")
saveRDS(expressed, "output/expression/expressed.rds")
saveRDS(H3K27ac, "output/H3K27ac/H3K27ac.rds")
saveRDS(H3K27ac_specificity_ranked, "output/H3K27ac/H3K27ac_specificity_rank.rds")
