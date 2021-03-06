#!/usr/bin/Rscript

# load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))

# inputs
setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/")
dir.create("output/expression/", showWarnings = F)
dir.create("output/H3K27ac/", showWarnings = F)
metadata <- read.delim("output/metadata.tsv", header = T) %>% as_tibble

# functions
values_sum_fn <- function(x){sum(x, na.rm = T)}
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

##### IMPORT #####

# EXPRESSION #
# import TSSs
TSSs <- read.delim(gzfile("data/GENCODE/gencode.v34lift37.basic.tss.bed.gz"),
                   header = F) %>% as_tibble
names(TSSs) <- c("chrom", "start", "end", "ensg", "symbol", "enst")

# import expression data, combine and tidy
expression_files <- paste0("data/expression/",
                           filter(metadata, object == "expression")$accession,
                           ".tsv")
expression_raw <- expression_files %>%
  map(~ .x %>% read.delim %>% as_tibble %>% 
        rename(gene_id = 1, TPM = 2))
names(expression_raw) <- filter(metadata, object == "expression")$celltype
expression_raw <- expression_raw %>% 
  # bind rows
  bind_rows(.id = "celltype") %>%
  # get only ensgs, remove ensg extensions  
  filter(grepl("ENSG", gene_id)) %>%
  transmute(celltype, ensg = gene_id %>% gsub("\\.[0-9]*", "", .), TPM) %>%
  distinct %>%
  # spread each gene, sum duplicate gene entries
  pivot_wider(id_cols = "ensg", 
              names_from = celltype, values_from = TPM,
              values_fill = 0,
              values_fn = values_sum_fn) %>%
  # remove gene rows where all values are 0
  filter(rowSums(across(where(is.numeric))) > 0) %>%
  column_to_rownames("ensg") %>%
  # expression = log2(TPM + 1) 
  mutate(across(everything(), ~ log2(.x + 1))) %>%
  as.matrix

# H3K27ac #
# import H3K27ac matrix, replace codes with names
H3K27ac_raw <- read.delim("output/H3K27ac/H3K27ac_input_matrix.tsv") %>%
  mutate(region = paste0(X..chr., ":", X.start., "-", X.end.)) %>%
  select(region, everything(), -X..chr., -X.start., -X.end.) %>% 
  column_to_rownames("region") %>% as.matrix
acc_to_ct <- tibble(accession = colnames(H3K27ac_raw) %>% gsub("X\\.", "", .) %>% gsub("\\..*", "", .)) %>%
  full_join(metadata %>% filter(object == "H3K27ac")) 
if(any(is.na(acc_to_ct$celltype))){stop(filter(acc_to_ct, is.na(celltype))$accession, " does not have a metadata entry!")}
colnames(H3K27ac_raw) <- acc_to_ct$celltype

##### POST-PROCESSING #####
# create list
process <- list()
process$raw <- list(H3K27ac = H3K27ac_raw, expression = expression_raw) 
dims <- process$raw %>% map(dimnames) 

# signal = quantile normalise
process$qn <- process$raw %>% names %>%
  sapply(function(x){
    x_raw <- process$raw[[x]]
    # ppC::normalize.quantiles() has improper NA handling, convert to 0
    x_raw[is.na(x_raw)] <- 0   
    # get signal (qn)
    x_signal <- preprocessCore::normalize.quantiles(x_raw, copy = T) 
    dimnames(x_signal) <- dims[[x]]
    # get specificity
    x_specificity <- x_signal %>% apply(1, specificity_score) %>% t 
    dimnames(x_specificity) <- dims[[x]]
    # return list
    list(signal = x_signal, specificity = x_specificity) %>% 
      map(~.x %>% {
        if(x == "H3K27ac")          dplyr::as_tibble(., rownames = "DHS") 
        else if(x == "expression")  dplyr::as_tibble(., rownames = "ensg")
      })
}, USE.NAMES = T, simplify = F)

# bin
process$binned <- process$qn %>% names %>%
  sapply(function(x){
    process$qn[[x]] %>% names %>%
      sapply(function(y){
        process$qn[[x]][[y]] %>% mutate(across(where(is.numeric), ~ ntile(.x, 10)/10)) 
      }, USE.NAMES = T, simplify = F)
}, USE.NAMES = T, simplify = F)

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
H3K27ac_specificity_ranked <- process$qn$H3K27ac$specificity %>%
  mutate(across(where(is.numeric), rank))

##### SAVE #####
# intermediate files
process %>% names %>%
  sapply(function(x){ print(x)
    process[[x]] %>% names %>%
      sapply(function(y){ print(y)
        if(!is.list(process[[x]][[y]])){
          write.table(process[[x]][[y]],
                      paste0("output/", y, "/", y, "_", x, ".tsv"), sep = "\t", quote = F)
        } else {
          process[[x]][[y]] %>% names %>%
            sapply(function(z){ print(z)
                write.table(process[[x]][[y]][[z]],
                            paste0("output/", y, "/", y, "_", z, "_", x, ".tsv"), sep = "\t", quote = F)
            })
        }
      })
  })
# objects
H3K27ac <- process$binned$H3K27ac
expression <- process$binned$expression
saveRDS(expression,                 "output/expression.rds")
saveRDS(expressed,                  "output/expressed.rds")
saveRDS(H3K27ac,                    "output/H3K27ac.rds")
saveRDS(H3K27ac_specificity_ranked, "output/H3K27ac_specificity_rank.rds")

