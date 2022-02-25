#!/usr/bin/Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Hs.eg.db))

file <- commandArgs(trailingOnly = T)[1]      # file name
keep_ensg <- commandArgs(trailingOnly = T)[2] %>% as.logical() # keep ensg column?

# load in genes file
data <- read.delim(gzfile(file),
                   header = T,
                   stringsAsFactors = F)
cols_order <- colnames(data) 
new_cols_order <- cols_order %>% gsub("^ensg$", "symbol", .)

ensg <- data$ensg
symbol <- mapIds(org.Hs.eg.db, 
                  keys = ensg, 
                  keytype = "ENSEMBL",
                  column="SYMBOL") 

ensg_to_symbol <- tibble(ensg = names(symbol),
                         symbol = symbol) %>%
  distinct

lost <- ensg_to_symbol %>% filter(is.na(symbol)) %>% n_distinct
converted <- ensg_to_symbol %>% filter(!is.na(symbol)) %>% n_distinct
total <- ensg_to_symbol %>% dplyr::select(ensg) %>% n_distinct

cat(converted, "/", total,
    "gene ENSGs successfully converted to gene symbols.",
    lost, 
    "lost.\n")

# add symbols
data <- data %>% 
  left_join(ensg_to_symbol, by = "ensg")
if(!keep_ensg){
  # reorder to replace ensg column position with symbol column
  data <- data[,new_cols_order]
}

# save
outfile <- file %>% 
  gsub("\\.tsv.*", "", .) %>% 
  paste0("_with_symbols.tsv")
write.table(data,
            outfile,
            sep = "\t",
            row.names = F, col.names = T,
            quote = F)
cat("Output saved to", outfile,"\n")
