# Convert FitHiChIP to bedpe
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

option.list <- list(
  make_option("--in.FitHiChIP", type="character"),
  make_option("--out.bedpe", type="character"))
opt <- parse_args(OptionParser(option_list=option.list))

x <- read.delim(opt$in.FitHiChIP)
x2 <- x %>% mutate(score = -log10(Q.Value_Bias)) %>% select(chr1:e2, score)
x2 %>% write.table(opt$out.bedpe, quote = F, row.names = F, col.names = F, sep = "\t")
