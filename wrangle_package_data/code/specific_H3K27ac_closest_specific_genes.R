setwd("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/") ; devtools::load_all("/working/lab_jonathb/alexandT/tgp/")

#inDir <- "/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/enhancer_activity_gene_target/correlate_enh_spec_expression/enh2gene_DHS_target/spec_signal/ENHbin8_9_10_EXPbin6_7_8_9_10/"
inDir <- "/working/lab_jonathb/jonathB/projects/jb_lab/target_gene_prediction/signal_matrix/build_matrix/output/enh2gene/enh2gene_DHS_target/spec_signal/ENHbin8_9_10_EXPbin6_7_8_9_10/"
# This dir contains bed files with cell type specific DHS sites (top decile) linked to nearest specifically expressed gene
# (top decile, ENSG because the gene expression data is at the level of gene rather than transcript).

specific_H3K27ac_closest_specific_genes <- list() ; for(file in list.files(inDir, pattern = "bed", full.names = T)){
  celltype <- file %>% basename %>% gsub("\\_.*", "", .) 
  # clean up cancer celltype names
  if(celltype == "BRST.MCF7"){celltype = "BRST.MCF7.CNCR"}
  if(celltype == "PRAD.LNCAP"){celltype = "PR.LNCAP.CNCR"}
  if(celltype == "GI.HCT116"){celltype = "GI.HCT116.CNCR"}
  print(celltype)
  specific_H3K27ac_closest_specific_genes[[celltype]] <- import_BED(file, metadata_cols = celltype)
}

specific_H3K27ac_closest_specific_genes <- specific_H3K27ac_closest_specific_genes %>%
  purrr::reduce(dplyr::full_join) %>%
  dplyr::select(chrom:end, everything())

# save
H3K27acDir <- "output/H3K27ac/" ; dir.create(H3K27acDir)
saveRDS(specific_H3K27ac_closest_specific_genes, file = paste0(H3K27acDir, "/specific_H3K27ac_closest_specific_genes.rds"))
