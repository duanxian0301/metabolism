suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

project_root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
knk_root <- file.path(project_root, "results", "24_knk_nonlipid8_F1_PD")
summary_csv <- file.path(knk_root, "summary", "knk_summary_nonlipid8_F1_PD.csv")

kegg_rdata_candidates <- c(
  "D:/scPagwas/scPagwas-main/scPagwas-main/data/Genes_by_pathway_kegg.RData",
  "D:/scPagwas/data/Genes_by_pathway_kegg.RData"
)
kegg_rdata <- kegg_rdata_candidates[file.exists(kegg_rdata_candidates)][1]
if (is.na(kegg_rdata)) stop("Genes_by_pathway_kegg.RData not found")
load(kegg_rdata)

term2gene <- rbindlist(lapply(names(Genes_by_pathway_kegg), function(pw) {
  data.table(pathway_id = rep(pw, length(Genes_by_pathway_kegg[[pw]])), gene = as.character(Genes_by_pathway_kegg[[pw]]))
}))

group_to_dir <- c(
  Pericyte = "Pericyte",
  Oligodendrocyte_precursor_cell = "Oligodendrocyte_precursor_cell",
  VIP_GABAergic_cortical_interneuron = "VIP_GABAergic_cortical_interneuron"
)

rank_and_gsea <- function(diff_df) {
  if (!nrow(diff_df) || !("gene" %in% colnames(diff_df)) || !("Z" %in% colnames(diff_df))) {
    return(list(df = data.table(), top = NA_character_, n_sig = 0L))
  }
  gene_list <- setNames(as.numeric(diff_df[["Z"]]), as.character(diff_df[["gene"]]))
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  if (length(gene_list) < 50) {
    return(list(df = data.table(), top = NA_character_, n_sig = 0L))
  }
  res <- tryCatch(
    GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(df = data.table(), top = NA_character_, n_sig = 0L, err = conditionMessage(res)))
  }
  rdf <- as.data.table(as.data.frame(res))
  sig_df <- rdf[p.adjust < 0.05]
  top_paths <- if (nrow(sig_df)) paste(head(sig_df$ID, 5), collapse = "; ") else NA_character_
  list(df = rdf, top = top_paths, n_sig = nrow(sig_df))
}

sum_df <- fread(summary_csv)
for (i in seq_len(nrow(sum_df))) {
  row <- sum_df[i]
  prefix <- file.path(knk_root, group_to_dir[[row$cell_type]], paste(row$disease_dataset, row$cell_type, row$gene, sep = "_"), paste(row$disease_dataset, row$cell_type, row$gene, sep = "_"))
  pert_path <- paste0(prefix, "_perturbed_genes.csv")
  pw_path <- paste0(prefix, "_pathway_enrichment.csv")
  if (!file.exists(pert_path)) next
  diff_df <- fread(pert_path)
  gsea_info <- rank_and_gsea(diff_df)
  fwrite(gsea_info$df, pw_path)
  sum_df[i, top_pathways := gsea_info$top]
}
fwrite(sum_df, summary_csv)
