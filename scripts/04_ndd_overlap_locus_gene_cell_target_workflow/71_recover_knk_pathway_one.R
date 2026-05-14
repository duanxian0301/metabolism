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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 71_recover_knk_pathway_one.R <perturbed_genes.csv> <pathway_enrichment.csv>")
}

pert_path <- args[[1]]
pw_path <- args[[2]]

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

diff_df <- fread(pert_path)
gene_list <- setNames(as.numeric(diff_df[["Z"]]), as.character(diff_df[["gene"]]))
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

if (length(gene_list) < 50) {
  fwrite(data.table(), pw_path)
  quit(save = "no", status = 0)
}

res <- GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

fwrite(as.data.table(as.data.frame(res)), pw_path)
