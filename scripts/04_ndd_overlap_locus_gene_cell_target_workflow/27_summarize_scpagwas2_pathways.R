options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 27_summarize_scpagwas2_pathways.R <run_dir> <prefix>")
}

run_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
prefix <- args[[2]]

score_path <- file.path(run_dir, paste0(prefix, "_singlecell_scPagwas_score_pvalue.Result.csv"))
sig_celltype_path <- file.path(run_dir, paste0(prefix, "_significant_celltypes.csv"))
pathway_matrix_path <- file.path(run_dir, paste0(prefix, "_Pathway_singlecell_lm_results.txt"))
pagwas_rds_path <- file.path(run_dir, paste0("Pagwas_data_", prefix, ".rds"))
if (!file.exists(pagwas_rds_path)) {
  alt_rds_path <- file.path("D:/scPagwas/metabolic_scpagwas2", prefix, paste0("Pagwas_data_", prefix, ".rds"))
  if (file.exists(alt_rds_path)) {
    pagwas_rds_path <- alt_rds_path
  }
}

if (!file.exists(score_path)) stop("Missing score file: ", score_path)
if (!file.exists(sig_celltype_path)) stop("Missing significant celltypes file: ", sig_celltype_path)
if (!file.exists(pathway_matrix_path)) stop("Missing pathway matrix file: ", pathway_matrix_path)
if (!file.exists(pagwas_rds_path)) stop("Missing pagwas RDS: ", pagwas_rds_path)

score_df <- fread(score_path, data.table = FALSE)
cell_col <- names(score_df)[1]
rownames(score_df) <- score_df[[cell_col]]
score_df[[cell_col]] <- NULL

sig_cells <- fread(sig_celltype_path, data.table = FALSE)
if (!nrow(sig_cells)) {
  stop("No significant cell types found in: ", sig_celltype_path)
}

genes_by_pathway <- readRDS("D:/scPagwas/Genes_by_pathway_kegg.rds")
pathway_gene_n <- vapply(genes_by_pathway, length, integer(1))

pathway_mat <- fread(pathway_matrix_path, data.table = FALSE)
first_col <- names(pathway_mat)[1]
if (!(first_col %in% names(pathway_gene_n))) {
  rownames(pathway_mat) <- pathway_mat[[first_col]]
  pathway_mat[[first_col]] <- NULL
}
if (nrow(pathway_mat) != nrow(score_df)) {
  stop("Pathway matrix rows (", nrow(pathway_mat), ") do not match score rows (", nrow(score_df), ")")
}
if (is.null(rownames(pathway_mat)) || !all(rownames(pathway_mat) == rownames(score_df))) {
  rownames(pathway_mat) <- rownames(score_df)
}

pagwas_obj <- readRDS(pagwas_rds_path)
meta_df <- pagwas_obj@meta.data
if (!"cell_type" %in% colnames(meta_df)) {
  stop("Missing cell_type in Pagwas meta.data")
}

meta_df <- meta_df[rownames(score_df), , drop = FALSE]
if (!all(rownames(meta_df) == rownames(score_df))) {
  stop("Meta data cell order does not match score file")
}

trs <- score_df[["scPagwas.TRS.Score"]]
if (is.null(trs)) stop("Missing scPagwas.TRS.Score in score file")

ref_dir <- "D:/文章/GS/scpagwas/ALPS_F2_MSSM_AD/Pathway_TRS"
pathway_name_map <- NULL
if (dir.exists(ref_dir)) {
  ref_files <- list.files(ref_dir, pattern = "^Result_.*_Pathway_vs_TRS_all[.]csv$", full.names = TRUE)
  if (length(ref_files) > 0) {
    ref_map <- rbindlist(lapply(ref_files, function(p) {
      x <- tryCatch(fread(p), error = function(e) NULL)
      if (is.null(x) || !all(c("pathway_id", "pathway_name") %in% names(x))) return(NULL)
      unique(x[, .(pathway_id, pathway_name)])
    }), fill = TRUE)
    if (!is.null(ref_map) && nrow(ref_map)) {
      ref_map <- unique(ref_map[!is.na(pathway_id) & !is.na(pathway_name), ])
      pathway_name_map <- setNames(ref_map$pathway_name, ref_map$pathway_id)
    }
  }
}

sanitize_celltype <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", trimws(x))
}

cor_pvalue <- function(r, n) {
  r <- pmin(pmax(r, -0.999999999), 0.999999999)
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t_stat), df = n - 2)
}

pathway_out_dir <- file.path(run_dir, "Pathway_TRS")
dir.create(pathway_out_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()

for (i in seq_len(nrow(sig_cells))) {
  celltype_label <- sig_cells$celltype[[i]]
  keep_idx <- which(meta_df$cell_type == celltype_label & !is.na(trs))
  if (length(keep_idx) < 20) next

  cell_ids <- rownames(meta_df)[keep_idx]
  trs_vec <- trs[keep_idx]
  mat_sub <- pathway_mat[cell_ids, , drop = FALSE]

  r_vals <- vapply(mat_sub, function(col) suppressWarnings(cor(col, trs_vec, method = "pearson", use = "pairwise.complete.obs")), numeric(1))
  ok <- !is.na(r_vals)
  if (!any(ok)) next

  pathway_ids <- names(r_vals)[ok]
  r_vals <- r_vals[ok]
  n_cells <- length(trs_vec)
  p_vals <- cor_pvalue(r_vals, n_cells)
  fdr_vals <- p.adjust(p_vals, method = "BH")
  gene_n <- pathway_gene_n[pathway_ids]
  pathway_names <- if (!is.null(pathway_name_map)) pathway_name_map[pathway_ids] else pathway_ids

  out_df <- data.frame(
    pathway_id = pathway_ids,
    cell_type = celltype_label,
    n_genes = as.integer(gene_n),
    n_cells = n_cells,
    r = as.numeric(r_vals),
    p_value = as.numeric(p_vals),
    FDR = as.numeric(fdr_vals),
    pathway_name = pathway_names,
    abs_r = abs(as.numeric(r_vals)),
    stringsAsFactors = FALSE
  )
  out_df <- out_df[order(out_df$FDR, out_df$p_value, -out_df$abs_r), ]

  safe_cell <- sanitize_celltype(celltype_label)
  fwrite(out_df, file.path(pathway_out_dir, paste0("Result_", safe_cell, "_Pathway_vs_TRS_all.csv")))
  fwrite(out_df[out_df$FDR < 0.05, ], file.path(pathway_out_dir, paste0("Result_", safe_cell, "_Pathway_vs_TRS_FDRsig.csv")))

  sig_df <- out_df[out_df$FDR < 0.05, , drop = FALSE]
  if (nrow(sig_df)) {
    sig_df$analysis <- prefix
    summary_rows[[length(summary_rows) + 1L]] <- sig_df
  }
}

summary_path <- file.path(run_dir, paste0(prefix, "_pathway_summary_sig.csv"))
if (length(summary_rows)) {
  summary_df <- rbindlist(summary_rows, fill = TRUE)
  fwrite(summary_df, summary_path)
} else {
  fwrite(data.frame(
    pathway_id = character(),
    cell_type = character(),
    n_genes = integer(),
    n_cells = integer(),
    r = numeric(),
    p_value = numeric(),
    FDR = numeric(),
    pathway_name = character(),
    abs_r = numeric(),
    analysis = character(),
    stringsAsFactors = FALSE
  ), summary_path)
}

message("Wrote pathway summaries to: ", pathway_out_dir)
message("Wrote summary table to: ", summary_path)
