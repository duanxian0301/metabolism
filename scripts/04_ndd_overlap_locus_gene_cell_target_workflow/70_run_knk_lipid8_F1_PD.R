suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scTenifoldKnk)
  library(scTenifoldNet)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(data.table)
  library(ggplot2)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
shard_id <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
shard_n <- if (length(args) >= 2) as.integer(args[[2]]) else 1L
save_full_rds <- identical(Sys.getenv("KNK_SAVE_FULL_RDS", "0"), "1")

project_root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
input_root <- file.path(project_root, "work", "knk_inputs")
output_root <- file.path(project_root, "results", "25_knk_lipid8_F1_PD")
log_dir <- file.path(output_root, "logs")
sum_dir <- file.path(output_root, "summary")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

group_to_dir <- c(
  Pericyte = "Pericyte",
  Oligodendrocyte_precursor_cell = "Oligodendrocyte_precursor_cell"
)
for (sub in unique(unname(group_to_dir))) {
  dir.create(file.path(output_root, sub), recursive = TRUE, showWarnings = FALSE)
}

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

queue <- rbindlist(list(
  data.table(dataset = "MSSM_PD", group = "Pericyte", gene = c("TMEM175", "DGKQ", "LRRC37A2", "ARL17B"), priority = "tier1", rationale = "lipid8_F1-PD core genes tested in scPagwas-significant pericyte"),
  data.table(dataset = "MSSM_PD", group = "Oligodendrocyte_precursor_cell", gene = c("LRRC37A2", "ARL17B"), priority = "supplement", rationale = "lipid8_F1-PD OPC supplement genes")
))
if (shard_n > 1L) {
  keep_idx <- (((seq_len(nrow(queue)) - 1L) %% shard_n) + 1L) == shard_id
  queue <- queue[keep_idx, ]
}

log_msg <- function(logfile, ...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = logfile, append = TRUE)
}

read_group_to_seurat <- function(group_dir) {
  mtx <- Matrix::readMM(file.path(group_dir, "matrix.mtx"))
  barcodes <- fread(file.path(group_dir, "barcodes.tsv"), header = FALSE)$V1
  features <- fread(file.path(group_dir, "features.tsv"), header = FALSE)
  metadata <- fread(file.path(group_dir, "metadata.csv"), data.table = FALSE)
  rownames(metadata) <- metadata[[1]]
  metadata[[1]] <- NULL
  symbols <- features$V2
  symbols[is.na(symbols) | symbols == ""] <- features$V1[is.na(symbols) | symbols == ""]
  rownames(mtx) <- symbols
  colnames(mtx) <- barcodes
  dup <- duplicated(rownames(mtx))
  if (any(dup)) mtx <- mtx[!dup, ]
  keep <- Matrix::rowSums(mtx > 0) > (0.05 * ncol(mtx))
  mtx <- mtx[keep, ]
  metadata <- metadata[colnames(mtx), , drop = FALSE]
  CreateSeuratObject(counts = mtx, meta.data = metadata)
}

get_counts_mat <- function(obj) {
  tryCatch(
    GetAssayData(obj, layer = "counts"),
    error = function(e) GetAssayData(obj, slot = "counts")
  )
}

rank_and_gsea <- function(diff_df, out_prefix, logfile) {
  if (!nrow(diff_df) || !("gene" %in% colnames(diff_df)) || !("Z" %in% colnames(diff_df))) {
    fwrite(data.table(), paste0(out_prefix, "_pathway_enrichment.csv"))
    return(list(top = NA_character_, n_sig = 0))
  }
  gene_list <- setNames(as.numeric(diff_df[["Z"]]), as.character(diff_df[["gene"]]))
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  if (length(gene_list) < 50) {
    fwrite(data.table(), paste0(out_prefix, "_pathway_enrichment.csv"))
    return(list(top = NA_character_, n_sig = 0))
  }
  gsea_res <- tryCatch(
    GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(gsea_res, "error")) {
    log_msg(logfile, "GSEA failed: ", conditionMessage(gsea_res))
    fwrite(data.table(), paste0(out_prefix, "_pathway_enrichment.csv"))
    return(list(top = NA_character_, n_sig = 0))
  }
  gsea_df <- as.data.frame(gsea_res)
  fwrite(gsea_df, paste0(out_prefix, "_pathway_enrichment.csv"))
  sig_df <- gsea_df[gsea_df$p.adjust < 0.05, , drop = FALSE]
  top_paths <- if (nrow(sig_df)) paste(head(sig_df$ID, 5), collapse = "; ") else NA_character_
  list(top = top_paths, n_sig = nrow(sig_df))
}

plot_volcano <- function(diff_df, out_file) {
  if (!all(c("FC", "p.value", "p.adj", "Z") %in% colnames(diff_df))) return(invisible(NULL))
  diff_df$logFC <- log2(diff_df$FC + 1e-6)
  diff_df$negLogP <- -log10(diff_df$p.value + 1e-300)
  diff_df$significant <- diff_df$p.adj < 0.05 & abs(diff_df$Z) > 2
  p <- ggplot(diff_df, aes(x = logFC, y = negLogP)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("grey70", "#d62728")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    theme_classic(base_size = 12)
  ggsave(out_file, p, width = 6, height = 5, dpi = 200)
}

run_one <- function(dataset, group, gene, priority, rationale) {
  out_dir <- file.path(output_root, group_to_dir[[group]], paste(dataset, group, gene, sep = "_"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(log_dir, paste(dataset, group, gene, "log.txt", sep = "_"))
  prefix <- file.path(out_dir, paste(dataset, group, gene, sep = "_"))

  group_dir <- file.path(input_root, dataset, group)
  req_files <- file.path(group_dir, c("matrix.mtx", "barcodes.tsv", "features.tsv", "metadata.csv"))
  if (!all(file.exists(req_files))) {
    log_msg(logfile, "Skipping: missing prepared group input")
    return(data.table(disease_dataset = dataset, cell_type = group, gene = gene, input_cell_number = NA_integer_, ko_success = FALSE, perturbed_gene_n = NA_integer_, top_pathways = NA_character_, remarks = "Missing prepared group input"))
  }

  obj <- read_group_to_seurat(group_dir)
  counts0 <- get_counts_mat(obj)
  if (!(gene %in% rownames(counts0))) {
    log_msg(logfile, "Skipping: gene missing from expression matrix")
    return(data.table(disease_dataset = dataset, cell_type = group, gene = gene, input_cell_number = ncol(obj), ko_success = FALSE, perturbed_gene_n = NA_integer_, top_pathways = NA_character_, remarks = "Gene not in expression matrix"))
  }
  min_cells_required <- if (group == "Pericyte") 300L else 1000L
  if (ncol(obj) < min_cells_required) {
    log_msg(logfile, "Skipping: too few cells for stable KNK (need ", min_cells_required, ", have ", ncol(obj), ")")
    return(data.table(disease_dataset = dataset, cell_type = group, gene = gene, input_cell_number = ncol(obj), ko_success = FALSE, perturbed_gene_n = NA_integer_, top_pathways = NA_character_, remarks = paste("Too few cells:", ncol(obj))))
  }

  log_msg(logfile, "Running ", dataset, " / ", group, " / ", gene, " with ", ncol(obj), " cells")
  obj <- suppressWarnings(FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(6000, nrow(obj) - 1)))
  hvgs <- VariableFeatures(obj)
  counts <- get_counts_mat(obj)
  genes_use <- unique(c(gene, hvgs))
  genes_use <- intersect(genes_use, rownames(counts))
  data_knk <- as.matrix(counts[genes_use, ])
  if (ncol(data_knk) > 3000) {
    use_cells <- sample(colnames(data_knk), 3000)
    data_knk <- data_knk[, use_cells, drop = FALSE]
  }

  result <- tryCatch(
    scTenifoldKnk(countMatrix = data_knk, gKO = gene, qc = FALSE, nc_nNet = 10, nc_nCells = min(500, max(200, floor(ncol(data_knk) / 4))), nCores = 1),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    log_msg(logfile, "KO failed: ", conditionMessage(result))
    return(data.table(disease_dataset = dataset, cell_type = group, gene = gene, input_cell_number = ncol(data_knk), ko_success = FALSE, perturbed_gene_n = NA_integer_, top_pathways = NA_character_, remarks = paste("KO failed:", conditionMessage(result))))
  }

  diff_df <- as.data.frame(result$diffRegulation)
  fwrite(diff_df, paste0(prefix, "_perturbed_genes.csv"))
  if (save_full_rds) saveRDS(result, paste0(prefix, "_scTenifoldKnk.rds"))
  plot_volcano(diff_df, paste0(prefix, "_volcano.png"))
  gsea_info <- rank_and_gsea(diff_df, prefix, logfile)
  sig_n <- if ("p.adj" %in% colnames(diff_df)) sum(diff_df$p.adj < 0.05, na.rm = TRUE) else NA_integer_
  fwrite(data.table(dataset = dataset, group = group, gene = gene, priority = priority, rationale = rationale, input_cells = ncol(data_knk), gene_count = nrow(data_knk), perturbed_gene_n = sig_n, significant_pathway_n = gsea_info$n_sig), paste0(prefix, "_summary.csv"))

  data.table(disease_dataset = dataset, cell_type = group, gene = gene, input_cell_number = ncol(data_knk), ko_success = TRUE, perturbed_gene_n = sig_n, top_pathways = gsea_info$top, remarks = rationale)
}

make_convergence <- function(sum_df) {
  ok_df <- sum_df[ko_success == TRUE, ]
  if (!nrow(ok_df)) return(invisible(NULL))
  combos <- unique(ok_df[, .(disease_dataset, cell_type)])
  for (i in seq_len(nrow(combos))) {
    ds <- combos$disease_dataset[[i]]
    ct <- combos$cell_type[[i]]
    ss <- ok_df[disease_dataset == ds & cell_type == ct, ]
    if (nrow(ss) < 2) next
    gene_sets <- list()
    path_sets <- list()
    for (j in seq_len(nrow(ss))) {
      g <- ss$gene[[j]]
      prefix <- file.path(output_root, group_to_dir[[ct]], paste(ds, ct, g, sep = "_"), paste(ds, ct, g, sep = "_"))
      pert_path <- paste0(prefix, "_perturbed_genes.csv")
      pw_path <- paste0(prefix, "_pathway_enrichment.csv")
      if (file.exists(pert_path)) {
        df <- fread(pert_path)
        gene_sets[[g]] <- df$gene[df$p.adj < 0.05]
      }
      if (file.exists(pw_path)) {
        pdf <- fread(pw_path)
        path_sets[[g]] <- pdf$ID[pdf$p.adjust < 0.05]
      }
    }
    genes <- names(gene_sets)
    if (length(genes) >= 2) {
      mat <- matrix(NA_real_, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
      for (g1 in genes) for (g2 in genes) {
        u <- union(gene_sets[[g1]], gene_sets[[g2]])
        mat[g1, g2] <- if (length(u)) length(intersect(gene_sets[[g1]], gene_sets[[g2]])) / length(u) else NA_real_
      }
      fwrite(as.data.table(as.table(mat)), file.path(sum_dir, paste(ds, ct, "gene_overlap_jaccard.csv", sep = "_")))
    }
    pgenes <- names(path_sets)
    if (length(pgenes) >= 2) {
      mat <- matrix(NA_real_, nrow = length(pgenes), ncol = length(pgenes), dimnames = list(pgenes, pgenes))
      for (g1 in pgenes) for (g2 in pgenes) {
        u <- union(path_sets[[g1]], path_sets[[g2]])
        mat[g1, g2] <- if (length(u)) length(intersect(path_sets[[g1]], path_sets[[g2]])) / length(u) else NA_real_
      }
      fwrite(as.data.table(as.table(mat)), file.path(sum_dir, paste(ds, ct, "pathway_overlap_jaccard.csv", sep = "_")))
    }
  }
}

summary_rows <- rbindlist(lapply(seq_len(nrow(queue)), function(i) {
  row <- queue[i, ]
  run_one(row$dataset, row$group, row$gene, row$priority, row$rationale)
}), fill = TRUE)

suffix <- if (shard_n > 1L) paste0("_shard", shard_id, "of", shard_n) else ""
fwrite(summary_rows, file.path(sum_dir, paste0("knk_summary_lipid8_F1_PD", suffix, ".csv")))
make_convergence(summary_rows)
fwrite(queue, file.path(sum_dir, paste0("knk_execution_queue_lipid8_F1_PD", suffix, ".csv")))
