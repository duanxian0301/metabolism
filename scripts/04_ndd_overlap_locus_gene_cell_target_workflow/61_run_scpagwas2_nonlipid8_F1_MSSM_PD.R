options(stringsAsFactors = FALSE)
options(timeout = 600)
Sys.setenv(OPENBLAS_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(scPagwas)
  library(Seurat)
  library(data.table)
})

# Compatibility patch: current scPagwas still calls GetAssayData() with slot=,
# which is defunct under current SeuratObject. We keep the F2-AD workflow but
# patch this accessor to use LayerData() directly on the assay object.
patched_Single_data_input <- function(Pagwas, Single_data, nfeatures = NULL, Pathway_list = NULL,
                                      assay = "RNA", min_clustercells = 5) {
  message("Input single cell data!")
  if (!("Seurat" %in% class(Single_data))) {
    message("Single_data is not of class Seurat!")
    return(Pagwas)
  }
  Celltype_anno <- data.frame(
    cellnames = rownames(Single_data@meta.data),
    annotation = as.vector(SeuratObject::Idents(Single_data))
  )
  rownames(Celltype_anno) <- Celltype_anno$cellnames
  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(length(Afterre_cell_types), " cell types are remain, after filter!")
  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in% Afterre_cell_types, ]
  Single_data <- Single_data[, Celltype_anno$cellnames]
  Pagwas$Celltype_anno <- Celltype_anno
  Pagwas$data_mat <- SeuratObject::LayerData(Single_data[[assay]], layer = "data")
  merge_scexpr <- Seurat::AggregateExpression(Single_data, assays = assay)[[assay]]
  Pagwas$VariableFeatures <- rownames(Pagwas$data_mat)
  rm(Single_data)
  pagene <- intersect(unique(unlist(Pathway_list)), rownames(Pagwas$data_mat))
  if (length(pagene) < 100) {
    stop("There are little match between rownames of Single_data and pathway genes!")
  }
  Pagwas$VariableFeatures <- intersect(Pagwas$VariableFeatures, pagene)
  if (ncol(merge_scexpr) == 1) {
    a <- colnames(merge_scexpr)
    Pagwas$merge_scexpr <- data.matrix(merge_scexpr[Pagwas$VariableFeatures, ])
    rownames(Pagwas$merge_scexpr) <- Pagwas$VariableFeatures
    colnames(Pagwas$merge_scexpr) <- a
  } else {
    Pagwas$merge_scexpr <- merge_scexpr[Pagwas$VariableFeatures, ]
  }
  return(Pagwas)
}
assignInNamespace("Single_data_input", patched_Single_data_input, ns = "scPagwas")

patched_Get_CorrectBg_p <- function(Single_data, scPagwas.TRS.Score, iters_singlecell,
                                    n_topgenes, scPagwas_topgenes, assay = "RNA") {
  gene_matrix <- SeuratObject::LayerData(Single_data[[assay]], layer = "data")
  mat_ctrl_raw_score <- matrix(0, nrow = ncol(gene_matrix), ncol = iters_singlecell)
  dic_ctrl_list <- list()
  pb <- txtProgressBar(style = 3)
  for (i in 1:iters_singlecell) {
    set.seed(i)
    dic_ctrl_list[[i]] <- sample(rownames(Single_data), n_topgenes)
    Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(dic_ctrl_list[[i]]), name = c("contr_genes"))
    mat_ctrl_raw_score[, i] <- Single_data$contr_genes1
    Single_data$contr_genes1 <- NULL
    setTxtProgressBar(pb, i/iters_singlecell)
  }
  close(pb)
  genes <- intersect(rownames(Single_data), rownames(gene_matrix))
  scPagwas_topgenes <- intersect(scPagwas_topgenes, genes)
  gene_matrix <- gene_matrix[genes, ]
  df_gene <- data.frame(gene = genes, var = apply(gene_matrix, 1, var))
  rownames(df_gene) <- df_gene$gene
  v_var_ratio_c2t <- rep(1, iters_singlecell)
  for (i_ctrl in 1:iters_singlecell) {
    v_var_ratio_c2t[i_ctrl] <- sum(df_gene[dic_ctrl_list[[i_ctrl]], "var"])
  }
  v_var_ratio_c2t <- v_var_ratio_c2t / sum(df_gene[scPagwas_topgenes, "var"])
  correct_pdf <- scPagwas::correct_background(scPagwas.TRS.Score, mat_ctrl_raw_score, v_var_ratio_c2t)
  rownames(correct_pdf) <- colnames(Single_data)
  return(correct_pdf)
}
assignInNamespace("Get_CorrectBg_p", patched_Get_CorrectBg_p, ns = "scPagwas")

root_out <- "D:/scPagwas/metabolic_scpagwas2"
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)
setwd(root_out)

resource_dir <- "D:/scPagwas"
genes_by_pathway_kegg <- readRDS(file.path(resource_dir, "Genes_by_pathway_kegg.rds"))
block_annotation_hg37 <- readRDS(file.path(resource_dir, "block_annotation_hg37.rds"))
chrom_ld <- readRDS(file.path(resource_dir, "chrom_ld.rds"))

gwas_file <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/clean_factor_inputs/nonlipid8_F1_clean.txt"
sc_file <- "D:/scPagwas/MSSM_PD_20k_from_ALPS_PD.rds"
prefix <- "nonlipid8_F1_MSSM_PD"

read_gwas <- function(path) {
  x <- fread(path)
  setnames(
    x,
    old = c("SNP", "CHR", "BP", "BETA", "SE", "FREQ"),
    new = c("rsid", "chrom", "pos", "beta", "se", "maf"),
    skip_absent = FALSE
  )
  x <- x[
    !is.na(rsid) & !is.na(chrom) & !is.na(pos) &
      !is.na(beta) & !is.na(se) & se > 0 &
      !is.na(maf) & maf > 0 & maf < 1 &
      !is.na(P) & P > 0 & P <= 1 &
      !is.na(N) & N > 0
  ]
  x[, `:=`(
    chrom = as.character(chrom),
    pos = as.integer(pos),
    maf = as.numeric(maf),
    beta = as.numeric(beta),
    se = as.numeric(se)
  )]
  as.data.frame(x)
}

gwas <- read_gwas(gwas_file)
sce <- readRDS(sc_file)
if (!"cell_type" %in% colnames(sce@meta.data)) {
  stop("Missing cell_type metadata in MSSM_PD_20k.rds")
}
Idents(sce) <- "cell_type"

run_dir <- file.path(root_out, prefix)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

pagwas_data <- scPagwas_main2(
  Pagwas = NULL,
  gwas_data = gwas,
  Single_data = sce,
  output.prefix = prefix,
  output.dirs = prefix,
  block_annotation = block_annotation_hg37,
  assay = "RNA",
  Pathway_list = genes_by_pathway_kegg,
  n.cores = 1,
  iters_singlecell = 100,
  chrom_ld = chrom_ld,
  singlecell = TRUE,
  celltype = TRUE
)

saveRDS(pagwas_data, file = file.path(run_dir, paste0("Pagwas_data_", prefix, ".rds")))

merged_path <- file.path(run_dir, paste0(prefix, "_Merged_celltype_pvalue.csv"))
if (file.exists(merged_path)) {
  merged <- fread(merged_path)
  setnames(merged, old = names(merged)[1], new = "row_id")
  merged[, celltype_FDR := p.adjust(pvalue, method = "BH")]
  setorder(merged, celltype_FDR, pvalue)
  fwrite(merged, file.path(run_dir, paste0(prefix, "_Merged_celltype_pvalue_withFDR.csv")))
  fwrite(merged[celltype_FDR < 0.05 & !is.na(celltype)], file.path(run_dir, paste0(prefix, "_significant_celltypes.csv")))
}

message("scPagwas2 finished: ", run_dir)
