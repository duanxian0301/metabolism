suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

outdir <- "D:/codex/GenomicSEM/metabolic/figures/figure4_scpagwas_umap_data"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

configs <- list(
  list(
    key = "HDL_core_AD",
    label = "HDL-core / AD",
    rds = "D:/scPagwas/metabolic_scpagwas2/lipid8_F2_MSSM_AD/Pagwas_data_lipid8_F2_MSSM_AD.rds",
    score = "D:/scPagwas/metabolic_scpagwas2/lipid8_F2_MSSM_AD/lipid8_F2_MSSM_AD_singlecell_scPagwas_score_pvalue.Result.csv"
  ),
  list(
    key = "Ketone_core_PD",
    label = "Ketone-core / PD",
    rds = "D:/scPagwas/metabolic_scpagwas2/nonlipid8_F1_MSSM_PD/Pagwas_data_nonlipid8_F1_MSSM_PD.rds",
    score = "D:/scPagwas/metabolic_scpagwas2/nonlipid8_F1_MSSM_PD/nonlipid8_F1_MSSM_PD_singlecell_scPagwas_score_pvalue.Result.csv"
  ),
  list(
    key = "TG_VLDL_core_PD",
    label = "TG/VLDL-core / PD",
    rds = "D:/scPagwas/metabolic_scpagwas2/lipid8_F1_MSSM_PD/Pagwas_data_lipid8_F1_MSSM_PD.rds",
    score = "D:/scPagwas/metabolic_scpagwas2/lipid8_F1_MSSM_PD/lipid8_F1_MSSM_PD_singlecell_scPagwas_score_pvalue.Result.csv"
  )
)

short_celltype <- function(x) {
  map <- c(
    "oligodendrocyte" = "Oligodendrocyte",
    "microglial cell" = "Microglia",
    "lamp5 GABAergic cortical interneuron" = "LAMP5 GABA",
    "oligodendrocyte precursor cell" = "OPC",
    "L2/3 intratelencephalic projecting glutamatergic neuron" = "L2/3 IT",
    "VIP GABAergic cortical interneuron" = "VIP GABA",
    "astrocyte" = "Astrocyte",
    "vascular leptomeningeal cell" = "VLMC",
    "L2/3-6 intratelencephalic projecting glutamatergic neuron" = "L2/3-6 IT",
    "L6 intratelencephalic projecting glutamatergic neuron" = "L6 IT",
    "pvalb GABAergic cortical interneuron" = "PVALB GABA",
    "sst GABAergic cortical interneuron" = "SST GABA",
    "endothelial cell" = "Endothelial",
    "T cell" = "T cell",
    "GABAergic neuron" = "GABA neuron",
    "L6b glutamatergic neuron of the primary motor cortex" = "L6b Glu",
    "L6 corticothalamic-projecting glutamatergic cortical neuron" = "L6 CT",
    "L5/6 near-projecting glutamatergic neuron" = "L5/6 NP",
    "natural killer cell" = "NK cell",
    "pericyte" = "Pericyte",
    "perivascular macrophage" = "PVM",
    "smooth muscle cell" = "SMC"
  )
  y <- unname(map[as.character(x)])
  y[is.na(y)] <- as.character(x)[is.na(y)]
  y
}

build_umap <- function(obj) {
  if (length(obj@reductions) > 0 && "umap" %in% names(obj@reductions)) {
    return(obj)
  }
  DefaultAssay(obj) <- "RNA"
  if ("meta.features" %in% slotNames(obj[["RNA"]])) {
    obj[["RNA"]]@meta.features <- data.frame(row.names = rownames(obj[["RNA"]]))
  }
  mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  mu <- Matrix::rowMeans(mat)
  mu2 <- Matrix::rowMeans(mat ^ 2)
  vars <- mu2 - mu ^ 2
  vars <- sort(vars, decreasing = TRUE)
  feats <- names(vars)[seq_len(min(2000, length(vars)))]
  obj <- ScaleData(obj, features = feats, verbose = FALSE)
  obj <- RunPCA(obj, features = feats, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, n.neighbors = 30, min.dist = 0.35, seed.use = 20260509, verbose = FALSE)
  obj
}

for (cfg in configs) {
  message("Processing ", cfg$key)
  obj <- readRDS(cfg$rds)
  scores <- fread(cfg$score)
  setnames(scores, old = names(scores)[1], new = "cell_id")
  obj <- obj[, intersect(colnames(obj), scores$cell_id)]
  obj <- build_umap(obj)
  emb <- as.data.frame(Embeddings(obj, "umap"))
  emb$cell_id <- rownames(emb)
  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)
  keep <- c("cell_id", "cell_type", "class", "subclass", "subtype", "disease")
  keep <- intersect(keep, colnames(meta))
  dat <- merge(emb, meta[, keep, drop = FALSE], by = "cell_id", all.x = TRUE)
  dat <- merge(dat, scores[, .(cell_id, scPagwas.TRS.Score, Random_Correct_BG_p, Random_Correct_BG_adjp, Random_Correct_BG_z)], by = "cell_id", all.x = TRUE)
  dat$analysis <- cfg$label
  dat$analysis_key <- cfg$key
  dat$cell_type_short <- short_celltype(dat$cell_type)
  fwrite(dat, file.path(outdir, paste0(cfg$key, "_umap_trs.csv")))
  dat_dt <- as.data.table(dat)
  tab <- dat_dt[, .(
    n = .N,
    median_trs = median(scPagwas.TRS.Score, na.rm = TRUE),
    mean_trs = mean(scPagwas.TRS.Score, na.rm = TRUE),
    min_adj_p = min(Random_Correct_BG_adjp, na.rm = TRUE)
  ), by = .(analysis, analysis_key, cell_type, cell_type_short)]
  fwrite(tab[order(analysis_key, -median_trs)], file.path(outdir, paste0(cfg$key, "_celltype_trs_summary.csv")))
  rm(obj, dat, scores)
  gc()
}
