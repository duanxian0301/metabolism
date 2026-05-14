suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

set.seed(123)

project_root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
source_rds <- "D:/codex/GenomicSEM/data_100k/MSSM_PD_final.rds"
out_root <- file.path(project_root, "work", "knk_inputs", "MSSM_PD")
manifest_path <- file.path(project_root, "work", "knk_inputs", "pd_knk_input_manifest.csv")

targets <- list(
  Pericyte = list(cell_types = c("pericyte"), max_cells = 4000L),
  Oligodendrocyte_precursor_cell = list(cell_types = c("oligodendrocyte precursor cell"), max_cells = 4000L),
  VIP_GABAergic_cortical_interneuron = list(cell_types = c("VIP GABAergic cortical interneuron"), max_cells = 4000L)
)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

message("[prepare_pd_knk] reading ", source_rds)
obj <- readRDS(source_rds)
meta <- obj@meta.data
meta$cell_type <- as.character(meta$cell_type)

get_counts_mat <- function(x) {
  tryCatch(
    GetAssayData(x, layer = "counts"),
    error = function(e) GetAssayData(x, slot = "counts")
  )
}

counts <- get_counts_mat(obj)
feature_ids <- rownames(counts)
feature_symbols <- feature_ids

manifest_rows <- list()
for (group_name in names(targets)) {
  spec <- targets[[group_name]]
  group_dir <- file.path(out_root, group_name)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

  keep_cells <- rownames(meta)[meta$cell_type %in% spec$cell_types]
  available_n <- length(keep_cells)
  sampled_cells <- keep_cells
  if (available_n > spec$max_cells) {
    sampled_cells <- sample(keep_cells, size = spec$max_cells, replace = FALSE)
  }
  sampled_n <- length(sampled_cells)

  status <- "OK"
  note <- ""
  if (available_n == 0) {
    status <- "SKIP_NO_CELLTYPE"
    note <- "No matching cells in MSSM_PD_final.rds"
  } else {
    group_counts <- counts[, sampled_cells, drop = FALSE]
    group_meta <- meta[sampled_cells, , drop = FALSE]
    Matrix::writeMM(group_counts, file.path(group_dir, "matrix.mtx"))
    fwrite(data.table(barcode = colnames(group_counts)), file.path(group_dir, "barcodes.tsv"), sep = "\t", col.names = FALSE)
    fwrite(data.table(feature_id = feature_ids, feature_name = feature_symbols), file.path(group_dir, "features.tsv"), sep = "\t", col.names = FALSE)
    fwrite(as.data.table(group_meta, keep.rownames = "barcode"), file.path(group_dir, "metadata.csv"))
  }

  manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
    dataset = "MSSM_PD",
    group = group_name,
    cell_types = paste(spec$cell_types, collapse = "|"),
    available_cells = available_n,
    sampled_cells = sampled_n,
    status = status,
    note = note
  )
  message("[prepare_pd_knk] ", group_name, ": available=", available_n, ", sampled=", sampled_n, ", status=", status)
}

manifest <- rbindlist(manifest_rows, fill = TRUE)
fwrite(manifest, manifest_path)
message("[prepare_pd_knk] wrote manifest: ", manifest_path)
