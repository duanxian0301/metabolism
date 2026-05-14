options(stringsAsFactors = FALSE)
options(timeout = 600)

suppressPackageStartupMessages(library(Seurat))

in_file <- "D:/codex/GenomicSEM/data_100k/MSSM_PD_final.rds"
out_file <- "D:/scPagwas/MSSM_PD_20k.rds"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

pick_celltype_col <- function(meta) {
  candidates <- c("cell_type", "celltype", "CellType", "broad_celltype", "broad.celltype")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) {
    stop("No supported cell-type metadata column found. Available columns: ",
         paste(colnames(meta), collapse = ", "))
  }
  hit[[1]]
}

message("Reading PD object: ", in_file)
obj <- readRDS(in_file)
ct_col <- pick_celltype_col(obj@meta.data)
message("Using cell-type column: ", ct_col)

obj$cell_type <- as.character(obj@meta.data[[ct_col]])
tab <- sort(table(obj$cell_type), decreasing = TRUE)
message("Original cells: ", ncol(obj))
print(head(tab, 25))

set.seed(20260501)
target_total <- 20000L
cells <- colnames(obj)
meta <- obj@meta.data

per_type <- split(rownames(meta), meta$cell_type)
n_by_type <- sapply(per_type, length)
prop_target <- floor(n_by_type / sum(n_by_type) * target_total)
prop_target[prop_target == 0] <- 1L

# Reconcile rounding while respecting per-type maxima.
while (sum(prop_target) > target_total) {
  idx <- which.max(prop_target)
  if (prop_target[idx] > 1) {
    prop_target[idx] <- prop_target[idx] - 1L
  } else {
    break
  }
}
while (sum(prop_target) < target_total) {
  room <- n_by_type - prop_target
  idx <- which.max(room)
  if (room[idx] > 0) {
    prop_target[idx] <- prop_target[idx] + 1L
  } else {
    break
  }
}
prop_target <- pmin(prop_target, n_by_type)

selected <- unlist(
  lapply(names(per_type), function(ct) {
    v <- per_type[[ct]]
    k <- as.integer(prop_target[[ct]])
    if (length(v) <= k) v else sample(v, k)
  }),
  use.names = FALSE
)

selected <- unique(selected)
if (length(selected) > target_total) {
  selected <- sample(selected, target_total)
}

message("Selected cells: ", length(selected))
obj20k <- subset(obj, cells = selected)
obj20k$cell_type <- as.character(obj20k$cell_type)

# Keep a lighter object while preserving RNA assay content needed by scPagwas.
DefaultAssay(obj20k) <- "RNA"
obj20k <- DietSeurat(obj20k, assays = "RNA", dimreducs = NULL, graphs = NULL, misc = TRUE)

saveRDS(obj20k, out_file)
message("Saved: ", out_file)
print(sort(table(obj20k$cell_type), decreasing = TRUE))
