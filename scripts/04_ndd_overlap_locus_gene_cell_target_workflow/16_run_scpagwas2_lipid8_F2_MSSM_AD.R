options(stringsAsFactors = FALSE)
options(timeout = 600)
Sys.setenv(OPENBLAS_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(scPagwas)
  library(Seurat)
  library(data.table)
})

root_out <- "/mnt/d/scPagwas/metabolic_scpagwas2"
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)
setwd(root_out)

resource_dir <- "/mnt/d/scPagwas"
genes_by_pathway_kegg <- readRDS(file.path(resource_dir, "Genes_by_pathway_kegg.rds"))
block_annotation_hg37 <- readRDS(file.path(resource_dir, "block_annotation_hg37.rds"))
chrom_ld <- readRDS(file.path(resource_dir, "chrom_ld.rds"))

gwas_file <- "/mnt/d/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/clean_factor_inputs/lipid8_F2_clean.txt"
sc_file <- "/mnt/d/scPagwas/MSSM_AD_20k.rds"
prefix <- "lipid8_F2_MSSM_AD"

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
  stop("Missing cell_type metadata in MSSM_AD_20k.rds")
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
