options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 600)

suppressPackageStartupMessages(library(scPagwas))
suppressPackageStartupMessages(library(data.table))

gwas_file <- "/mnt/d/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/clean_factor_inputs/lipid8_F2_clean.txt"
sc_rds_file <- "/mnt/d/scPagwas/MSSM_AD_20k.rds"
resource_dir <- "/mnt/d/scPagwas"
out_dir <- "metabolic_lipid8_F2_MSSM_AD"
out_prefix <- "lipid8_F2_MSSM_AD"

setwd(resource_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading GWAS: ", gwas_file)
gwas <- fread(gwas_file)

setnames(
  gwas,
  old = c("SNP", "CHR", "BP", "BETA", "SE", "FREQ"),
  new = c("rsid", "chrom", "pos", "beta", "se", "maf"),
  skip_absent = FALSE
)

need_cols <- c("rsid", "chrom", "pos", "beta", "se", "maf", "A1", "A2", "P", "N")
missing_cols <- setdiff(need_cols, colnames(gwas))
if (length(missing_cols) > 0) {
  stop("GWAS missing columns: ", paste(missing_cols, collapse = ", "))
}

gwas <- gwas[
  !is.na(rsid) &
    !is.na(chrom) &
    !is.na(pos) &
    !is.na(beta) &
    !is.na(se) & se > 0 &
    !is.na(maf) & maf > 0 & maf < 1 &
    !is.na(P) & P > 0 & P <= 1 &
    !is.na(N) & N > 0
]

gwas[, chrom := as.character(chrom)]
gwas[, pos := as.integer(pos)]
gwas[, maf := as.numeric(maf)]
gwas[, beta := as.numeric(beta)]
gwas[, se := as.numeric(se)]

message("GWAS rows after scPagwas QC: ", nrow(gwas))

Genes_by_pathway_kegg <- readRDS("Genes_by_pathway_kegg.rds")
block_annotation_hg37 <- readRDS("block_annotation_hg37.rds")
chrom_ld <- readRDS("chrom_ld.rds")

Pagwas_data <- scPagwas_main(
  Pagwas = NULL,
  gwas_data = as.data.frame(gwas),
  Single_data = sc_rds_file,
  output.prefix = out_prefix,
  output.dirs = out_dir,
  block_annotation = block_annotation_hg37,
  assay = "RNA",
  Pathway_list = Genes_by_pathway_kegg,
  n.cores = 1,
  iters_singlecell = 100,
  chrom_ld = chrom_ld,
  singlecell = FALSE,
  celltype = TRUE
)

saveRDS(Pagwas_data, file = file.path(out_dir, paste0(out_prefix, "_Pagwas_data.rds")))
message("Done. Results saved in: ", file.path(resource_dir, out_dir))
