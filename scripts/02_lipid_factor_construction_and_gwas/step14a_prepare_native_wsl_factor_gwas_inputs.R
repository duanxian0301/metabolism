options(stringsAsFactors = FALSE)
.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))

library(GenomicSEM)
library(data.table)

root_dir <- "/mnt/d/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
raw_gwas_dir <- "/mnt/d/metabolic/GWAS"
panel_file <- file.path(root_dir, "final_model_lipid_module_final8", "final8_trait_manifest.tsv")
output_dir <- file.path(root_dir, "step14_native_wsl_factor_gwas_inputs")
input_dir <- file.path(output_dir, "input")
capped_dir <- file.path(output_dir, "capped_raw")
p_floor <- 1e-300

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(capped_dir, showWarnings = FALSE, recursive = TRUE)

panel_dt <- fread(panel_file)
if (!nrow(panel_dt)) {
  stop("No traits found in final8 panel file: ", panel_file)
}

build_reference_from_raw <- function(panel_dt, n_ref_max = NULL) {
  accession <- panel_dt$study_accession[[1]]
  in_file <- file.path(raw_gwas_dir, paste0(accession, ".txt"))
  ref_dt <- fread(
    in_file,
    select = c("SNP", "A1", "A2", "MAF"),
    nrows = if (is.null(n_ref_max)) Inf else n_ref_max,
    showProgress = FALSE
  )
  ref_dt <- unique(ref_dt[, .(SNP, A1, A2, MAF)])
  ref_dt[!is.na(MAF) & MAF > 0 & MAF < 1]
}

cap_extreme_pvalues <- function(in_file, out_file, p_floor = 1e-300) {
  dt <- fread(in_file, showProgress = FALSE)
  if (!"P" %in% names(dt)) {
    stop("Missing P column in raw GWAS file: ", in_file)
  }
  n_capped <- sum(!is.na(dt$P) & dt$P < p_floor)
  if (n_capped > 0) {
    dt[P < p_floor, P := p_floor]
  }
  fwrite(dt, out_file, sep = "\t")
  n_capped
}

message("Building internal reference file from raw GWAS...")
reference_dt <- build_reference_from_raw(panel_dt)
ref_file <- file.path(output_dir, "final8_internal_reference.tsv.gz")
fwrite(reference_dt, ref_file, sep = "\t")

files <- file.path(raw_gwas_dir, paste0(panel_dt$study_accession, ".txt"))
trait_names <- panel_dt$trait_code

message("Preparing raw GWAS files with p-value floor = ", format(p_floor, scientific = TRUE), " ...")
capped_files <- file.path(capped_dir, paste0(panel_dt$study_accession, "_p_capped.txt"))
capped_counts <- integer(length(files))
for (i in seq_along(files)) {
  capped_counts[i] <- cap_extreme_pvalues(files[i], capped_files[i], p_floor = p_floor)
}

capped_diag <- data.table(
  study_accession = panel_dt$study_accession,
  trait_code = trait_names,
  raw_file = files,
  capped_file = capped_files,
  n_pvalues_capped = capped_counts
)
fwrite(capped_diag, file.path(output_dir, "pvalue_capping_diagnostics.tsv"), sep = "\t")

available_cores <- parallel::detectCores(logical = TRUE)
use_cores <- 1L
chunk_size <- 50000L

message("Running GenomicSEM::sumstats() under WSL...")
sumstats_obj <- sumstats(
  files = capped_files,
  ref = ref_file,
  trait.names = trait_names,
  se.logit = rep(FALSE, length(files)),
  OLS = rep(FALSE, length(files)),
  linprob = rep(FALSE, length(files)),
  N = rep(NA, length(files)),
  betas = NULL,
  info.filter = 0.6,
  maf.filter = 0.01,
  keep.indel = FALSE,
  parallel = FALSE,
  cores = use_cores
)

saveRDS(sumstats_obj, file.path(output_dir, "final8_factorGWAS_sumstats.rds"))
fwrite(as.data.frame(sumstats_obj), file.path(output_dir, "final8_factorGWAS_sumstats.tsv.gz"), sep = "\t")

n_rows <- nrow(sumstats_obj)
chunk_index <- rep(seq_len(ceiling(n_rows / chunk_size)), each = chunk_size)[seq_len(n_rows)]
sumstats_list <- split(sumstats_obj, chunk_index)

chunk_manifest <- data.frame(
  chunk_id = seq_along(sumstats_list),
  file_name = paste0("GenSem_sub_", seq_along(sumstats_list), "_final8_lipid3F.tsv"),
  n_rows = vapply(sumstats_list, nrow, integer(1))
)

for (i in seq_along(sumstats_list)) {
  out_file <- file.path(input_dir, chunk_manifest$file_name[i])
  fwrite(sumstats_list[[i]], file = out_file, sep = "\t")
}

fwrite(chunk_manifest, file.path(output_dir, "chunk_manifest.tsv"), sep = "\t")

diag_dt <- data.table(
  metric = c("n_traits", "n_reference_snps", "n_sumstats_rows", "n_chunks", "cores_used", "p_floor", "n_total_pvalues_capped"),
  value = c(nrow(panel_dt), nrow(reference_dt), n_rows, nrow(chunk_manifest), use_cores, p_floor, sum(capped_counts))
)
fwrite(diag_dt, file.path(output_dir, "prepare_diagnostics.tsv"), sep = "\t")

message("Native WSL factor GWAS input preparation finished.")
