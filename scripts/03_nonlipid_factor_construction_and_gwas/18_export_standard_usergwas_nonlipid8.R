library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
native_dir <- file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results")
merged_dir <- file.path(native_dir, "merged_nonlipid_final8")
std_dir <- file.path(merged_dir, "standard_txt")
sumstats_lookup_file <- file.path(root_dir, "step21_native_wsl_factor_gwas_inputs_nonlipid8", "nonlipid8_factorGWAS_sumstats.tsv.gz")

dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(std_dir, showWarnings = FALSE, recursive = TRUE)

factor_files <- c(
  F1 = file.path(merged_dir, "nonlipid_final8_F1_userGWAS_merged.tsv.gz"),
  F2 = file.path(merged_dir, "nonlipid_final8_F2_userGWAS_merged.tsv.gz"),
  F3 = file.path(merged_dir, "nonlipid_final8_F3_userGWAS_merged.tsv.gz")
)

missing_files <- factor_files[!file.exists(factor_files)]
if (length(missing_files) > 0) {
  stop("Merged userGWAS files are required before export:\n", paste(missing_files, collapse = "\n"))
}
if (!file.exists(sumstats_lookup_file)) {
  stop("Missing factor GWAS SNP lookup file: ", sumstats_lookup_file)
}

n_value <- 599249L
lookup_dt <- fread(sumstats_lookup_file, select = c("SNP", "A1", "A2", "MAF"))
lookup_dt <- unique(lookup_dt, by = "SNP")

export_one <- function(infile, outfile, n_value) {
  dt <- fread(infile)
  dt <- merge(lookup_dt, dt[, .(SNP, est, SE, Pval_Estimate)], by = "SNP", all.x = FALSE, all.y = TRUE)
  out <- data.table(
    SNP = dt$SNP,
    A1 = dt$A1,
    A2 = dt$A2,
    FRQ = dt$MAF,
    BETA = dt$est,
    SE = dt$SE,
    P = dt$Pval_Estimate,
    N = as.integer(n_value)
  )
  fwrite(out, outfile, sep = "\t")
}

export_one(
  factor_files[["F1"]],
  file.path(std_dir, "nonlipid_final8_F1_standard.txt"),
  n_value
)

export_one(
  factor_files[["F2"]],
  file.path(std_dir, "nonlipid_final8_F2_standard.txt"),
  n_value
)

export_one(
  factor_files[["F3"]],
  file.path(std_dir, "nonlipid_final8_F3_standard.txt"),
  n_value
)

meta <- data.table(
  factor = c("F1", "F2", "F3"),
  merged_file = unname(factor_files),
  standard_file = c(
    file.path(std_dir, "nonlipid_final8_F1_standard.txt"),
    file.path(std_dir, "nonlipid_final8_F2_standard.txt"),
    file.path(std_dir, "nonlipid_final8_F3_standard.txt")
  ),
  n_used = n_value,
  notes = "Native WSL GenomicSEM userGWAS rerun on nonlipid final8 3-factor model; FREQ uses MAF."
)

q_path <- file.path(merged_dir, "nonlipid_final8_QSNP_userGWAS_merged.tsv.gz")
q_status <- if (file.exists(q_path)) "present" else "not_generated"
meta[, q_snp_status := q_status]

fwrite(meta, file.path(std_dir, "nonlipid_final8_standard_txt_metadata.tsv"), sep = "\t")
