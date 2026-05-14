library(GenomicSEM)
library(data.table)

# -----------------------------------------
# Nonlipid module compact panel:
# 1) read the compact nonlipid trait panel
# 2) reuse existing .sumstats.gz files
# 3) run module-specific multivariate LDSC
# 4) export S / V / rg matrices for downstream EFA/ESEM
# -----------------------------------------

panel_path <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/nonlipid_module_from_full_manifest/compact_panel_review/nonlipid_module_compact_kept.tsv"
gwas_dir <- "D:/metabolic/GWAS"
sumstats_dir <- file.path(gwas_dir, "sumstats")
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
output_dir <- file.path(gwas_dir, "genomicgem_main_zgt4_nonproportion", "step16_ldsc_nonlipid_module_compact15")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c("study_accession", "trait_code", "biomarker_name", "group", "sumstats_file")

if (!file.exists(panel_path)) {
  stop("Compact nonlipid panel not found: ", panel_path)
}
if (!dir.exists(sumstats_dir)) {
  stop("sumstats directory not found: ", sumstats_dir)
}
if (!dir.exists(ld_ref_dir)) {
  stop("LDSC reference directory not found: ", ld_ref_dir)
}

panel_df <- fread(panel_path)
missing_cols <- setdiff(required_cols, names(panel_df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s) in compact nonlipid panel: ",
    paste(missing_cols, collapse = ", ")
  )
}

panel_df <- unique(panel_df[, .(
  study_accession,
  trait_code,
  biomarker_name,
  group,
  proposed_module,
  marker_role,
  sumstats_file
)])

panel_df[, study_accession := as.character(study_accession)]
panel_df[, trait_code := as.character(trait_code)]
panel_df[, sumstats_file := as.character(sumstats_file)]

if (nrow(panel_df) < 2) {
  stop("Need at least two traits to run multivariate LDSC.")
}

missing_sumstats <- panel_df[!file.exists(sumstats_file)]
if (nrow(missing_sumstats) > 0) {
  fwrite(missing_sumstats, file.path(output_dir, "missing_sumstats.tsv"), sep = "\t")
  stop(
    "Missing .sumstats.gz files for ",
    nrow(missing_sumstats),
    " trait(s). See: ",
    file.path(output_dir, "missing_sumstats.tsv")
  )
}

trait_manifest <- copy(panel_df)
fwrite(trait_manifest, file.path(output_dir, "trait_manifest.tsv"), sep = "\t")

message("Running module-specific multivariate LDSC for ", nrow(panel_df), " nonlipid traits...")
ldsc_out <- ldsc(
  traits = panel_df$sumstats_file,
  sample.prev = rep(NA, nrow(panel_df)),
  population.prev = rep(NA, nrow(panel_df)),
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  trait.names = panel_df$trait_code,
  ldsc.log = file.path(output_dir, "nonlipid_module_compact15_multivariate_ldsc")
)

saveRDS(ldsc_out, file.path(output_dir, "nonlipid_module_compact15_multivariate_ldsc.rds"))

S_matrix <- ldsc_out$S
V_matrix <- ldsc_out$V
I_matrix <- ldsc_out$I
rg_matrix <- cov2cor(S_matrix)

write.csv(S_matrix, file.path(output_dir, "nonlipid_module_compact15_S_matrix.csv"), row.names = TRUE)
write.csv(V_matrix, file.path(output_dir, "nonlipid_module_compact15_V_matrix.csv"), row.names = TRUE)
write.csv(I_matrix, file.path(output_dir, "nonlipid_module_compact15_I_matrix.csv"), row.names = TRUE)
write.csv(rg_matrix, file.path(output_dir, "nonlipid_module_compact15_rg_matrix.csv"), row.names = TRUE)

summary_dt <- data.table(
  trait = colnames(S_matrix),
  h2 = diag(S_matrix),
  intercept = diag(I_matrix)
)
fwrite(summary_dt, file.path(output_dir, "nonlipid_module_compact15_ldsc_summary.tsv"), sep = "\t")

message("Finished module-specific multivariate LDSC.")
message("Results written to: ", output_dir)
