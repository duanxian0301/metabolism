library(GenomicSEM)
library(data.table)

panel_path <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/nonlipid_module_from_full_manifest/compact_panel_review/ultra_pure_3factor_review/nonlipid_module_ultrapure3_kept.tsv"
gwas_dir <- "D:/metabolic/GWAS"
sumstats_dir <- file.path(gwas_dir, "sumstats")
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
output_dir <- file.path(gwas_dir, "genomicgem_main_zgt4_nonproportion", "step19_ldsc_nonlipid_module_ultrapure8")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c("study_accession", "trait_code", "biomarker_name", "group", "sumstats_file")
panel_df <- fread(panel_path)
missing_cols <- setdiff(required_cols, names(panel_df))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

panel_df <- unique(panel_df[, .(
  study_accession,
  trait_code,
  biomarker_name,
  group,
  ultra_pure_factor,
  marker_role,
  sumstats_file
)])

missing_sumstats <- panel_df[!file.exists(sumstats_file)]
if (nrow(missing_sumstats) > 0) {
  fwrite(missing_sumstats, file.path(output_dir, "missing_sumstats.tsv"), sep = "\t")
  stop("Missing sumstats file(s) for ultrapure8 panel.")
}

fwrite(panel_df, file.path(output_dir, "trait_manifest.tsv"), sep = "\t")

ldsc_out <- ldsc(
  traits = panel_df$sumstats_file,
  sample.prev = rep(NA, nrow(panel_df)),
  population.prev = rep(NA, nrow(panel_df)),
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  trait.names = panel_df$trait_code,
  ldsc.log = file.path(output_dir, "nonlipid_module_ultrapure8_multivariate_ldsc")
)

saveRDS(ldsc_out, file.path(output_dir, "nonlipid_module_ultrapure8_multivariate_ldsc.rds"))
write.csv(ldsc_out$S, file.path(output_dir, "nonlipid_module_ultrapure8_S_matrix.csv"), row.names = TRUE)
write.csv(ldsc_out$V, file.path(output_dir, "nonlipid_module_ultrapure8_V_matrix.csv"), row.names = TRUE)
write.csv(ldsc_out$I, file.path(output_dir, "nonlipid_module_ultrapure8_I_matrix.csv"), row.names = TRUE)
write.csv(cov2cor(ldsc_out$S), file.path(output_dir, "nonlipid_module_ultrapure8_rg_matrix.csv"), row.names = TRUE)

summary_dt <- data.table(
  trait = colnames(ldsc_out$S),
  h2 = diag(ldsc_out$S),
  intercept = diag(ldsc_out$I)
)
fwrite(summary_dt, file.path(output_dir, "nonlipid_module_ultrapure8_ldsc_summary.tsv"), sep = "\t")
