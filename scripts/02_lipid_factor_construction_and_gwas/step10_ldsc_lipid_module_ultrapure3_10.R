library(GenomicSEM)
library(data.table)

panel_path <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/lipid_module_from_full_manifest/compact_panel_review/ultra_pure_3factor_review/lipid_module_ultrapure3_kept.tsv"
gwas_dir <- "D:/metabolic/GWAS"
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
output_dir <- file.path(gwas_dir, "genomicgem_main_zgt4_nonproportion", "step10_ldsc_lipid_module_ultrapure3_10")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

panel_df <- fread(panel_path)
required_cols <- c("study_accession", "trait_code", "sumstats_file", "biomarker_name", "group", "ultra_pure_factor")
missing_cols <- setdiff(required_cols, names(panel_df))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

panel_df <- unique(panel_df[, .(
  study_accession,
  trait_code,
  biomarker_name,
  group,
  ultra_pure_factor,
  sumstats_file
)])

missing_sumstats <- panel_df[!file.exists(sumstats_file)]
if (nrow(missing_sumstats) > 0) {
  fwrite(missing_sumstats, file.path(output_dir, "missing_sumstats.tsv"), sep = "\t")
  stop("Missing sumstats for one or more ultra-pure markers.")
}

fwrite(panel_df, file.path(output_dir, "trait_manifest.tsv"), sep = "\t")

message("Running multivariate LDSC for ultra-pure 3-factor lipid panel...")
ldsc_out <- ldsc(
  traits = panel_df$sumstats_file,
  sample.prev = rep(NA, nrow(panel_df)),
  population.prev = rep(NA, nrow(panel_df)),
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  trait.names = panel_df$trait_code,
  ldsc.log = file.path(output_dir, "lipid_module_ultrapure3_10_multivariate_ldsc")
)

saveRDS(ldsc_out, file.path(output_dir, "lipid_module_ultrapure3_10_multivariate_ldsc.rds"))

S_matrix <- ldsc_out$S
V_matrix <- ldsc_out$V
I_matrix <- ldsc_out$I
rg_matrix <- cov2cor(S_matrix)

write.csv(S_matrix, file.path(output_dir, "lipid_module_ultrapure3_10_S_matrix.csv"), row.names = TRUE)
write.csv(V_matrix, file.path(output_dir, "lipid_module_ultrapure3_10_V_matrix.csv"), row.names = TRUE)
write.csv(I_matrix, file.path(output_dir, "lipid_module_ultrapure3_10_I_matrix.csv"), row.names = TRUE)
write.csv(rg_matrix, file.path(output_dir, "lipid_module_ultrapure3_10_rg_matrix.csv"), row.names = TRUE)

summary_dt <- data.table(
  trait = colnames(S_matrix),
  h2 = diag(S_matrix),
  intercept = diag(I_matrix)
)
fwrite(summary_dt, file.path(output_dir, "lipid_module_ultrapure3_10_ldsc_summary.tsv"), sep = "\t")

message("Finished multivariate LDSC.")
message("Results written to: ", output_dir)
