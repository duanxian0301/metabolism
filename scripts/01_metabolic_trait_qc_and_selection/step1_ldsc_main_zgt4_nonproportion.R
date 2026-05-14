library(GenomicSEM)
library(data.table)
library(readxl)

# -----------------------------------------
# Metabolic GenomicGEM Step 1
# Main_Zgt4_nonproportion:
# 1) read the selected traits from Excel
# 2) reuse existing .sumstats.gz files
# 3) run multivariate LDSC
# 4) export S / V / rg matrices for downstream EFA/CFA/GEM
# -----------------------------------------

excel_path <- "D:/metabolic/metabolite_FGWAS_selection_lists.xlsx"
sheet_name <- "Main_Zgt4_nonproportion"
gwas_dir <- "D:/metabolic/GWAS"
sumstats_dir <- file.path(gwas_dir, "sumstats")
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
output_dir <- file.path(gwas_dir, "genomicgem_main_zgt4_nonproportion", "step1_ldsc_results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_excel_cols <- c("study_accession", "trait")

if (!file.exists(excel_path)) {
  stop("Excel file not found: ", excel_path)
}
if (!dir.exists(sumstats_dir)) {
  stop("sumstats directory not found: ", sumstats_dir)
}
if (!dir.exists(ld_ref_dir)) {
  stop("LDSC reference directory not found: ", ld_ref_dir)
}

selection_df <- as.data.table(read_excel(excel_path, sheet = sheet_name))
missing_cols <- setdiff(required_excel_cols, names(selection_df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s) in sheet '", sheet_name, "': ",
    paste(missing_cols, collapse = ", ")
  )
}

selection_df <- unique(selection_df[, .(study_accession, trait, `biomarker name`, group)])
selection_df[, study_accession := as.character(study_accession)]
selection_df[, trait := as.character(trait)]

if (nrow(selection_df) < 2) {
  stop("Need at least two traits to run multivariate LDSC.")
}

selection_df[, sumstats_file := file.path(sumstats_dir, paste0(study_accession, ".sumstats.gz"))]
missing_sumstats <- selection_df[!file.exists(sumstats_file)]
if (nrow(missing_sumstats) > 0) {
  fwrite(missing_sumstats, file.path(output_dir, "missing_sumstats.tsv"), sep = "\t")
  stop(
    "Missing .sumstats.gz files for ",
    nrow(missing_sumstats),
    " trait(s). See: ",
    file.path(output_dir, "missing_sumstats.tsv")
  )
}

trait_manifest <- copy(selection_df)
setnames(
  trait_manifest,
  old = c("study_accession", "trait", "biomarker name", "group"),
  new = c("study_accession", "trait_code", "biomarker_name", "group")
)
fwrite(trait_manifest, file.path(output_dir, "trait_manifest.tsv"), sep = "\t")

message("Running multivariate LDSC for ", nrow(selection_df), " selected traits...")
ldsc_out <- ldsc(
  traits = selection_df$sumstats_file,
  sample.prev = rep(NA, nrow(selection_df)),
  population.prev = rep(NA, nrow(selection_df)),
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  trait.names = selection_df$trait,
  ldsc.log = file.path(output_dir, "Main_Zgt4_nonproportion_multivariate_ldsc")
)

saveRDS(ldsc_out, file.path(output_dir, "Main_Zgt4_nonproportion_multivariate_ldsc.rds"))

S_matrix <- ldsc_out$S
V_matrix <- ldsc_out$V
I_matrix <- ldsc_out$I
rg_matrix <- cov2cor(S_matrix)

write.csv(S_matrix, file.path(output_dir, "Main_Zgt4_nonproportion_S_matrix.csv"), row.names = TRUE)
write.csv(V_matrix, file.path(output_dir, "Main_Zgt4_nonproportion_V_matrix.csv"), row.names = TRUE)
write.csv(I_matrix, file.path(output_dir, "Main_Zgt4_nonproportion_I_matrix.csv"), row.names = TRUE)
write.csv(rg_matrix, file.path(output_dir, "Main_Zgt4_nonproportion_rg_matrix.csv"), row.names = TRUE)

summary_dt <- data.table(
  trait = colnames(S_matrix),
  h2 = diag(S_matrix),
  intercept = diag(I_matrix)
)
fwrite(summary_dt, file.path(output_dir, "Main_Zgt4_nonproportion_ldsc_summary.tsv"), sep = "\t")

message("Finished multivariate LDSC.")
message("Results written to: ", output_dir)
