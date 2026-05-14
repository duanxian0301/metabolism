library(data.table)

tasks <- list(
  list(
    standard = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step14_native_wsl_usergwas_final8_results/merged_lipid_final8/standard_txt/lipid_final8_F1_standard.txt",
    merged = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step14_native_wsl_usergwas_final8_results/merged_lipid_final8/lipid_final8_F1_userGWAS_merged.tsv.gz"
  ),
  list(
    standard = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step14_native_wsl_usergwas_final8_results/merged_lipid_final8/standard_txt/lipid_final8_F2_standard.txt",
    merged = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step14_native_wsl_usergwas_final8_results/merged_lipid_final8/lipid_final8_F2_userGWAS_merged.tsv.gz"
  ),
  list(
    standard = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step14_native_wsl_usergwas_final8_results/merged_lipid_final8/standard_txt/lipid_final8_F3_standard.txt",
    merged = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step14_native_wsl_usergwas_final8_results/merged_lipid_final8/lipid_final8_F3_userGWAS_merged.tsv.gz"
  ),
  list(
    standard = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step22_native_wsl_usergwas_nonlipid8_results/merged_nonlipid_final8/standard_txt/nonlipid_final8_F1_standard.txt",
    merged = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step22_native_wsl_usergwas_nonlipid8_results/merged_nonlipid_final8/nonlipid_final8_F1_userGWAS_merged.tsv.gz"
  ),
  list(
    standard = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step22_native_wsl_usergwas_nonlipid8_results/merged_nonlipid_final8/standard_txt/nonlipid_final8_F2_standard.txt",
    merged = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step22_native_wsl_usergwas_nonlipid8_results/merged_nonlipid_final8/nonlipid_final8_F2_userGWAS_merged.tsv.gz"
  ),
  list(
    standard = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step22_native_wsl_usergwas_nonlipid8_results/merged_nonlipid_final8/standard_txt/nonlipid_final8_F3_standard.txt",
    merged = "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/step22_native_wsl_usergwas_nonlipid8_results/merged_nonlipid_final8/nonlipid_final8_F3_userGWAS_merged.tsv.gz"
  )
)

summary_file <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/factor_txt_gz_finalization_summary.tsv"

finalize_one <- function(standard_file, merged_file) {
  message("Processing standard txt: ", standard_file)
  txt <- fread(standard_file)
  txt_before <- nrow(txt)
  txt_na_before <- txt[is.na(CHR) | is.na(BP), .N]
  txt_clean <- txt[!is.na(CHR) & !is.na(BP)]
  txt_clean[, CHR := as.integer(CHR)]
  txt_clean[, BP := as.integer(BP)]

  txt_tmp <- paste0(standard_file, ".tmp")
  fwrite(txt_clean, txt_tmp, sep = "\t")
  file.remove(standard_file)
  file.rename(txt_tmp, standard_file)

  lookup <- txt_clean[, .(SNP, CHR, BP, MAF = FRQ, A1, A2)]
  setkey(lookup, SNP)

  message("Processing merged gz: ", merged_file)
  merged <- fread(merged_file)
  merged_before <- nrow(merged)
  wrong_prefix <- intersect(c("CHR", "BP", "MAF", "A1", "A2"), names(merged))
  suffix_cols <- setdiff(names(merged), c("SNP", wrong_prefix))
  merged_core <- merged[, c("SNP", suffix_cols), with = FALSE]
  fixed <- lookup[merged_core, on = "SNP", nomatch = 0]
  fixed <- fixed[, c("SNP", "CHR", "BP", "MAF", "A1", "A2", suffix_cols), with = FALSE]

  merged_tmp <- paste0(merged_file, ".tmp")
  fwrite(fixed, merged_tmp, sep = "\t")
  file.remove(merged_file)
  file.rename(merged_tmp, merged_file)

  data.table(
    standard_file = standard_file,
    merged_file = merged_file,
    txt_rows_before = txt_before,
    txt_na_chr_bp_before = txt_na_before,
    txt_rows_after = nrow(txt_clean),
    merged_rows_before = merged_before,
    merged_rows_after = nrow(fixed),
    rows_removed = txt_before - nrow(txt_clean)
  )
}

summary_dt <- rbindlist(
  lapply(tasks, function(x) finalize_one(x$standard, x$merged)),
  use.names = TRUE,
  fill = TRUE
)

fwrite(summary_dt, summary_file, sep = "\t")
message("Wrote finalization summary to: ", summary_file)
