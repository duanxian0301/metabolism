library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
native_dir <- file.path(root_dir, "step14_native_wsl_usergwas_final8_results")
merged_dir <- file.path(native_dir, "merged_lipid_final8")

dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

merge_one <- function(subdir, outfile) {
  files <- list.files(file.path(native_dir, subdir), pattern = "\\.tsv$", full.names = TRUE)
  files <- sort(files)
  if (!length(files)) {
    message("No files found in ", file.path(native_dir, subdir), "; skipping.")
    return(invisible(NULL))
  }
  dt <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
  fwrite(dt, outfile, sep = "\t")
}

merge_one("F1", file.path(merged_dir, "lipid_final8_F1_userGWAS_merged.tsv.gz"))
merge_one("F2", file.path(merged_dir, "lipid_final8_F2_userGWAS_merged.tsv.gz"))
merge_one("F3", file.path(merged_dir, "lipid_final8_F3_userGWAS_merged.tsv.gz"))
merge_one("Q_SNP", file.path(merged_dir, "lipid_final8_QSNP_userGWAS_merged.tsv.gz"))

