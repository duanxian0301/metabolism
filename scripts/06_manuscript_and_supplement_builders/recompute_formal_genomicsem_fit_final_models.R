.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))
library(GenomicSEM)
library(data.table)

root_dir <- "/mnt/d/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
out_dir <- file.path(root_dir, "manuscript_formal_genomicsem_fit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fit_and_export <- function(name, rds_path, model_syntax, output_prefix) {
  covstruc <- readRDS(rds_path)
  fit <- usermodel(
    covstruc = covstruc,
    model = model_syntax,
    estimation = "DWLS",
    CFIcalc = TRUE,
    std.lv = TRUE
  )

  modelfit_dt <- as.data.table(as.data.frame(fit$modelfit), keep.rownames = "metric")
  fwrite(modelfit_dt, file.path(out_dir, paste0(output_prefix, "_formal_genomicsem_modelfit.tsv")), sep = "\t")
  fwrite(as.data.table(fit$results), file.path(out_dir, paste0(output_prefix, "_formal_genomicsem_results.tsv")), sep = "\t")
  saveRDS(fit, file.path(out_dir, paste0(output_prefix, "_formal_genomicsem_fit.rds")))

  warnings_vec <- warnings()
  if (!is.null(warnings_vec)) {
    writeLines(capture.output(warnings_vec), file.path(out_dir, paste0(output_prefix, "_formal_genomicsem_warnings.txt")))
  }

  invisible(fit)
}

lipid_model <- "
F1 =~ M_HDL_TG + VLDL_size + MUFA + S_VLDL_TG
F2 =~ ApoA1 + HDL_CE
F3 =~ XS_VLDL_FC + VLDL_CE
F1 ~~ F2
F1 ~~ F3
F2 ~~ F3
"

nonlipid_model <- "
F1 =~ Acetoacetate + bOHbutyrate
F2 =~ Val + Leu + Phe
F3 =~ Acetate + Lactate + Glucose
F1 ~~ F2
F1 ~~ F3
F2 ~~ F3
"

fit_and_export(
  name = "lipid_final8",
  rds_path = file.path(root_dir, "step12_ldsc_lipid_module_final8", "lipid_module_final8_multivariate_ldsc.rds"),
  model_syntax = lipid_model,
  output_prefix = "lipid_final8"
)

fit_and_export(
  name = "nonlipid_final8",
  rds_path = file.path(root_dir, "step19_ldsc_nonlipid_module_ultrapure8", "nonlipid_module_ultrapure8_multivariate_ldsc.rds"),
  model_syntax = nonlipid_model,
  output_prefix = "nonlipid_final8"
)
