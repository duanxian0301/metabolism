library(GenomicSEM)
library(data.table)
library(Matrix)
library(lavaan)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
model_dir <- file.path(root_dir, "step13_efa_esem_lipid_module_final8")
panel_path <- file.path(root_dir, "lipid_module_from_full_manifest", "compact_panel_review", "ultra_pure_3factor_review", "final8_review", "lipid_module_ultrapure3_final8_kept.tsv")
out_dir <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/results/22_supplement_tables_lipid8_F2_AD/loo_sensitivity"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

panel_dt <- fread(panel_path)
target_dt <- fread(file.path(model_dir, "EFA_minres_3factor_loadings.tsv"))
S_all <- as.matrix(read.csv(file.path(model_dir, "lipid_module_final8_ALL_S_matrix_smoothed.csv"), row.names = 1, check.names = FALSE))
S_odd <- as.matrix(read.csv(file.path(model_dir, "lipid_module_final8_ODD_S_matrix_smoothed.csv"), row.names = 1, check.names = FALSE))

dataset_list <- list(
  ALL = S_all,
  ODD = S_odd
)

safe_num <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(x[1])
}

fit_one_drop <- function(S, drop_trait, dataset_name) {
  keep_traits <- setdiff(colnames(S), drop_trait)
  S_sub <- S[keep_traits, keep_traits, drop = FALSE]
  target_sub <- as.matrix(target_dt[trait %in% keep_traits, .(F1, F2, F3)])
  rownames(target_sub) <- target_dt[trait %in% keep_traits, trait]
  target_sub <- target_sub[keep_traits, , drop = FALSE]
  target_bin <- ifelse(abs(target_sub) >= 0.25, 1, 0)

  fit <- tryCatch(
    efa(
      sample.cov = S_sub,
      sample.nobs = 200,
      nfactors = 3,
      rotation = "target",
      rotation.args = list(target = target_bin, rstarts = 50),
      output = "lavaan"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.table(
      scenario = paste0("drop_", drop_trait),
      dataset = dataset_name,
      converged = FALSE,
      chisq = NA_real_,
      df = NA_real_,
      p_chisq = NA_real_,
      cfi = NA_real_,
      srmr = NA_real_,
      negative_residuals = NA_integer_,
      boundary_residuals = NA_integer_,
      status = paste("failed:", conditionMessage(fit))
    ))
  }

  pe <- parameterEstimates(fit, standardized = TRUE)
  obs_resid <- pe[pe$op == "~~" & pe$lhs == pe$rhs & !(pe$lhs %in% c("f1", "f2", "f3")), ]
  neg_resid <- sum(obs_resid$est < 0, na.rm = TRUE)
  boundary_resid <- sum(abs(obs_resid$est) < 1e-6, na.rm = TRUE)

  data.table(
    scenario = paste0("drop_", drop_trait),
    dataset = dataset_name,
    converged = isTRUE(lavInspect(fit, "converged")),
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    negative_residuals = neg_resid,
    boundary_residuals = boundary_resid,
    status = "ok"
  )
}

rows <- list()
for (dataset_name in names(dataset_list)) {
  S <- dataset_list[[dataset_name]]
  for (trait in colnames(S)) {
    rows[[length(rows) + 1L]] <- fit_one_drop(S, trait, dataset_name)
  }
}

summary_dt <- rbindlist(rows, fill = TRUE)
fwrite(summary_dt, file.path(out_dir, "final8_loo_sensitivity_summary.tsv"), sep = "\t")

message("LOO sensitivity finished: ", file.path(out_dir, "final8_loo_sensitivity_summary.tsv"))
