library(data.table)
library(lavaan)
library(GenomicSEM)
library(Matrix)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
step8_dir <- file.path(root_dir, "step8_efa_esem_lipid_module_compact20")
panel_dir <- file.path(root_dir, "lipid_module_from_full_manifest", "compact_panel_review")
output_dir <- file.path(root_dir, "step9_esem_lipid_module_compact20_2factor")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

panel_dt <- fread(file.path(panel_dir, "lipid_module_compact_kept.tsv"))

extract_loading_table <- function(fit) {
  pe <- parameterEstimates(fit, standardized = TRUE)
  loadings <- pe[pe$op == "=~", c("lhs", "rhs", "est", "se", "z", "pvalue", "std.all")]
  setDT(loadings)
  setnames(loadings, c("lhs", "rhs"), c("factor", "trait"))
  loadings
}

build_target_matrix <- function(loadings_dt, factor_cols, salient_cutoff = 0.25) {
  target <- ifelse(abs(as.matrix(loadings_dt[, ..factor_cols])) >= salient_cutoff, 1, 0)
  colnames(target) <- factor_cols
  rownames(target) <- loadings_dt$trait
  target
}

fit_target_esem <- function(S, loadings_file, dataset_label, salient_cutoff = 0.25) {
  efa_loadings <- fread(loadings_file)
  factor_cols <- c("F1", "F2")
  target <- build_target_matrix(efa_loadings, factor_cols, salient_cutoff = salient_cutoff)

  fit <- efa(
    sample.cov = S,
    sample.nobs = 200,
    nfactors = 2,
    rotation = "target",
    rotation.args = list(target = target, rstarts = 50),
    output = "lavaan"
  )

  tag <- paste0(dataset_label, "_2factor")
  write.csv(S, file.path(output_dir, paste0(tag, "_S_matrix.csv")), row.names = TRUE)
  write.csv(target, file.path(output_dir, paste0(tag, "_target_matrix.csv")))
  saveRDS(fit, file.path(output_dir, paste0(tag, "_fit.rds")))

  loading_table <- extract_loading_table(fit)
  loading_table <- merge(
    loading_table,
    panel_dt[, .(trait_code, biomarker_name, group, proposed_module, marker_role)],
    by.x = "trait",
    by.y = "trait_code",
    all.x = TRUE
  )
  loading_table[, abs_std_all := abs(std.all)]
  setorder(loading_table, factor, -abs_std_all, trait)
  fwrite(loading_table, file.path(output_dir, paste0(tag, "_loadings.tsv")), sep = "\t")

  class_counts <- loading_table[, .(
    n_salient = sum(abs(std.all) >= salient_cutoff)
  ), by = factor]
  fwrite(class_counts, file.path(output_dir, paste0(tag, "_salient_counts.tsv")), sep = "\t")

  data.table(
    dataset = dataset_label,
    nfactors = 2,
    converged = lavInspect(fit, "converged"),
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    rmsea = fitMeasures(fit, "rmsea"),
    tli = fitMeasures(fit, "tli"),
    aic = fitMeasures(fit, "aic"),
    bic = fitMeasures(fit, "bic"),
    salient_cutoff = salient_cutoff
  )
}

S_odd <- as.matrix(
  read.csv(
    file.path(step8_dir, "lipid_module_compact20_ODD_S_matrix_smoothed.csv"),
    row.names = 1,
    check.names = FALSE
  )
)
S_all <- as.matrix(
  read.csv(
    file.path(step8_dir, "lipid_module_compact20_ALL_S_matrix_smoothed.csv"),
    row.names = 1,
    check.names = FALSE
  )
)

results <- list(
  fit_target_esem(
    S = S_odd,
    loadings_file = file.path(step8_dir, "EFA_minres_2factor_loadings.tsv"),
    dataset_label = "ODD",
    salient_cutoff = 0.25
  ),
  fit_target_esem(
    S = S_all,
    loadings_file = file.path(step8_dir, "EFA_minres_2factor_loadings.tsv"),
    dataset_label = "ALL",
    salient_cutoff = 0.25
  )
)

summary_dt <- rbindlist(results, fill = TRUE)
fwrite(summary_dt, file.path(output_dir, "lipid_module_compact20_2factor_esem_summary.tsv"), sep = "\t")

message("Lipid-module compact20 2-factor target-rotation ESEM finished.")
message("Results written to: ", output_dir)
