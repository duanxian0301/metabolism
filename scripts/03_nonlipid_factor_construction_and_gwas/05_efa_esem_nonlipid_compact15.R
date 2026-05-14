library(GenomicSEM)
library(Matrix)
library(data.table)
library(psych)
library(lavaan)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
panel_dir <- file.path(root_dir, "nonlipid_module_from_full_manifest", "compact_panel_review")
ldsc_all_dir <- file.path(root_dir, "step16_ldsc_nonlipid_module_compact15")
output_dir <- file.path(root_dir, "step17_efa_esem_nonlipid_module_compact15")
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

panel_dt <- fread(file.path(panel_dir, "nonlipid_module_compact_kept.tsv"))
kept_traits <- panel_dt$trait_code

safe_metric <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(x[1])
}

load_all_covariance <- function() {
  S <- as.matrix(
    read.csv(
      file.path(ldsc_all_dir, "nonlipid_module_compact15_S_matrix.csv"),
      row.names = 1,
      check.names = FALSE
    )
  )
  dimnames(S) <- list(kept_traits, kept_traits)
  S <- S[kept_traits, kept_traits]
  eig_min <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 0) {
    S <- as.matrix(nearPD(S, corr = FALSE)$mat)
    dimnames(S) <- list(kept_traits, kept_traits)
  }
  S
}

run_or_load_odd_ldsc <- function() {
  odd_rds <- file.path(output_dir, "nonlipid_module_compact15_ODD_EFA.rds")
  if (file.exists(odd_rds)) {
    message("Reusing existing ODD nonlipid-module LDSC result...")
    return(readRDS(odd_rds))
  }

  message("Running ODD-chromosome LDSC for nonlipid-module EFA...")
  odd_out <- ldsc(
    traits = panel_dt$sumstats_file,
    sample.prev = rep(NA, nrow(panel_dt)),
    population.prev = rep(NA, nrow(panel_dt)),
    ld = ld_ref_dir,
    wld = ld_ref_dir,
    trait.names = panel_dt$trait_code,
    ldsc.log = file.path(output_dir, "nonlipid_module_compact15_ODD_EFA"),
    select = "ODD"
  )
  saveRDS(odd_out, odd_rds)
  odd_out
}

extract_psych_loadings <- function(fa_obj, nfactors) {
  mat <- as.matrix(unclass(fa_obj$loadings))
  mat <- mat[, seq_len(nfactors), drop = FALSE]
  rownames(mat) <- rownames(fa_obj$loadings)
  colnames(mat) <- paste0("F", seq_len(nfactors))
  dt <- as.data.table(mat, keep.rownames = "trait")
  merge(
    dt,
    panel_dt[, .(trait_code, biomarker_name, group, proposed_module, marker_role)],
    by.x = "trait",
    by.y = "trait_code",
    all.x = TRUE
  )
}

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

fit_target_esem <- function(S, loadings_file, nfactors, dataset_label, salient_cutoff = 0.25) {
  factor_cols <- paste0("F", seq_len(nfactors))
  efa_loadings <- fread(loadings_file)
  target <- build_target_matrix(efa_loadings, factor_cols, salient_cutoff = salient_cutoff)

  fit <- efa(
    sample.cov = S,
    sample.nobs = 200,
    nfactors = nfactors,
    rotation = "target",
    rotation.args = list(target = target, rstarts = 50),
    output = "lavaan"
  )

  model_tag <- paste0(dataset_label, "_", nfactors, "factor")
  write.csv(S, file.path(output_dir, paste0(model_tag, "_S_matrix.csv")), row.names = TRUE)
  write.csv(target, file.path(output_dir, paste0(model_tag, "_target_matrix.csv")))
  saveRDS(fit, file.path(output_dir, paste0(model_tag, "_fit.rds")))

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
  fwrite(loading_table, file.path(output_dir, paste0(model_tag, "_loadings.tsv")), sep = "\t")

  class_counts <- loading_table[, .(
    n_salient = sum(abs(std.all) >= salient_cutoff)
  ), by = factor]
  fwrite(class_counts, file.path(output_dir, paste0(model_tag, "_salient_counts.tsv")), sep = "\t")

  data.table(
    dataset = dataset_label,
    nfactors = nfactors,
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

odd_ldsc <- run_or_load_odd_ldsc()
raw_S_odd <- as.matrix(odd_ldsc$S)
smoothed_S_odd <- as.matrix(nearPD(raw_S_odd, corr = FALSE)$mat)
raw_rg_odd <- cov2cor(raw_S_odd)
smoothed_rg_odd <- cov2cor(smoothed_S_odd)

dimnames(raw_S_odd) <- list(kept_traits, kept_traits)
dimnames(smoothed_S_odd) <- list(kept_traits, kept_traits)
dimnames(raw_rg_odd) <- list(kept_traits, kept_traits)
dimnames(smoothed_rg_odd) <- list(kept_traits, kept_traits)

write.csv(raw_S_odd, file.path(output_dir, "nonlipid_module_compact15_ODD_S_matrix_raw.csv"), row.names = TRUE)
write.csv(smoothed_S_odd, file.path(output_dir, "nonlipid_module_compact15_ODD_S_matrix_smoothed.csv"), row.names = TRUE)
write.csv(raw_rg_odd, file.path(output_dir, "nonlipid_module_compact15_ODD_rg_matrix_raw.csv"), row.names = TRUE)
write.csv(smoothed_rg_odd, file.path(output_dir, "nonlipid_module_compact15_ODD_rg_matrix_smoothed.csv"), row.names = TRUE)

S_all <- load_all_covariance()
write.csv(S_all, file.path(output_dir, "nonlipid_module_compact15_ALL_S_matrix_smoothed.csv"), row.names = TRUE)

diagnostics_dt <- data.table(
  metric = c("n_traits", "raw_min_eigenvalue_odd", "smoothed_min_eigenvalue_odd", "smoothed_min_eigenvalue_all"),
  value = c(
    length(kept_traits),
    min(eigen(raw_S_odd, symmetric = TRUE, only.values = TRUE)$values),
    min(eigen(smoothed_S_odd, symmetric = TRUE, only.values = TRUE)$values),
    min(eigen(S_all, symmetric = TRUE, only.values = TRUE)$values)
  )
)
fwrite(diagnostics_dt, file.path(output_dir, "efa_input_matrix_diagnostics.tsv"), sep = "\t")

R_mat <- cov2cor(smoothed_S_odd)
parallel_obj <- fa.parallel(
  R_mat,
  n.obs = 200,
  fm = "minres",
  fa = "fa",
  plot = FALSE,
  error.bars = FALSE
)
eigs <- eigen(R_mat, symmetric = TRUE, only.values = TRUE)$values
criteria_dt <- data.table(
  criterion = c("eigenvalues", "kaiser_n", "fa_parallel_suggested"),
  value = c(
    paste(round(eigs, 6), collapse = "; "),
    sum(eigs > 1),
    parallel_obj$nfact
  )
)
fwrite(criteria_dt, file.path(output_dir, "nonlipid_module_compact15_factor_criteria.tsv"), sep = "\t")

fit_rows <- list()
for (nf in 1:5) {
  efa_fit <- tryCatch(
    fa(
      r = R_mat,
      nfactors = nf,
      rotate = if (nf == 1) "none" else "promax",
      fm = "minres",
      SMC = TRUE
    ),
    error = function(e) e
  )

  if (inherits(efa_fit, "error")) {
    fit_rows[[nf]] <- data.table(
      nfactors = nf,
      objective = NA_real_,
      rms = NA_real_,
      rmsea = NA_real_,
      tli = NA_real_,
      bic = NA_real_,
      fit_off = NA_real_,
      status = paste("failed:", conditionMessage(efa_fit))
    )
    next
  }

  loadings_dt <- extract_psych_loadings(efa_fit, nf)
  fwrite(loadings_dt, file.path(output_dir, sprintf("EFA_minres_%dfactor_loadings.tsv", nf)), sep = "\t")

  uniqueness_dt <- data.table(
    trait = rownames(efa_fit$loadings),
    uniqueness = efa_fit$uniquenesses
  )
  fwrite(uniqueness_dt, file.path(output_dir, sprintf("EFA_minres_%dfactor_uniqueness.tsv", nf)), sep = "\t")

  fit_rows[[nf]] <- data.table(
    nfactors = nf,
    objective = safe_metric(efa_fit$criteria["objective"]),
    rms = safe_metric(efa_fit$rms),
    rmsea = safe_metric(efa_fit$RMSEA),
    tli = safe_metric(efa_fit$TLI),
    bic = safe_metric(efa_fit$BIC),
    fit_off = safe_metric(efa_fit$fit.off),
    status = "ok"
  )
}

fit_summary_dt <- rbindlist(fit_rows, fill = TRUE)
fwrite(fit_summary_dt, file.path(output_dir, "EFA_minres_fit_summary.tsv"), sep = "\t")

esem_results <- list(
  fit_target_esem(
    S = smoothed_S_odd,
    loadings_file = file.path(output_dir, "EFA_minres_2factor_loadings.tsv"),
    nfactors = 2,
    dataset_label = "ODD"
  ),
  fit_target_esem(
    S = smoothed_S_odd,
    loadings_file = file.path(output_dir, "EFA_minres_3factor_loadings.tsv"),
    nfactors = 3,
    dataset_label = "ODD"
  ),
  fit_target_esem(
    S = S_all,
    loadings_file = file.path(output_dir, "EFA_minres_2factor_loadings.tsv"),
    nfactors = 2,
    dataset_label = "ALL"
  ),
  fit_target_esem(
    S = S_all,
    loadings_file = file.path(output_dir, "EFA_minres_3factor_loadings.tsv"),
    nfactors = 3,
    dataset_label = "ALL"
  )
)

summary_dt <- rbindlist(esem_results, fill = TRUE)
fwrite(summary_dt, file.path(output_dir, "nonlipid_module_compact15_esem_target_summary.tsv"), sep = "\t")

message("Nonlipid-module compact15 EFA + target-rotation ESEM finished.")
message("Results written to: ", output_dir)
