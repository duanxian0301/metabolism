library(GenomicSEM)
library(Matrix)
library(data.table)
library(psych)
library(lavaan)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
panel_dir <- file.path(root_dir, "lipid_module_from_full_manifest", "compact_panel_review", "ultra_pure_3factor_review")
ldsc_all_dir <- file.path(root_dir, "step10_ldsc_lipid_module_ultrapure3_10")
output_dir <- file.path(root_dir, "step11_efa_esem_lipid_module_ultrapure3_10")
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

panel_dt <- fread(file.path(panel_dir, "lipid_module_ultrapure3_kept.tsv"))
kept_traits <- panel_dt$trait_code

safe_metric <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(x[1])
}

run_or_load_odd_ldsc <- function() {
  odd_rds <- file.path(output_dir, "lipid_module_ultrapure3_10_ODD_EFA.rds")
  if (file.exists(odd_rds)) return(readRDS(odd_rds))
  odd_out <- ldsc(
    traits = panel_dt$sumstats_file,
    sample.prev = rep(NA, nrow(panel_dt)),
    population.prev = rep(NA, nrow(panel_dt)),
    ld = ld_ref_dir,
    wld = ld_ref_dir,
    trait.names = panel_dt$trait_code,
    ldsc.log = file.path(output_dir, "lipid_module_ultrapure3_10_ODD_EFA"),
    select = "ODD"
  )
  saveRDS(odd_out, odd_rds)
  odd_out
}

load_all_covariance <- function() {
  S <- as.matrix(read.csv(file.path(ldsc_all_dir, "lipid_module_ultrapure3_10_S_matrix.csv"), row.names = 1, check.names = FALSE))
  dimnames(S) <- list(kept_traits, kept_traits)
  if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    S <- as.matrix(nearPD(S, corr = FALSE)$mat)
    dimnames(S) <- list(kept_traits, kept_traits)
  }
  S
}

extract_psych_loadings <- function(fa_obj, nfactors) {
  mat <- as.matrix(unclass(fa_obj$loadings))
  mat <- mat[, seq_len(nfactors), drop = FALSE]
  rownames(mat) <- rownames(fa_obj$loadings)
  colnames(mat) <- paste0("F", seq_len(nfactors))
  dt <- as.data.table(mat, keep.rownames = "trait")
  merge(
    dt,
    panel_dt[, .(trait_code, biomarker_name, group, ultra_pure_factor)],
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

fit_target_esem <- function(S, loadings_file, dataset_label, salient_cutoff = 0.25) {
  efa_loadings <- fread(loadings_file)
  factor_cols <- c("F1", "F2", "F3")
  target <- build_target_matrix(efa_loadings, factor_cols, salient_cutoff = salient_cutoff)
  fit <- efa(
    sample.cov = S,
    sample.nobs = 200,
    nfactors = 3,
    rotation = "target",
    rotation.args = list(target = target, rstarts = 50),
    output = "lavaan"
  )

  tag <- paste0(dataset_label, "_3factor")
  write.csv(S, file.path(output_dir, paste0(tag, "_S_matrix.csv")), row.names = TRUE)
  write.csv(target, file.path(output_dir, paste0(tag, "_target_matrix.csv")))
  saveRDS(fit, file.path(output_dir, paste0(tag, "_fit.rds")))

  loading_table <- extract_loading_table(fit)
  loading_table <- merge(
    loading_table,
    panel_dt[, .(trait_code, biomarker_name, group, ultra_pure_factor)],
    by.x = "trait",
    by.y = "trait_code",
    all.x = TRUE
  )
  loading_table[, abs_std_all := abs(std.all)]
  setorder(loading_table, factor, -abs_std_all, trait)
  fwrite(loading_table, file.path(output_dir, paste0(tag, "_loadings.tsv")), sep = "\t")

  data.table(
    dataset = dataset_label,
    nfactors = 3,
    converged = lavInspect(fit, "converged"),
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    rmsea = fitMeasures(fit, "rmsea"),
    tli = fitMeasures(fit, "tli"),
    aic = fitMeasures(fit, "aic"),
    bic = fitMeasures(fit, "bic")
  )
}

odd_ldsc <- run_or_load_odd_ldsc()
raw_S_odd <- as.matrix(odd_ldsc$S)
smoothed_S_odd <- as.matrix(nearPD(raw_S_odd, corr = FALSE)$mat)
dimnames(raw_S_odd) <- list(kept_traits, kept_traits)
dimnames(smoothed_S_odd) <- list(kept_traits, kept_traits)
S_all <- load_all_covariance()

write.csv(raw_S_odd, file.path(output_dir, "lipid_module_ultrapure3_10_ODD_S_matrix_raw.csv"), row.names = TRUE)
write.csv(smoothed_S_odd, file.path(output_dir, "lipid_module_ultrapure3_10_ODD_S_matrix_smoothed.csv"), row.names = TRUE)
write.csv(S_all, file.path(output_dir, "lipid_module_ultrapure3_10_ALL_S_matrix_smoothed.csv"), row.names = TRUE)

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
parallel_obj <- fa.parallel(R_mat, n.obs = 200, fm = "minres", fa = "fa", plot = FALSE, error.bars = FALSE)
eigs <- eigen(R_mat, symmetric = TRUE, only.values = TRUE)$values
criteria_dt <- data.table(
  criterion = c("eigenvalues", "kaiser_n", "fa_parallel_suggested"),
  value = c(paste(round(eigs, 6), collapse = "; "), sum(eigs > 1), parallel_obj$nfact)
)
fwrite(criteria_dt, file.path(output_dir, "lipid_module_ultrapure3_10_factor_criteria.tsv"), sep = "\t")

fit_rows <- list()
for (nf in 1:4) {
  efa_fit <- tryCatch(
    fa(r = R_mat, nfactors = nf, rotate = if (nf == 1) "none" else "promax", fm = "minres", SMC = TRUE),
    error = function(e) e
  )
  if (inherits(efa_fit, "error")) {
    fit_rows[[nf]] <- data.table(nfactors = nf, objective = NA_real_, rms = NA_real_, rmsea = NA_real_, tli = NA_real_, bic = NA_real_, fit_off = NA_real_, status = paste("failed:", conditionMessage(efa_fit)))
    next
  }
  fwrite(extract_psych_loadings(efa_fit, nf), file.path(output_dir, sprintf("EFA_minres_%dfactor_loadings.tsv", nf)), sep = "\t")
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
fwrite(rbindlist(fit_rows, fill = TRUE), file.path(output_dir, "EFA_minres_fit_summary.tsv"), sep = "\t")

results <- list(
  fit_target_esem(smoothed_S_odd, file.path(output_dir, "EFA_minres_3factor_loadings.tsv"), "ODD"),
  fit_target_esem(S_all, file.path(output_dir, "EFA_minres_3factor_loadings.tsv"), "ALL")
)
fwrite(rbindlist(results, fill = TRUE), file.path(output_dir, "lipid_module_ultrapure3_10_esem_summary.tsv"), sep = "\t")

message("Ultra-pure 3-factor lipid EFA + ESEM finished.")
message("Results written to: ", output_dir)
