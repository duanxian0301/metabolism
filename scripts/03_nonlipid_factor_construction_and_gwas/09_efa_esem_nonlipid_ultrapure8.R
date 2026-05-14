library(GenomicSEM)
library(Matrix)
library(data.table)
library(psych)
library(lavaan)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
panel_dir <- file.path(root_dir, "nonlipid_module_from_full_manifest", "compact_panel_review", "ultra_pure_3factor_review")
ldsc_all_dir <- file.path(root_dir, "step19_ldsc_nonlipid_module_ultrapure8")
output_dir <- file.path(root_dir, "step20_efa_esem_nonlipid_module_ultrapure8")
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

panel_dt <- fread(file.path(panel_dir, "nonlipid_module_ultrapure3_kept.tsv"))
kept_traits <- panel_dt$trait_code

safe_metric <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(x[1])
}

load_all_covariance <- function() {
  S <- as.matrix(
    read.csv(
      file.path(ldsc_all_dir, "nonlipid_module_ultrapure8_S_matrix.csv"),
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
  odd_rds <- file.path(output_dir, "nonlipid_module_ultrapure8_ODD_EFA.rds")
  if (file.exists(odd_rds)) return(readRDS(odd_rds))
  odd_out <- ldsc(
    traits = panel_dt$sumstats_file,
    sample.prev = rep(NA, nrow(panel_dt)),
    population.prev = rep(NA, nrow(panel_dt)),
    ld = ld_ref_dir,
    wld = ld_ref_dir,
    trait.names = panel_dt$trait_code,
    ldsc.log = file.path(output_dir, "nonlipid_module_ultrapure8_ODD_EFA"),
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
    panel_dt[, .(trait_code, biomarker_name, group, ultra_pure_factor, marker_role)],
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
    panel_dt[, .(trait_code, biomarker_name, group, ultra_pure_factor, marker_role)],
    by.x = "trait",
    by.y = "trait_code",
    all.x = TRUE
  )
  loading_table[, abs_std_all := abs(std.all)]
  setorder(loading_table, factor, -abs_std_all, trait)
  fwrite(loading_table, file.path(output_dir, paste0(model_tag, "_loadings.tsv")), sep = "\t")

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

write.csv(raw_S_odd, file.path(output_dir, "nonlipid_module_ultrapure8_ODD_S_matrix_raw.csv"), row.names = TRUE)
write.csv(smoothed_S_odd, file.path(output_dir, "nonlipid_module_ultrapure8_ODD_S_matrix_smoothed.csv"), row.names = TRUE)
write.csv(raw_rg_odd, file.path(output_dir, "nonlipid_module_ultrapure8_ODD_rg_matrix_raw.csv"), row.names = TRUE)
write.csv(smoothed_rg_odd, file.path(output_dir, "nonlipid_module_ultrapure8_ODD_rg_matrix_smoothed.csv"), row.names = TRUE)

S_all <- load_all_covariance()
write.csv(S_all, file.path(output_dir, "nonlipid_module_ultrapure8_ALL_S_matrix_smoothed.csv"), row.names = TRUE)

R_mat <- cov2cor(smoothed_S_odd)
fit_rows <- list()
for (nf in 1:4) {
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
    fit_rows[[nf]] <- data.table(nfactors = nf, status = paste("failed:", conditionMessage(efa_fit)))
    next
  }
  loadings_dt <- extract_psych_loadings(efa_fit, nf)
  fwrite(loadings_dt, file.path(output_dir, sprintf("EFA_minres_%dfactor_loadings.tsv", nf)), sep = "\t")
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

esem_results <- list(
  fit_target_esem(smoothed_S_odd, file.path(output_dir, "EFA_minres_3factor_loadings.tsv"), 3, "ODD"),
  fit_target_esem(S_all, file.path(output_dir, "EFA_minres_3factor_loadings.tsv"), 3, "ALL")
)
fwrite(rbindlist(esem_results, fill = TRUE), file.path(output_dir, "nonlipid_module_ultrapure8_esem_summary.tsv"), sep = "\t")
