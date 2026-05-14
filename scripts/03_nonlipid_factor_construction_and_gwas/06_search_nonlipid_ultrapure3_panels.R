library(GenomicSEM)
library(Matrix)
library(data.table)
library(psych)
library(lavaan)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
compact_dir <- file.path(root_dir, "step17_efa_esem_nonlipid_module_compact15")
output_dir <- file.path(root_dir, "step18_nonlipid_ultrapure3_search")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_cov <- function(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}

S_all_full <- read_cov(file.path(compact_dir, "nonlipid_module_compact15_ALL_S_matrix_smoothed.csv"))
S_odd_full <- read_cov(file.path(compact_dir, "nonlipid_module_compact15_ODD_S_matrix_smoothed.csv"))
meta <- fread(file.path(root_dir, "nonlipid_module_from_full_manifest", "compact_panel_review", "nonlipid_module_compact_kept.tsv"))
meta <- unique(meta[, .(trait_code, biomarker_name, group, proposed_module, marker_role)])

amino_sets <- list(
  c("Val", "Leu", "Phe"),
  c("Val", "Leu", "Tyr"),
  c("Val", "Leu", "Phe", "Tyr"),
  c("Val", "Leu", "Phe", "His"),
  c("Val", "Leu", "Phe", "Gln"),
  c("Val", "Leu", "Phe", "Tyr", "His")
)

ketone_sets <- list(
  c("Acetoacetate", "bOHbutyrate", "Acetate"),
  c("Acetoacetate", "bOHbutyrate"),
  c("Acetoacetate", "Acetate", "bOHbutyrate", "Lactate")
)

bridge_sets <- list(
  c("Lactate", "GlycA"),
  c("Lactate", "Glucose"),
  c("Lactate", "Albumin"),
  c("Lactate", "Glucose", "GlycA"),
  c("Lactate", "Albumin", "GlycA"),
  c("Lactate", "Glucose", "Albumin"),
  c("Lactate", "Glucose", "Albumin", "GlycA"),
  c("GlycA", "Glucose", "Albumin"),
  c("Lactate", "Citrate", "Glucose")
)

safe_metric <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(x[1])
}

extract_psych_loadings <- function(fa_obj, nfactors, traits) {
  mat <- as.matrix(unclass(fa_obj$loadings))
  mat <- mat[, seq_len(nfactors), drop = FALSE]
  rownames(mat) <- traits
  colnames(mat) <- paste0("F", seq_len(nfactors))
  as.data.table(mat, keep.rownames = "trait")
}

build_target_matrix <- function(loadings_dt, factor_cols, salient_cutoff = 0.25) {
  target <- ifelse(abs(as.matrix(loadings_dt[, ..factor_cols])) >= salient_cutoff, 1, 0)
  colnames(target) <- factor_cols
  rownames(target) <- loadings_dt$trait
  target
}

fit_esem_from_subset <- function(S, traits, nfactors = 3, salient_cutoff = 0.25) {
  S_sub <- S[traits, traits, drop = FALSE]
  eig_min <- min(eigen(S_sub, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 0) {
    S_sub <- as.matrix(nearPD(S_sub, corr = FALSE)$mat)
    dimnames(S_sub) <- list(traits, traits)
  }
  R_mat <- cov2cor(S_sub)
  efa_fit <- fa(
    r = R_mat,
    nfactors = nfactors,
    rotate = "promax",
    fm = "minres",
    SMC = TRUE
  )
  loadings_dt <- extract_psych_loadings(efa_fit, nfactors, traits)
  factor_cols <- paste0("F", seq_len(nfactors))
  target <- build_target_matrix(loadings_dt, factor_cols, salient_cutoff = salient_cutoff)
  fit <- efa(
    sample.cov = S_sub,
    sample.nobs = 200,
    nfactors = nfactors,
    rotation = "target",
    rotation.args = list(target = target, rstarts = 30),
    output = "lavaan"
  )
  list(
    fit = fit,
    S = S_sub,
    target = target
  )
}

results <- list()
ix <- 0L

for (amino in amino_sets) {
  for (ketone in ketone_sets) {
    for (bridge in bridge_sets) {
      traits <- unique(c(amino, ketone, bridge))
      if (length(traits) < 7 || length(traits) > 10) next
      ix <- ix + 1L

      odd_res <- tryCatch(fit_esem_from_subset(S_odd_full, traits), error = function(e) e)
      all_res <- tryCatch(fit_esem_from_subset(S_all_full, traits), error = function(e) e)

      if (inherits(odd_res, "error") || inherits(all_res, "error")) {
        results[[ix]] <- data.table(
          panel_id = sprintf("panel_%03d", ix),
          traits = paste(traits, collapse = ","),
          n_traits = length(traits),
          status = "failed",
          error_message = paste(
            if (inherits(odd_res, "error")) conditionMessage(odd_res) else "",
            if (inherits(all_res, "error")) conditionMessage(all_res) else ""
          )
        )
        next
      }

      odd_fit <- odd_res$fit
      all_fit <- all_res$fit

      results[[ix]] <- data.table(
        panel_id = sprintf("panel_%03d", ix),
        traits = paste(traits, collapse = ","),
        n_traits = length(traits),
        status = "ok",
        odd_cfi = fitMeasures(odd_fit, "cfi"),
        odd_srmr = fitMeasures(odd_fit, "srmr"),
        odd_rmsea = fitMeasures(odd_fit, "rmsea"),
        all_cfi = fitMeasures(all_fit, "cfi"),
        all_srmr = fitMeasures(all_fit, "srmr"),
        all_rmsea = fitMeasures(all_fit, "rmsea"),
        odd_chisq = fitMeasures(odd_fit, "chisq"),
        all_chisq = fitMeasures(all_fit, "chisq"),
        odd_df = fitMeasures(odd_fit, "df"),
        all_df = fitMeasures(all_fit, "df"),
        score = fitMeasures(all_fit, "cfi") + fitMeasures(odd_fit, "cfi") -
          fitMeasures(all_fit, "srmr") - fitMeasures(odd_fit, "srmr") -
          fitMeasures(all_fit, "rmsea") - fitMeasures(odd_fit, "rmsea")
      )
    }
  }
}

results_dt <- rbindlist(results, fill = TRUE)
setorder(results_dt, -score, -all_cfi, all_srmr, all_rmsea, -odd_cfi, odd_srmr, odd_rmsea)
fwrite(results_dt, file.path(output_dir, "nonlipid_ultrapure3_panel_search.tsv"), sep = "\t")

top_ok <- results_dt[status == "ok"][1:min(30, .N)]
fwrite(top_ok, file.path(output_dir, "nonlipid_ultrapure3_panel_search_top30.tsv"), sep = "\t")

message("Finished nonlipid ultrapure3 panel search. Results written to: ", output_dir)
