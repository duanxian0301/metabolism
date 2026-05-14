library(GenomicSEM)
library(data.table)
library(lavaan)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
gwas_dir <- "D:/metabolic/GWAS"
model_dir <- file.path(root_dir, "step13_efa_esem_lipid_module_final8")
ldsc_dir <- file.path(root_dir, "step12_ldsc_lipid_module_final8")
panel_path <- file.path(root_dir, "final_model_lipid_module_final8", "final8_trait_manifest.tsv")
output_dir <- file.path(root_dir, "step14_usergwas_lipid_module_final8")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[[1]] else "test"
test_n <- if (length(args) >= 2) as.integer(args[[2]]) else 200L
test_n <- ifelse(is.na(test_n) || test_n <= 0, 200L, test_n)

panel_dt <- fread(panel_path)
fit_all <- readRDS(file.path(model_dir, "ALL_3factor_fit.rds"))
covstruc <- readRDS(file.path(ldsc_dir, "lipid_module_final8_multivariate_ldsc.rds"))
pt <- as.data.table(parTable(fit_all))

format_num <- function(x, digits = 10) {
  formatC(as.numeric(x), digits = digits, format = "fg", flag = "#")
}

build_fixed_model_syntax <- function(pt, traits) {
  loading_dt <- pt[op == "=~", .(lhs, rhs, est)]
  resid_dt <- pt[op == "~~" & lhs == rhs & lhs %in% traits, .(lhs, est)]
  latent_var_dt <- pt[op == "~~" & lhs == rhs & lhs %in% unique(loading_dt$lhs), .(lhs, est)]
  latent_cov_dt <- pt[op == "~~" & lhs != rhs & lhs %in% unique(loading_dt$lhs) & rhs %in% unique(loading_dt$lhs), .(lhs, rhs, est)]

  factor_levels <- unique(loading_dt$lhs)
  lines <- character()

  for (fac in factor_levels) {
    sub_dt <- loading_dt[lhs == fac]
    rhs_terms <- paste0(format_num(sub_dt$est), "*", sub_dt$rhs)
    lines <- c(lines, paste0(fac, " =~ ", paste(rhs_terms, collapse = " + ")))
  }

  for (fac in factor_levels) {
    est_val <- latent_var_dt[lhs == fac, est][1]
    lines <- c(lines, paste0(fac, " ~~ ", format_num(est_val), "*", fac))
  }

  for (i in seq_len(nrow(latent_cov_dt))) {
    lines <- c(
      lines,
      paste0(
        latent_cov_dt$lhs[i],
        " ~~ ",
        format_num(latent_cov_dt$est[i]),
        "*",
        latent_cov_dt$rhs[i]
      )
    )
  }

  for (i in seq_len(nrow(resid_dt))) {
    lines <- c(lines, paste0(resid_dt$lhs[i], " ~~ ", format_num(resid_dt$est[i]), "*", resid_dt$lhs[i]))
  }

  for (fac in factor_levels) {
    lines <- c(lines, paste0(fac, " ~ SNP"))
  }

  for (trait in traits) {
    lines <- c(lines, paste0(trait, " ~ 0*SNP"))
  }

  paste(lines, collapse = "\n")
}

build_marker_model_syntax <- function(panel_dt) {
  factor_levels <- unique(panel_dt$final_factor_name)
  lines <- character()

  for (fac in factor_levels) {
    sub_traits <- panel_dt[final_factor_name == fac, trait_code]
    lines <- c(lines, paste0(fac, " =~ ", paste(sub_traits, collapse = " + ")))
  }

  if (length(factor_levels) >= 2) {
    for (i in seq_len(length(factor_levels) - 1)) {
      for (j in (i + 1):length(factor_levels)) {
        lines <- c(lines, paste0(factor_levels[[i]], " ~~ ", factor_levels[[j]]))
      }
    }
  }

  for (fac in factor_levels) {
    lines <- c(lines, paste0(fac, " ~ SNP"))
  }

  for (trait in panel_dt$trait_code) {
    lines <- c(lines, paste0(trait, " ~ 0*SNP"))
  }

  paste(lines, collapse = "\n")
}

read_one_gwas <- function(accession, trait_code, n_max = NULL) {
  path <- file.path(gwas_dir, paste0(accession, ".txt"))
  dt <- fread(
    path,
    select = c("SNP", "A1", "A2", "MAF", "BETA", "SE", "N"),
    nrows = n_max,
    showProgress = FALSE
  )
  setnames(
    dt,
    old = c("A1", "A2", "MAF", "BETA", "SE", "N"),
    new = c(
      paste0("A1__", trait_code),
      paste0("A2__", trait_code),
      paste0("MAF__", trait_code),
      paste0("beta.", trait_code),
      paste0("se.", trait_code),
      paste0("N__", trait_code)
    )
  )
  dt
}

merge_gwas_panel <- function(panel_dt, n_max = NULL) {
  trait_codes <- panel_dt$trait_code
  anchor_trait <- trait_codes[[1]]

  merged_dt <- read_one_gwas(panel_dt$study_accession[[1]], anchor_trait, n_max = n_max)
  setnames(
    merged_dt,
    old = c(paste0("A1__", anchor_trait), paste0("A2__", anchor_trait), paste0("MAF__", anchor_trait), paste0("N__", anchor_trait)),
    new = c("A1", "A2", "MAF_anchor", "N_anchor")
  )

  for (i in 2:nrow(panel_dt)) {
    accession <- panel_dt$study_accession[[i]]
    trait_code <- panel_dt$trait_code[[i]]
    dt_i <- read_one_gwas(accession, trait_code, n_max = n_max)
    merged_dt <- merge(merged_dt, dt_i, by = "SNP", all = FALSE)

    a1_col <- paste0("A1__", trait_code)
    a2_col <- paste0("A2__", trait_code)
    beta_col <- paste0("beta.", trait_code)

    same_idx <- merged_dt[[a1_col]] == merged_dt[["A1"]] & merged_dt[[a2_col]] == merged_dt[["A2"]]
    flip_idx <- merged_dt[[a1_col]] == merged_dt[["A2"]] & merged_dt[[a2_col]] == merged_dt[["A1"]]
    keep_idx <- same_idx | flip_idx
    merged_dt <- merged_dt[keep_idx]

    if (any(flip_idx[keep_idx])) {
      merged_dt[flip_idx[keep_idx], (beta_col) := -get(beta_col)]
    }
  }

  maf_cols <- grep("^MAF__", names(merged_dt), value = TRUE)
  merged_dt[, MAF := apply(.SD, 1, median, na.rm = TRUE), .SDcols = maf_cols]

  out_cols <- c("SNP", "A1", "A2", "MAF")
  for (trait_code in trait_codes) {
    out_cols <- c(out_cols, paste0("beta.", trait_code), paste0("se.", trait_code))
  }

  out_dt <- merged_dt[, ..out_cols]
  out_dt
}

fixed_model_syntax <- build_fixed_model_syntax(pt, panel_dt$trait_code)
marker_model_syntax <- build_marker_model_syntax(panel_dt)
writeLines(fixed_model_syntax, file.path(output_dir, "final8_fixed_measurement_model.txt"))
writeLines(marker_model_syntax, file.path(output_dir, "final8_usergwas_marker_model.txt"))
fwrite(pt, file.path(output_dir, "ALL_3factor_parTable.tsv"), sep = "\t")

sub_targets <- c("f1~SNP", "f2~SNP", "f3~SNP")

if (identical(run_mode, "test")) {
  message("Building test SNP input from first ", test_n, " rows of each raw GWAS file...")
  test_dt <- merge_gwas_panel(panel_dt, n_max = test_n)
  saveRDS(test_dt, file.path(output_dir, "test_input_subset.rds"))
  fwrite(test_dt, file.path(output_dir, "final8_usergwas_test_input.tsv"), sep = "\t")

  diag_dt <- data.table(
    metric = c("n_traits", "n_snps", "run_mode", "test_n"),
    value = c(nrow(panel_dt), nrow(test_dt), run_mode, test_n)
  )
  fwrite(diag_dt, file.path(output_dir, "input_diagnostics.tsv"), sep = "\t")

  message("Running userGWAS test on first ", nrow(test_dt), " SNPs...")
  test_res <- userGWAS(
    covstruc = covstruc,
    SNPs = test_dt,
    estimation = "DWLS",
    model = marker_model_syntax,
    sub = sub_targets,
    parallel = FALSE,
    std.lv = TRUE,
    fix_measurement = TRUE,
    Q_SNP = TRUE
  )

  for (i in seq_along(test_res)) {
    out_name <- paste0("usergwas_test_", gsub("~", "_to_", sub_targets[[i]]), ".tsv")
    fwrite(as.data.table(test_res[[i]]), file.path(output_dir, out_name), sep = "\t")
  }

  message("Test run finished.")
} else if (identical(run_mode, "prep_full")) {
  message("Building full merged SNP input from raw GWAS txt files...")
  snp_dt <- merge_gwas_panel(panel_dt)
  fwrite(snp_dt, file.path(output_dir, "final8_usergwas_input.tsv.gz"), sep = "\t")

  diag_dt <- data.table(
    metric = c("n_traits", "n_snps", "run_mode", "test_n"),
    value = c(nrow(panel_dt), nrow(snp_dt), run_mode, test_n)
  )
  fwrite(diag_dt, file.path(output_dir, "input_diagnostics.tsv"), sep = "\t")

  message("Full merged SNP input written.")
} else if (identical(run_mode, "full")) {
  message("Building full merged SNP input from raw GWAS txt files...")
  snp_dt <- merge_gwas_panel(panel_dt)
  fwrite(snp_dt, file.path(output_dir, "final8_usergwas_input.tsv.gz"), sep = "\t")

  diag_dt <- data.table(
    metric = c("n_traits", "n_snps", "run_mode", "test_n"),
    value = c(nrow(panel_dt), nrow(snp_dt), run_mode, test_n)
  )
  fwrite(diag_dt, file.path(output_dir, "input_diagnostics.tsv"), sep = "\t")

  message("Running full userGWAS on all SNPs...")
  full_res <- userGWAS(
    covstruc = covstruc,
    SNPs = snp_dt,
    estimation = "DWLS",
    model = marker_model_syntax,
    sub = sub_targets,
    parallel = TRUE,
    cores = max(1L, parallel::detectCores() - 1L),
    std.lv = TRUE,
    fix_measurement = TRUE,
    Q_SNP = TRUE
  )

  for (i in seq_along(full_res)) {
    out_name <- paste0("usergwas_full_", gsub("~", "_to_", sub_targets[[i]]), ".tsv.gz")
    fwrite(as.data.table(full_res[[i]]), file.path(output_dir, out_name), sep = "\t")
  }

  message("Full run finished.")
} else {
  stop("Unknown run_mode: ", run_mode, ". Use 'test', 'prep_full', or 'full'.")
}

message("Results written to: ", output_dir)
