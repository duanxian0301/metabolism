library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
supp_dir <- file.path(root_dir, "supplement_nonlipid_final8")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

manifest_file <- file.path(root_dir, "final_model_nonlipid_module_final8", "final8_trait_manifest.tsv")
esem_summary_file <- file.path(root_dir, "step20_efa_esem_nonlipid_module_ultrapure8", "nonlipid_module_ultrapure8_esem_summary.tsv")
loadings_file <- file.path(root_dir, "step20_efa_esem_nonlipid_module_ultrapure8", "ALL_3factor_loadings.tsv")
std_meta_file <- file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results", "merged_nonlipid_final8", "standard_txt", "nonlipid_final8_standard_txt_metadata.tsv")
ldsc_log_file <- file.path(root_dir, "step23_ldsc_validation_nonlipid_final8", "nonlipid_final8_factor_ldsc_ldsc.log")
ldsc_summary_file <- file.path(root_dir, "step23_ldsc_validation_nonlipid_final8", "nonlipid_final8_factor_ldsc_summary.tsv")
rg_matrix_file <- file.path(root_dir, "step23_ldsc_validation_nonlipid_final8", "nonlipid_final8_factor_rg_matrix.tsv")

manifest_dt <- fread(manifest_file)
esem_summary_dt <- fread(esem_summary_file)
loadings_dt <- fread(loadings_file)
std_meta_dt <- fread(std_meta_file)
ldsc_summary_dt <- fread(ldsc_summary_file)
rg_matrix_dt <- fread(rg_matrix_file)

fwrite(manifest_dt, file.path(supp_dir, "Supplementary_Table_S1_nonlipid_final8_trait_manifest.tsv"), sep = "\t")
fwrite(esem_summary_dt, file.path(supp_dir, "Supplementary_Table_S2_nonlipid_final8_esem_fit.tsv"), sep = "\t")
fwrite(loadings_dt, file.path(supp_dir, "Supplementary_Table_S3_nonlipid_final8_all3factor_loadings.tsv"), sep = "\t")
fwrite(std_meta_dt, file.path(supp_dir, "Supplementary_Table_S4_nonlipid_final8_factor_gwas_file_manifest.tsv"), sep = "\t")

log_lines <- readLines(ldsc_log_file, warn = FALSE)

trim_path <- function(x) sub("^.*[/\\\\]([^/\\\\]+)\\.sumstats\\.gz$", "\\1", x)

parse_univariate <- function(lines) {
  idx <- grep("^Heritability Results for trait:", lines)
  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    block <- lines[idx[i]:min(length(lines), idx[i] + 8L)]
    trait_line <- block[1]
    trait <- trim_path(sub("^Heritability Results for trait:\\s*", "", trait_line))
    get_pair <- function(pattern) {
      hit <- grep(pattern, block, value = TRUE)
      if (!length(hit)) return(c(NA_real_, NA_real_))
      m <- regexec(pattern, hit[1], perl = TRUE)
      reg <- regmatches(hit[1], m)[[1]]
      as.numeric(reg[2:3])
    }
    get_single <- function(pattern) {
      hit <- grep(pattern, block, value = TRUE)
      if (!length(hit)) return(NA_real_)
      m <- regexec(pattern, hit[1], perl = TRUE)
      reg <- regmatches(hit[1], m)[[1]]
      as.numeric(reg[2])
    }
    mean_chi2 <- get_single("^Mean Chi\\^2 across remaining SNPs:\\s*(-?[0-9.]+)$")
    lambda_gc <- get_single("^Lambda GC:\\s*(-?[0-9.]+)$")
    intercept <- get_pair("^Intercept:\\s*(-?[0-9.]+) \\((-?[0-9.]+)\\)$")
    ratio <- get_pair("^Ratio:\\s*(-?[0-9.]+) \\((-?[0-9.]+)\\)$")
    h2 <- get_pair("^Total Observed Scale h2:\\s*(-?[0-9.]+) \\((-?[0-9.]+)\\)$")
    h2_z <- get_single("^h2 Z:\\s*(-?[0-9.]+)$")
    out[[i]] <- data.table(
      trait = trait,
      mean_chi2 = mean_chi2,
      lambda_gc = lambda_gc,
      intercept = intercept[1],
      intercept_se = intercept[2],
      ratio = ratio[1],
      ratio_se = ratio[2],
      h2 = h2[1],
      h2_se = h2[2],
      h2_z = h2_z
    )
  }
  rbindlist(out, fill = TRUE)
}

parse_bivariate <- function(lines) {
  idx <- grep("^Results for genetic covariance between:", lines)
  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    start <- idx[i]
    block <- lines[start:min(length(lines), start + 8L)]
    pair_txt <- sub("^Results for genetic covariance between:\\s*", "", block[1])
    traits <- strsplit(pair_txt, " and ", fixed = TRUE)[[1]]
    trait1 <- trim_path(traits[1])
    trait2 <- trim_path(traits[2])
    get_pair <- function(pattern) {
      hit <- grep(pattern, block, value = TRUE)
      if (!length(hit)) return(c(NA_real_, NA_real_))
      m <- regexec(pattern, hit[1], perl = TRUE)
      reg <- regmatches(hit[1], m)[[1]]
      as.numeric(reg[2:3])
    }
    get_single <- function(pattern) {
      hit <- grep(pattern, block, value = TRUE)
      if (!length(hit)) return(NA_real_)
      m <- regexec(pattern, hit[1], perl = TRUE)
      reg <- regmatches(hit[1], m)[[1]]
      as.numeric(reg[2])
    }
    mean_zz <- get_single("^Mean Z\\*Z:\\s*(-?[0-9.]+)$")
    cross_intercept <- get_pair("^Cross trait Intercept:\\s*(-?[0-9.]+) \\((-?[0-9.]+)\\)$")
    gcov <- get_pair("^Total Observed Scale Genetic Covariance \\(g_cov\\):\\s*(-?[0-9.]+) \\((-?[0-9.]+)\\)$")
    gcov_z <- get_single("^g_cov Z:\\s*(-?[0-9.]+)$")
    gcov_p <- get_single("^g_cov P-value:\\s*([-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?)$")
    out[[i]] <- data.table(
      trait1 = trait1,
      trait2 = trait2,
      mean_z_z = mean_zz,
      cross_trait_intercept = cross_intercept[1],
      cross_trait_intercept_se = cross_intercept[2],
      g_cov = gcov[1],
      g_cov_se = gcov[2],
      g_cov_z = gcov_z,
      g_cov_p = gcov_p
    )
  }
  rbindlist(out, fill = TRUE)
}

parse_rg <- function(lines) {
  idx <- grep("^Genetic Correlation between ", lines)
  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    hit <- lines[idx[i]]
    m <- regexec("^Genetic Correlation between ([^ ]+) and ([^:]+):\\s*(-?[0-9.]+) \\(([0-9.]+)\\)$", hit)
    reg <- regmatches(hit, m)[[1]]
    if (!length(reg)) next
    out[[i]] <- data.table(
      trait1 = reg[2],
      trait2 = reg[3],
      rg = as.numeric(reg[4]),
      rg_se = as.numeric(reg[5])
    )
  }
  rbindlist(out, fill = TRUE)
}

parse_qc <- function(lines) {
  read_idx <- grep("^Read in summary statistics \\[[0-9]+/[0-9]+\\] from:", lines)
  out <- vector("list", length(read_idx))
  for (i in seq_along(read_idx)) {
    start <- read_idx[i]
    block <- lines[start:min(length(lines), start + 3L)]
    trait <- trim_path(sub("^Read in summary statistics \\[[0-9]+/[0-9]+\\] from:\\s*", "", block[1]))
    m1 <- regexec("^Out of ([0-9]+) SNPs, ([0-9]+) remain after merging with LD-score files$", block[2], perl = TRUE)
    reg1 <- regmatches(block[2], m1)[[1]]
    remove_line <- block[3]
    m2 <- regexec("^Removing ([0-9]+) SNPs with Chi\\^2 > [0-9.]+; ([0-9]+) remain$", remove_line, perl = TRUE)
    reg2 <- regmatches(remove_line, m2)[[1]]
    out[[i]] <- data.table(
      trait = trait,
      snps_post_munge = as.numeric(reg1[2]),
      snps_after_ld_merge = as.numeric(reg1[3]),
      snps_removed_chisq = as.numeric(reg2[2]),
      snps_after_chisq = as.numeric(reg2[3])
    )
  }
  rbindlist(out, fill = TRUE)
}

uni_dt <- parse_univariate(log_lines)
biv_dt <- parse_bivariate(log_lines)
rg_dt <- parse_rg(log_lines)
qc_dt <- parse_qc(log_lines)

supp_s5 <- merge(qc_dt, uni_dt, by = "trait", all = TRUE)
supp_s5 <- merge(supp_s5, ldsc_summary_dt, by = "trait", all.x = TRUE, suffixes = c("", "_summary"))
fwrite(supp_s5, file.path(supp_dir, "Supplementary_Table_S5_nonlipid_final8_univariate_ldsc_detailed.tsv"), sep = "\t")

supp_s6 <- merge(biv_dt, rg_dt, by = c("trait1", "trait2"), all = TRUE)
fwrite(supp_s6, file.path(supp_dir, "Supplementary_Table_S6_nonlipid_final8_bivariate_ldsc_detailed.tsv"), sep = "\t")
fwrite(rg_matrix_dt, file.path(supp_dir, "Supplementary_Table_S7_nonlipid_final8_rg_matrix.tsv"), sep = "\t")

overview_lines <- c(
  "Supplementary overview for nonlipid_final8",
  paste0("Final model traits: ", nrow(manifest_dt)),
  paste0("Final factors: ", paste(unique(manifest_dt$final_factor_name), collapse = ", ")),
  paste0("ALL 3-factor ESEM CFI: ", sprintf("%.4f", esem_summary_dt[dataset == "ALL" & nfactors == 3, cfi])),
  paste0("ALL 3-factor ESEM SRMR: ", sprintf("%.4f", esem_summary_dt[dataset == "ALL" & nfactors == 3, srmr])),
  paste0("ALL 3-factor ESEM RMSEA: ", sprintf("%.4f", esem_summary_dt[dataset == "ALL" & nfactors == 3, rmsea])),
  paste0("Factor GWAS outputs merged under: ", file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results", "merged_nonlipid_final8")),
  "Detailed LDSC metrics are in Supplementary_Table_S5 and Supplementary_Table_S6."
)
writeLines(overview_lines, file.path(supp_dir, "nonlipid_final8_supplement_overview.txt"))
