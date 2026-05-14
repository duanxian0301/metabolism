options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

project_root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
work_dir <- file.path(project_root, "work", "03_ldsc_metabolic_factors_vs_ndd")
out_dir <- file.path(project_root, "results", "03_ldsc_metabolic_factors_vs_ndd")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

clean_dir <- file.path(project_root, "work", "clean_factor_inputs")
hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

trait_inputs <- data.table(
  trait = c(
    "lipid8_F1", "lipid8_F2", "lipid8_F3",
    "nonlipid8_F1", "nonlipid8_F2", "nonlipid8_F3",
    "AD", "PD", "LBD"
  ),
  source_file = c(
    file.path(clean_dir, "lipid8_F1_clean.txt"),
    file.path(clean_dir, "lipid8_F2_clean.txt"),
    file.path(clean_dir, "lipid8_F3_clean.txt"),
    file.path(clean_dir, "nonlipid8_F1_clean.txt"),
    file.path(clean_dir, "nonlipid8_F2_clean.txt"),
    file.path(clean_dir, "nonlipid8_F3_clean.txt"),
    "D:/文章/4NDD/NDDGWAS/AD.txt",
    "D:/文章/4NDD/NDDGWAS/PD.txt",
    "D:/文章/4NDD/NDDGWAS/LBD.txt"
  ),
  N = c(
    276215, 532584, 236386,
    101905, 123638, 34657,
    NA_real_, NA_real_, NA_real_
  ),
  N_source = c(
    rep("factor_neff_mean", 6),
    rep("file_N_column", 3)
  )
)

if (!all(file.exists(trait_inputs$source_file))) {
  print(trait_inputs[!file.exists(source_file)])
  stop("Missing LDSC input file(s).")
}

trait_inputs[, work_txt := file.path(work_dir, paste0(trait, ".txt"))]
trait_inputs[, sumstats := file.path(work_dir, paste0(trait, ".sumstats.gz"))]

for (i in seq_len(nrow(trait_inputs))) {
  file.copy(trait_inputs$source_file[i], trait_inputs$work_txt[i], overwrite = TRUE)
}

need_munge <- !file.exists(trait_inputs$sumstats)
if (any(need_munge)) {
  munge(
    files = trait_inputs$work_txt[need_munge],
    hm3 = hm3,
    trait.names = trait_inputs$trait[need_munge],
    N = trait_inputs$N[need_munge]
  )
}

ldsc_out <- ldsc(
  traits = trait_inputs$sumstats,
  sample.prev = rep(NA, nrow(trait_inputs)),
  population.prev = rep(NA, nrow(trait_inputs)),
  ld = ld,
  wld = wld,
  trait.names = trait_inputs$trait,
  ldsc.log = file.path(out_dir, "metabolic_factors_vs_ndd_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "metabolic_factors_vs_ndd_ldsc.rds"))

S <- as.matrix(ldsc_out$S)
rownames(S) <- colnames(S)
I_mat <- as.matrix(ldsc_out$I)
rownames(I_mat) <- colnames(S)
colnames(I_mat) <- colnames(S)

ratio <- tcrossprod(1 / sqrt(diag(S)))
R <- S * ratio
rownames(R) <- rownames(S)
colnames(R) <- colnames(S)

se_cov <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S))
se_cov[lower.tri(se_cov, diag = TRUE)] <- sqrt(diag(ldsc_out$V))
se_cov[upper.tri(se_cov)] <- t(se_cov)[upper.tri(se_cov)]

ratio_vec <- ratio[lower.tri(ratio, diag = TRUE)]
V_std <- ldsc_out$V * tcrossprod(ratio_vec)

se_rg <- matrix(0, nrow(R), ncol(R), dimnames = dimnames(R))
se_rg[lower.tri(se_rg, diag = TRUE)] <- sqrt(diag(V_std))
se_rg[upper.tri(se_rg)] <- t(se_rg)[upper.tri(se_rg)]

diag_dt <- data.table(
  trait = colnames(S),
  h2 = diag(S),
  h2_se = diag(se_cov),
  h2_z = diag(S) / diag(se_cov),
  intercept = diag(ldsc_out$I)
)

factor_traits <- trait_inputs$trait[1:6]
ndd_traits <- c("AD", "PD", "LBD")

pair_grid <- as.data.table(expand.grid(
  trait1 = factor_traits,
  trait2 = ndd_traits,
  stringsAsFactors = FALSE
))

pair_grid[, `:=`(
  covariance = mapply(function(x, y) S[x, y], trait1, trait2),
  covariance_se = mapply(function(x, y) se_cov[x, y], trait1, trait2),
  rg = mapply(function(x, y) R[x, y], trait1, trait2),
  rg_se = mapply(function(x, y) se_rg[x, y], trait1, trait2),
  intercept = mapply(function(x, y) I_mat[x, y], trait1, trait2)
)]

pair_grid[, z_cov := covariance / covariance_se]
pair_grid[, p_cov := 2 * pnorm(abs(z_cov), lower.tail = FALSE)]
pair_grid[, z_rg := rg / rg_se]
pair_grid[, p_rg := 2 * pnorm(abs(z_rg), lower.tail = FALSE)]
pair_grid[, fdr_rg := p.adjust(p_rg, method = "BH")]
pair_grid[, fdr_cov := p.adjust(p_cov, method = "BH")]
setorder(pair_grid, p_rg)

input_manifest <- trait_inputs[, .(trait, source_file, copied_txt = work_txt, sumstats, N, N_source)]

fwrite(input_manifest, file.path(out_dir, "input_manifest.tsv"), sep = "\t")
fwrite(diag_dt, file.path(out_dir, "metabolic_factors_vs_ndd_ldsc_h2.tsv"), sep = "\t")
fwrite(as.data.table(R, keep.rownames = "trait"), file.path(out_dir, "metabolic_factors_vs_ndd_ldsc_rg_matrix.tsv"), sep = "\t")
fwrite(as.data.table(S, keep.rownames = "trait"), file.path(out_dir, "metabolic_factors_vs_ndd_ldsc_cov_matrix.tsv"), sep = "\t")
fwrite(pair_grid, file.path(out_dir, "metabolic_factors_vs_ndd_requested_pairs.tsv"), sep = "\t")

summary_lines <- c(
  "# Metabolic factors vs AD/PD/LBD LDSC summary",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Top requested factor-disease pairs by p_rg:",
  capture.output(print(pair_grid[, .(trait1, trait2, rg, rg_se, p_rg, fdr_rg, covariance, covariance_se, p_cov, fdr_cov, intercept)][1:min(.N, 18)]))
)
writeLines(summary_lines, file.path(out_dir, "LDSC_summary.md"))

message("Wrote LDSC outputs to: ", out_dir)
