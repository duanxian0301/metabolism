options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

project_root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
work_dir <- file.path(project_root, "work", "75_ldsc_factor_vs_component_ndd")
out_dir <- file.path(project_root, "results", "28_ldsc_factor_vs_component_ndd")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

# NOTE:
# Manuscript-facing factor labels follow the downstream project convention already
# used in the AD/PD results writing. In particular, lipid8_F2 is treated as the
# TG/VLDL-rich axis because that is how the final manuscript-facing factor was
# interpreted and used throughout the downstream analyses.
selected_traits <- data.table(
  panel = c(
    rep("lipid8_F2_vs_AD", 3),
    rep("nonlipid8_F1_vs_PD", 3),
    rep("lipid8_F1_vs_PD", 5)
  ),
  trait = c(
    "lipid8_F2", "ApoA1", "HDL_CE",
    "nonlipid8_F1", "Acetoacetate", "bOHbutyrate",
    "lipid8_F1", "VLDL_size", "M_HDL_TG", "S_VLDL_TG", "MUFA"
  ),
  role = c(
    "factor", rep("component", 2),
    "factor", rep("component", 2),
    "factor", rep("component", 4)
  ),
  target_disease = c(
    rep("AD", 3),
    rep("PD", 3),
    rep("PD", 5)
  ),
  rationale = c(
    "Formal GenomicSEM factor for AD line",
    "Direct indicator of HDL-core factor",
    "Direct indicator of HDL-core factor",
    "Manuscript main factor for PD ketone line",
    "Direct indicator of ketone-axis factor",
    "Direct indicator of ketone-axis factor",
    "Formal GenomicSEM factor for PD TG/VLDL-rich line",
    "Direct indicator of TG/VLDL-rich factor",
    "Direct indicator of TG/VLDL-rich factor",
    "Direct indicator of TG/VLDL-rich factor",
    "Direct indicator of TG/VLDL-rich factor"
  )
)

factor_inputs <- data.table(
  trait = c(
    "lipid8_F1", "lipid8_F2", "lipid8_F3",
    "nonlipid8_F1", "nonlipid8_F2", "nonlipid8_F3"
  ),
  source_file = c(
    file.path(project_root, "work", "clean_factor_inputs", "lipid8_F1_clean.txt"),
    file.path(project_root, "work", "clean_factor_inputs", "lipid8_F2_clean.txt"),
    file.path(project_root, "work", "clean_factor_inputs", "lipid8_F3_clean.txt"),
    file.path(project_root, "work", "clean_factor_inputs", "nonlipid8_F1_clean.txt"),
    file.path(project_root, "work", "clean_factor_inputs", "nonlipid8_F2_clean.txt"),
    file.path(project_root, "work", "clean_factor_inputs", "nonlipid8_F3_clean.txt")
  ),
  N = c(276215, 532584, 236386, 101905, 123638, 34657),
  source_type = "factor_clean_txt"
)

lipid_manifest <- fread("D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/final_model_lipid_module_final8/final8_trait_manifest.tsv")
nonlipid_manifest <- fread("D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/final_model_nonlipid_module_final8/final8_trait_manifest.tsv")
component_manifest <- rbindlist(list(lipid_manifest, nonlipid_manifest), fill = TRUE)
component_manifest <- component_manifest[, .(
  trait = trait_code,
  source_file = sumstats_file,
  N = as.numeric(NA),
  source_type = "component_sumstats_gz",
  biomarker_name,
  group,
  final_factor_name
)]

disease_inputs <- data.table(
  trait = c("AD", "PD"),
  source_file = c(
    file.path(project_root, "work", "03_ldsc_metabolic_factors_vs_ndd", "AD.txt"),
    file.path(project_root, "work", "03_ldsc_metabolic_factors_vs_ndd", "PD.txt")
  ),
  N = c(as.numeric(NA), as.numeric(NA)),
  source_type = "ndd_sumstats_txt"
)

all_inputs <- unique(rbindlist(list(
  factor_inputs[, .(trait, source_file, N, source_type)],
  component_manifest[, .(trait, source_file, N, source_type)],
  disease_inputs
), fill = TRUE), by = "trait")

trait_inputs <- merge(
  unique(selected_traits[, .(trait)]),
  all_inputs,
  by = "trait",
  all.x = TRUE
)
trait_inputs <- unique(rbindlist(list(
  trait_inputs,
  disease_inputs[trait %in% unique(selected_traits$target_disease)]
), fill = TRUE), by = "trait")

if (!all(file.exists(trait_inputs$source_file))) {
  print(trait_inputs[!file.exists(source_file)])
  stop("Missing LDSC input file(s) for factor-vs-component comparison.")
}

trait_inputs[, copied_input := source_file]
trait_inputs[, sumstats := file.path(work_dir, paste0(trait, ".sumstats.gz"))]

component_rows <- trait_inputs[source_type == "component_sumstats_gz"]
if (nrow(component_rows) > 0) {
  for (i in seq_len(nrow(component_rows))) {
    raw_path <- component_rows$source_file[i]
    out_path <- file.path(work_dir, paste0(component_rows$trait[i], "_ldsc_input.txt"))
    dt <- fread(raw_path)
    if (!all(c("SNP", "A1", "A2", "Z", "N") %in% names(dt))) {
      stop("Component sumstats missing required columns for trait: ", component_rows$trait[i])
    }
    dt <- dt[, .(SNP, A1, A2, Z, N)]
    dt[, P := 2 * pnorm(abs(Z), lower.tail = FALSE)]
    fwrite(dt, out_path, sep = "\t")
    trait_inputs[trait == component_rows$trait[i], copied_input := out_path]
  }
}

need_munge <- !file.exists(trait_inputs$sumstats)
if (any(need_munge)) {
  munge(
    files = trait_inputs$copied_input[need_munge],
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
  ldsc.log = file.path(out_dir, "factor_vs_component_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "factor_vs_component_ldsc.rds"))

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

pair_grid <- merge(
  selected_traits,
  unique(selected_traits[, .(target_disease)]),
  by = "target_disease",
  allow.cartesian = TRUE
)
pair_grid <- pair_grid[, .(
  panel,
  focal_trait = trait,
  role,
  target_disease,
  rationale
)]

pair_grid[, `:=`(
  covariance = mapply(function(x, y) S[x, y], focal_trait, target_disease),
  covariance_se = mapply(function(x, y) se_cov[x, y], focal_trait, target_disease),
  rg = mapply(function(x, y) R[x, y], focal_trait, target_disease),
  rg_se = mapply(function(x, y) se_rg[x, y], focal_trait, target_disease),
  intercept = mapply(function(x, y) I_mat[x, y], focal_trait, target_disease)
)]
pair_grid[, z_cov := covariance / covariance_se]
pair_grid[, p_cov := 2 * pnorm(abs(z_cov), lower.tail = FALSE)]
pair_grid[, z_rg := rg / rg_se]
pair_grid[, p_rg := 2 * pnorm(abs(z_rg), lower.tail = FALSE)]
pair_grid[, fdr_panel_rg := ave(p_rg, panel, FUN = function(x) p.adjust(x, method = "BH"))]
pair_grid[, fdr_global_rg := p.adjust(p_rg, method = "BH")]
pair_grid[, trait_rank_within_panel := frank(p_rg, ties.method = "min"), by = panel]
setorder(pair_grid, panel, role, p_rg)

summary_rows <- pair_grid[, .(
  factor_trait = focal_trait[role == "factor"][1],
  factor_rg = rg[role == "factor"][1],
  factor_p = p_rg[role == "factor"][1],
  factor_fdr_panel = fdr_panel_rg[role == "factor"][1],
  best_component = focal_trait[role == "component"][which.min(p_rg[role == "component"])],
  best_component_rg = min(rg[role == "component"], na.rm = TRUE),
  best_component_p = min(p_rg[role == "component"], na.rm = TRUE),
  n_components_nominal = sum(role == "component" & p_rg < 0.05, na.rm = TRUE),
  factor_more_significant_than_all_components = p_rg[role == "factor"][1] < min(p_rg[role == "component"], na.rm = TRUE)
), by = .(panel, target_disease)]

fwrite(trait_inputs, file.path(out_dir, "input_manifest.tsv"), sep = "\t")
fwrite(selected_traits, file.path(out_dir, "selected_traits.tsv"), sep = "\t")
fwrite(pair_grid, file.path(out_dir, "factor_vs_component_requested_pairs.tsv"), sep = "\t")
fwrite(summary_rows, file.path(out_dir, "factor_vs_component_summary.tsv"), sep = "\t")

summary_lines <- c(
  "# LDSC comparison: latent factors versus representative component metabolites",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Panels analyzed:",
  "- lipid8_F2 vs AD (components: ApoA1, HDL_CE)",
  "- nonlipid8_F1 vs PD (components: Acetoacetate, bOHbutyrate)",
  "- lipid8_F1 vs PD (components: VLDL_size, M_HDL_TG, S_VLDL_TG, MUFA)",
  "",
  "Panel-level summary:",
  capture.output(print(summary_rows))
)
writeLines(summary_lines, file.path(out_dir, "LDSC_factor_vs_component_summary.md"))

message("Wrote LDSC factor-vs-component outputs to: ", out_dir)
