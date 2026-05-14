#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
region_file <- file.path(root, "results/08_coloc_pwcoco_loci_lipid8_F2_AD/lipid8_F2_AD_regions_500kb.tsv")
out_root <- file.path(root, "results/10_pwcoco_lipid8_F2_AD")
input_root <- file.path(out_root, "inputs")
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)

ref_dir <- "D:/mixer/mixer/reference/1000G_EUR_Phase3_plink"
lipid_file <- file.path(root, "work/clean_factor_inputs/lipid8_F2_clean.txt")
ad_file <- "D:/文章/4NDD/NDDGWAS/AD.txt"

regions <- fread(region_file)
lipid <- fread(lipid_file)
ad <- fread(ad_file)

comp_map <- c(A = "T", T = "A", C = "G", G = "C")
comp_vec <- function(x) {
  y <- comp_map[toupper(x)]
  y[is.na(y)] <- NA_character_
  unname(y)
}

prep <- function(dt) {
  dt <- dt[!is.na(SNP) & !is.na(CHR) & !is.na(BP) & !is.na(A1) & !is.na(A2)]
  dt <- dt[!is.na(FREQ) & FREQ > 0 & FREQ < 1 & !is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(P) & P > 0 & P <= 1 & !is.na(N) & N > 0]
  dt[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  unique(dt, by = "SNP")
}

lipid <- prep(lipid)
ad <- prep(ad)

harmonize_to_ref <- function(dt, ref, flip_freq = TRUE) {
  merged <- merge(dt, ref[, .(SNP, ref_A1 = A1, ref_A2 = A2)], by = "SNP")
  if (!nrow(merged)) return(merged[0])
  merged[, `:=`(A1_comp = comp_vec(A1), A2_comp = comp_vec(A2))]
  same <- merged$A1 == merged$ref_A1 & merged$A2 == merged$ref_A2
  swap <- merged$A1 == merged$ref_A2 & merged$A2 == merged$ref_A1
  comp_same <- merged$A1_comp == merged$ref_A1 & merged$A2_comp == merged$ref_A2
  comp_swap <- merged$A1_comp == merged$ref_A2 & merged$A2_comp == merged$ref_A1
  keep <- same | swap | comp_same | comp_swap
  out <- merged[keep]
  if (!nrow(out)) return(out)
  flip <- (swap | comp_swap)[keep]
  out[flip, BETA := -BETA]
  if (flip_freq) out[flip, FREQ := 1 - FREQ]
  out[, `:=`(A1 = ref_A1, A2 = ref_A2)]
  out[!is.na(FREQ) & FREQ > 0 & FREQ < 1]
}

format_pwcoco <- function(dt, with_case = FALSE, n_case = NA_real_) {
  if (!nrow(dt)) return(dt)
  if (with_case) {
    dt[, .(SNP, A1, A2, A1_freq = FREQ, beta = BETA, se = SE, p = P, n = N, case = n_case)]
  } else {
    dt[, .(SNP, A1, A2, A1_freq = FREQ, beta = BETA, se = SE, p = P, n = N)]
  }
}

manifest <- list()
for (i in seq_len(nrow(regions))) {
  r <- regions[i]
  region_dir <- file.path(input_root, r$region_id)
  dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
  ref_prefix <- file.path(ref_dir, sprintf("1000G.EUR.QC.%s", r$chrnum))
  ref_bim <- paste0(ref_prefix, ".bim")
  if (!file.exists(ref_bim)) stop("Missing BIM: ", ref_bim)
  ref <- fread(ref_bim, header = FALSE, showProgress = FALSE)
  setnames(ref, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  ref_region <- ref[BP >= r$region_start & BP <= r$region_end, .(SNP, BP, A1, A2)]
  if (!nrow(ref_region)) next

  lipid_region <- lipid[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end & SNP %in% ref_region$SNP]
  ad_region <- ad[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end & SNP %in% ref_region$SNP]

  lipid_out <- format_pwcoco(harmonize_to_ref(lipid_region, ref_region, flip_freq = TRUE), with_case = FALSE)
  ad_out <- format_pwcoco(harmonize_to_ref(ad_region, ref_region, flip_freq = TRUE), with_case = TRUE, n_case = 85934)

  lipid_path <- file.path(region_dir, sprintf("lipid8_F2_%s_pwcoco.tsv", r$region_id))
  ad_path <- file.path(region_dir, sprintf("AD_%s_pwcoco.tsv", r$region_id))
  fwrite(lipid_out, lipid_path, sep = "\t", quote = FALSE, na = "NA")
  fwrite(ad_out, ad_path, sep = "\t", quote = FALSE, na = "NA")

  manifest[[length(manifest) + 1L]] <- data.table(
    region_id = r$region_id,
    chrnum = r$chrnum,
    region_start = r$region_start,
    region_end = r$region_end,
    sentinel_snps = r$sentinel_snps,
    bfile = ref_prefix,
    sum_stats1 = lipid_path,
    sum_stats2 = ad_path,
    nsnps_lipid8_F2 = nrow(lipid_out),
    nsnps_AD = nrow(ad_out)
  )
}

manifest_dt <- rbindlist(manifest, fill = TRUE)
fwrite(manifest_dt, file.path(out_root, "pwcoco_region_input_manifest.tsv"), sep = "\t")
cat("Wrote PWCoCo inputs to ", out_root, "\n", sep = "")
