#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 45_prepare_pwcoco_region_inputs_pair.R <trait1> <trait2> <trait2_gwas_file>")
}

trait1 <- args[[1]]
trait2 <- args[[2]]
trait2_file <- args[[3]]

root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
pair <- paste0(trait1, "_", trait2)
region_file <- file.path(root, sprintf("results/08_coloc_pwcoco_loci_%s/%s_regions_500kb.tsv", pair, pair))
out_root <- file.path(root, sprintf("results/10_pwcoco_%s", pair))
input_root <- file.path(out_root, "inputs")
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)

ref_dir <- "D:/mixer/mixer/reference/1000G_EUR_Phase3_plink"
trait1_file <- file.path(root, sprintf("work/clean_factor_inputs/%s_clean.txt", trait1))

regions <- fread(region_file)
trait1_dt <- fread(trait1_file)
trait2_dt <- fread(trait2_file)

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

trait1_dt <- prep(trait1_dt)
trait2_dt <- prep(trait2_dt)

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

  trait1_region <- trait1_dt[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end & SNP %in% ref_region$SNP]
  trait2_region <- trait2_dt[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end & SNP %in% ref_region$SNP]

  trait1_out <- format_pwcoco(harmonize_to_ref(trait1_region, ref_region, flip_freq = TRUE), with_case = FALSE)
  trait2_out <- format_pwcoco(harmonize_to_ref(trait2_region, ref_region, flip_freq = TRUE), with_case = TRUE, n_case = 56156)

  trait1_path <- file.path(region_dir, sprintf("%s_%s_pwcoco.tsv", trait1, r$region_id))
  trait2_path <- file.path(region_dir, sprintf("%s_%s_pwcoco.tsv", trait2, r$region_id))
  fwrite(trait1_out, trait1_path, sep = "\t", quote = FALSE, na = "NA")
  fwrite(trait2_out, trait2_path, sep = "\t", quote = FALSE, na = "NA")

  manifest[[length(manifest) + 1L]] <- data.table(
    region_id = r$region_id,
    chrnum = r$chrnum,
    region_start = r$region_start,
    region_end = r$region_end,
    sentinel_snps = r$sentinel_snps,
    bfile = ref_prefix,
    sum_stats1 = trait1_path,
    sum_stats2 = trait2_path,
    trait1 = trait1,
    trait2 = trait2,
    nsnps_trait1 = nrow(trait1_out),
    nsnps_trait2 = nrow(trait2_out)
  )
}

manifest_dt <- rbindlist(manifest, fill = TRUE)
fwrite(manifest_dt, file.path(out_root, "pwcoco_region_input_manifest.tsv"), sep = "\t")
cat("Wrote PWCoCo inputs to ", out_root, "\n", sep = "")
