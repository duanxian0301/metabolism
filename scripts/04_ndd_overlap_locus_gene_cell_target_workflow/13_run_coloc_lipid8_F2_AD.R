#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
region_file <- file.path(root, "results/08_coloc_pwcoco_loci_lipid8_F2_AD/lipid8_F2_AD_regions_500kb.tsv")
out_dir <- file.path(root, "results/09_coloc_lipid8_F2_AD")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lipid_file <- file.path(root, "work/clean_factor_inputs/lipid8_F2_clean.txt")
ad_file <- "D:/文章/4NDD/NDDGWAS/AD.txt"

regions <- fread(region_file)
lipid <- fread(lipid_file)
ad <- fread(ad_file)

prep <- function(dt) {
  dt <- dt[!is.na(SNP) & !is.na(CHR) & !is.na(BP) & !is.na(A1) & !is.na(A2)]
  dt <- dt[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(P) & P > 0 & P <= 1 & !is.na(N) & N > 0]
  dt[, `:=`(
    A1 = toupper(A1),
    A2 = toupper(A2),
    MAF = pmin(FREQ, 1 - FREQ)
  )]
  dt <- dt[!is.na(MAF) & MAF > 0 & MAF < 0.5]
  unique(dt, by = "SNP")
}

lipid <- prep(lipid)
ad <- prep(ad)

comp_map <- c(A = "T", T = "A", C = "G", G = "C")
comp <- function(x) unname(comp_map[toupper(x)])

harmonize_pair <- function(x, y) {
  m <- merge(
    x[, .(SNP, CHR, BP, A1_l = A1, A2_l = A2, BETA_l = BETA, SE_l = SE, P_l = P, N_l = N, MAF_l = MAF)],
    y[, .(SNP, A1_a = A1, A2_a = A2, BETA_a = BETA, SE_a = SE, P_a = P, N_a = N, MAF_a = MAF)],
    by = "SNP"
  )
  if (!nrow(m)) return(m)
  m[, `:=`(A1_ac = comp(A1_a), A2_ac = comp(A2_a))]
  same <- (m$A1_l == m$A1_a & m$A2_l == m$A2_a) | (m$A1_l == m$A1_ac & m$A2_l == m$A2_ac)
  swap <- (m$A1_l == m$A2_a & m$A2_l == m$A1_a) | (m$A1_l == m$A2_ac & m$A2_l == m$A1_ac)
  m <- m[same | swap]
  if (!nrow(m)) return(m)
  m[swap[same | swap], BETA_a := -BETA_a]
  m
}

results <- list()
idx <- 0L

for (i in seq_len(nrow(regions))) {
  r <- regions[i]
  x <- lipid[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end]
  y <- ad[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end]
  m <- harmonize_pair(x, y)
  if (nrow(m) < 50) {
    idx <- idx + 1L
    results[[idx]] <- data.table(
      region_id = r$region_id,
      chrnum = r$chrnum,
      region_start = r$region_start,
      region_end = r$region_end,
      sentinel_snps = r$sentinel_snps,
      nsnps = nrow(m),
      PP.H0 = NA_real_,
      PP.H1 = NA_real_,
      PP.H2 = NA_real_,
      PP.H3 = NA_real_,
      PP.H4 = NA_real_,
      status = "too_few_snps"
    )
    next
  }

  d1 <- list(
    snp = m$SNP,
    beta = m$BETA_l,
    varbeta = m$SE_l^2,
    N = stats::median(m$N_l, na.rm = TRUE),
    MAF = m$MAF_l,
    type = "quant",
    sdY = 1
  )
  d2 <- list(
    snp = m$SNP,
    beta = m$BETA_a,
    varbeta = m$SE_a^2,
    N = stats::median(m$N_a, na.rm = TRUE),
    MAF = m$MAF_a,
    type = "cc",
    s = 0.35
  )

  res <- tryCatch(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), error = function(e) e)
  if (inherits(res, "error")) {
    idx <- idx + 1L
    results[[idx]] <- data.table(
      region_id = r$region_id,
      chrnum = r$chrnum,
      region_start = r$region_start,
      region_end = r$region_end,
      sentinel_snps = r$sentinel_snps,
      nsnps = nrow(m),
      PP.H0 = NA_real_,
      PP.H1 = NA_real_,
      PP.H2 = NA_real_,
      PP.H3 = NA_real_,
      PP.H4 = NA_real_,
      status = paste0("error: ", conditionMessage(res))
    )
  } else {
    s <- res$summary
    idx <- idx + 1L
    results[[idx]] <- data.table(
      region_id = r$region_id,
      chrnum = r$chrnum,
      region_start = r$region_start,
      region_end = r$region_end,
      sentinel_snps = r$sentinel_snps,
      nsnps = unname(s["nsnps"]),
      PP.H0 = unname(s["PP.H0.abf"]),
      PP.H1 = unname(s["PP.H1.abf"]),
      PP.H2 = unname(s["PP.H2.abf"]),
      PP.H3 = unname(s["PP.H3.abf"]),
      PP.H4 = unname(s["PP.H4.abf"]),
      status = "ok"
    )
  }
}

out <- rbindlist(results, fill = TRUE)
setorder(out, -PP.H4, PP.H3)
out[, coloc_class := fifelse(PP.H4 >= 0.8, "strong_H4",
  fifelse(PP.H4 >= 0.5, "moderate_H4",
    fifelse(PP.H3 >= 0.8, "distinct_signals_H3", "inconclusive")
  )
)]
fwrite(out, file.path(out_dir, "coloc_lipid8_F2_AD_regions.tsv"), sep = "\t")
fwrite(out[, .N, by = coloc_class][order(coloc_class)], file.path(out_dir, "coloc_lipid8_F2_AD_class_counts.tsv"), sep = "\t")
cat("Wrote coloc results to ", out_dir, "\n", sep = "")
