#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 44_run_coloc_pair.R <trait1> <trait2> <trait2_gwas_file>")
}

trait1 <- args[[1]]
trait2 <- args[[2]]
trait2_file <- args[[3]]

root <- "D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd"
pair <- paste0(trait1, "_", trait2)
region_file <- file.path(root, sprintf("results/08_coloc_pwcoco_loci_%s/%s_regions_500kb.tsv", pair, pair))
out_dir <- file.path(root, sprintf("results/09_coloc_%s", pair))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trait1_file <- file.path(root, sprintf("work/clean_factor_inputs/%s_clean.txt", trait1))

regions <- fread(region_file)
trait1_dt <- fread(trait1_file)
trait2_dt <- fread(trait2_file)

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

trait1_dt <- prep(trait1_dt)
trait2_dt <- prep(trait2_dt)

comp_map <- c(A = "T", T = "A", C = "G", G = "C")
comp <- function(x) unname(comp_map[toupper(x)])

harmonize_pair <- function(x, y) {
  m <- merge(
    x[, .(SNP, CHR, BP, A1_1 = A1, A2_1 = A2, BETA_1 = BETA, SE_1 = SE, P_1 = P, N_1 = N, MAF_1 = MAF)],
    y[, .(SNP, A1_2 = A1, A2_2 = A2, BETA_2 = BETA, SE_2 = SE, P_2 = P, N_2 = N, MAF_2 = MAF)],
    by = "SNP"
  )
  if (!nrow(m)) return(m)
  m[, `:=`(A1_2c = comp(A1_2), A2_2c = comp(A2_2))]
  same <- (m$A1_1 == m$A1_2 & m$A2_1 == m$A2_2) | (m$A1_1 == m$A1_2c & m$A2_1 == m$A2_2c)
  swap <- (m$A1_1 == m$A2_2 & m$A2_1 == m$A1_2) | (m$A1_1 == m$A2_2c & m$A2_1 == m$A1_2c)
  m <- m[same | swap]
  if (!nrow(m)) return(m)
  m[swap[same | swap], BETA_2 := -BETA_2]
  m
}

results <- list()
idx <- 0L

for (i in seq_len(nrow(regions))) {
  r <- regions[i]
  x <- trait1_dt[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end]
  y <- trait2_dt[CHR == r$chrnum & BP >= r$region_start & BP <= r$region_end]
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
    beta = m$BETA_1,
    varbeta = m$SE_1^2,
    N = stats::median(m$N_1, na.rm = TRUE),
    MAF = m$MAF_1,
    type = "quant",
    sdY = 1
  )
  d2 <- list(
    snp = m$SNP,
    beta = m$BETA_2,
    varbeta = m$SE_2^2,
    N = stats::median(m$N_2, na.rm = TRUE),
    MAF = m$MAF_2,
    type = "cc",
    s = 0.33
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
fwrite(out, file.path(out_dir, sprintf("coloc_%s_regions.tsv", pair)), sep = "\t")
fwrite(out[, .N, by = coloc_class][order(coloc_class)], file.path(out_dir, sprintf("coloc_%s_class_counts.tsv", pair)), sep = "\t")
cat("Wrote coloc results to ", out_dir, "\n", sep = "")
