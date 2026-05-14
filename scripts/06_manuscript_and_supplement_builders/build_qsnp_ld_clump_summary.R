library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
plink_bin <- "D:/SMR/g1000/plink.exe"
bfile <- "D:/SMR/g1000/g1000_eur"
out_dir <- file.path(root_dir, "q_snp_ld_clump_final8")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

configs <- list(
  list(module = "lipid", factor = "F1", path = file.path(root_dir, "step14_native_wsl_usergwas_final8_results", "merged_lipid_final8", "lipid_final8_F1_userGWAS_merged.tsv.gz")),
  list(module = "lipid", factor = "F2", path = file.path(root_dir, "step14_native_wsl_usergwas_final8_results", "merged_lipid_final8", "lipid_final8_F2_userGWAS_merged.tsv.gz")),
  list(module = "lipid", factor = "F3", path = file.path(root_dir, "step14_native_wsl_usergwas_final8_results", "merged_lipid_final8", "lipid_final8_F3_userGWAS_merged.tsv.gz")),
  list(module = "nonlipid", factor = "F1", path = file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results", "merged_nonlipid_final8", "nonlipid_final8_F1_userGWAS_merged.tsv.gz")),
  list(module = "nonlipid", factor = "F2", path = file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results", "merged_nonlipid_final8", "nonlipid_final8_F2_userGWAS_merged.tsv.gz")),
  list(module = "nonlipid", factor = "F3", path = file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results", "merged_nonlipid_final8", "nonlipid_final8_F3_userGWAS_merged.tsv.gz"))
)

clump_p <- 5e-8
clump_r2 <- 0.1
clump_kb <- 1000L

stopifnot(file.exists(plink_bin))
stopifnot(file.exists(paste0(bfile, ".bed")))
stopifnot(file.exists(paste0(bfile, ".bim")))
stopifnot(file.exists(paste0(bfile, ".fam")))

run_plink <- function(args) {
  result <- system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE)
  invisible(result)
}

read_clumped <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.table(SNP = character(), CHR = integer(), BP = integer(), P = numeric()))
  }
  dt <- fread(path, fill = TRUE)
  if (!"SNP" %in% names(dt)) {
    return(data.table(SNP = character(), CHR = integer(), BP = integer(), P = numeric()))
  }
  dt
}

make_assoc <- function(dt, pcol, outfile) {
  assoc <- dt[, .(SNP, P = get(pcol))]
  assoc <- assoc[!is.na(P) & is.finite(P) & P > 0 & P <= 1]
  fwrite(assoc, outfile, sep = "\t")
  assoc
}

clump_one <- function(assoc_file, out_prefix) {
  args <- c(
    "--bfile", bfile,
    "--clump", assoc_file,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", format(clump_p, scientific = TRUE),
    "--clump-p2", format(clump_p, scientific = TRUE),
    "--clump-r2", as.character(clump_r2),
    "--clump-kb", as.character(clump_kb),
    "--allow-extra-chr",
    "--threads", "8",
    "--out", out_prefix
  )
  run_plink(args)
  read_clumped(paste0(out_prefix, ".clumped"))
}

ld_overlap_detail <- function(factor_leads, qsnp_leads, prefix) {
  if (nrow(factor_leads) == 0 || nrow(qsnp_leads) == 0) {
    return(list(detail = data.table(), unique_factor = factor_leads[, .(SNP, CHR, BP, factor_P = P)]))
  }
  factor_list <- file.path(out_dir, paste0(prefix, "_factor_lead_snps.txt"))
  qsnp_list <- file.path(out_dir, paste0(prefix, "_qsnp_lead_snps.txt"))
  union_list <- file.path(out_dir, paste0(prefix, "_union_snps.txt"))
  fwrite(factor_leads[, .(SNP)], factor_list, col.names = FALSE, sep = "\t")
  fwrite(qsnp_leads[, .(SNP)], qsnp_list, col.names = FALSE, sep = "\t")
  fwrite(data.table(SNP = unique(c(factor_leads$SNP, qsnp_leads$SNP))), union_list, col.names = FALSE, sep = "\t")

  ld_prefix <- file.path(out_dir, paste0(prefix, "_factor_vs_qsnp_ld"))
  args <- c(
    "--bfile", bfile,
    "--extract", union_list,
    "--ld-snp-list", factor_list,
    "--r2", "yes-really",
    "--ld-window", "999999",
    "--ld-window-kb", as.character(clump_kb),
    "--ld-window-r2", as.character(clump_r2),
    "--allow-extra-chr",
    "--threads", "8",
    "--out", ld_prefix
  )
  run_plink(args)
  ld_file <- paste0(ld_prefix, ".ld")
  if (!file.exists(ld_file) || file.info(ld_file)$size == 0) {
    return(list(detail = data.table(), unique_factor = factor_leads[, .(SNP, CHR, BP, factor_P = P)]))
  }
  ld_dt <- fread(ld_file)
  if (!all(c("SNP_A", "SNP_B", "R2") %in% names(ld_dt))) {
    return(list(detail = data.table(), unique_factor = factor_leads[, .(SNP, CHR, BP, factor_P = P)]))
  }
  ld_dt <- ld_dt[SNP_B %in% qsnp_leads$SNP]
  if (nrow(ld_dt) == 0) {
    return(list(detail = data.table(), unique_factor = factor_leads[, .(SNP, CHR, BP, factor_P = P)]))
  }
  detail <- merge(
    ld_dt,
    factor_leads[, .(SNP_A = SNP, factor_CHR = CHR, factor_BP = BP, factor_P = P)],
    by = "SNP_A",
    all.x = TRUE
  )
  detail <- merge(
    detail,
    qsnp_leads[, .(SNP_B = SNP, qsnp_CHR = CHR, qsnp_BP = BP, qsnp_P = P)],
    by = "SNP_B",
    all.x = TRUE
  )
  unique_factor <- factor_leads[!SNP %in% unique(detail$SNP_A), .(SNP, CHR, BP, factor_P = P)]
  list(detail = detail, unique_factor = unique_factor)
}

summary_rows <- list()
factor_lead_rows <- list()
qsnp_lead_rows <- list()
ld_overlap_rows <- list()
unique_factor_rows <- list()

for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  message("Processing ", cfg$module, " ", cfg$factor)
  dt <- fread(cfg$path, select = c("SNP", "CHR", "BP", "Pval_Estimate", "Q_SNP_pval"))
  dt <- dt[!is.na(SNP) & !is.na(CHR) & !is.na(BP)]

  prefix <- paste(cfg$module, cfg$factor, sep = "_")
  factor_assoc_file <- file.path(out_dir, paste0(prefix, "_factor_assoc.tsv"))
  qsnp_assoc_file <- file.path(out_dir, paste0(prefix, "_qsnp_assoc.tsv"))

  factor_assoc <- make_assoc(dt, "Pval_Estimate", factor_assoc_file)
  qsnp_assoc <- make_assoc(dt, "Q_SNP_pval", qsnp_assoc_file)

  factor_clumped <- clump_one(factor_assoc_file, file.path(out_dir, paste0(prefix, "_factor")))
  qsnp_clumped <- clump_one(qsnp_assoc_file, file.path(out_dir, paste0(prefix, "_qsnp")))

  factor_leads <- merge(
    factor_clumped[, .(SNP, CHR, BP, P)],
    dt[, .(SNP, factor_P = Pval_Estimate)],
    by = "SNP",
    all.x = TRUE
  )
  qsnp_leads <- merge(
    qsnp_clumped[, .(SNP, CHR, BP, P)],
    dt[, .(SNP, qsnp_P = Q_SNP_pval)],
    by = "SNP",
    all.x = TRUE
  )

  factor_leads[, `:=`(module = cfg$module, factor = cfg$factor)]
  qsnp_leads[, `:=`(module = cfg$module, factor = cfg$factor)]
  factor_lead_rows[[length(factor_lead_rows) + 1L]] <- factor_leads
  qsnp_lead_rows[[length(qsnp_lead_rows) + 1L]] <- qsnp_leads

  overlap_obj <- ld_overlap_detail(factor_leads, qsnp_leads, prefix)
  overlap_dt <- overlap_obj$detail
  if (nrow(overlap_dt)) {
    overlap_dt[, `:=`(module = cfg$module, factor = cfg$factor)]
    ld_overlap_rows[[length(ld_overlap_rows) + 1L]] <- overlap_dt
  }
  unique_dt <- overlap_obj$unique_factor
  if (nrow(unique_dt)) {
    unique_dt[, `:=`(module = cfg$module, factor = cfg$factor)]
    unique_factor_rows[[length(unique_factor_rows) + 1L]] <- unique_dt
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    module = cfg$module,
    factor = cfg$factor,
    clump_reference = "g1000_eur",
    clump_p_threshold = clump_p,
    clump_r2 = clump_r2,
    clump_kb = clump_kb,
    n_factor_gws = sum(dt$Pval_Estimate < clump_p, na.rm = TRUE),
    n_qsnp_gws = sum(dt$Q_SNP_pval < clump_p, na.rm = TRUE),
    n_factor_lead_loci_ldclump = nrow(factor_leads),
    n_qsnp_lead_loci_ldclump = nrow(qsnp_leads),
    n_exact_overlap_lead_snps = sum(factor_leads$SNP %in% qsnp_leads$SNP),
    n_factor_leads_excluding_ld_with_qsnp = nrow(unique_dt)
  )
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
factor_leads_dt <- rbindlist(factor_lead_rows, use.names = TRUE, fill = TRUE)
qsnp_leads_dt <- rbindlist(qsnp_lead_rows, use.names = TRUE, fill = TRUE)
ld_overlap_dt <- rbindlist(ld_overlap_rows, use.names = TRUE, fill = TRUE)
unique_factor_dt <- rbindlist(unique_factor_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(out_dir, "final8_qsnp_ldclump_summary.tsv"), sep = "\t")
fwrite(factor_leads_dt, file.path(out_dir, "final8_factor_lead_loci_ldclump.tsv"), sep = "\t")
fwrite(qsnp_leads_dt, file.path(out_dir, "final8_qsnp_lead_loci_ldclump.tsv"), sep = "\t")
fwrite(ld_overlap_dt, file.path(out_dir, "final8_factor_qsnp_ld_overlap_pairs.tsv"), sep = "\t")
fwrite(unique_factor_dt, file.path(out_dir, "final8_unique_factor_leads_excluding_qsnp_ld.tsv"), sep = "\t")

manifest_dt <- data.table(
  file = c(
    "final8_qsnp_ldclump_summary.tsv",
    "final8_factor_lead_loci_ldclump.tsv",
    "final8_qsnp_lead_loci_ldclump.tsv",
    "final8_factor_qsnp_ld_overlap_pairs.tsv",
    "final8_unique_factor_leads_excluding_qsnp_ld.tsv"
  ),
  description = c(
    "Module-factor summary of genome-wide significant SNP counts, LD-clumped lead loci, exact lead-SNP overlap, and unique factor lead loci after excluding LD with Q_SNP lead loci.",
    "LD-clumped lead loci for factor GWAS results.",
    "LD-clumped lead loci for Q_SNP results.",
    "Pairwise LD overlap records between factor lead SNPs and Q_SNP lead SNPs at r2 >= 0.1 within 1000 kb.",
    "Factor lead loci remaining after excluding SNPs in LD with Q_SNP lead loci."
  )
)
fwrite(manifest_dt, file.path(out_dir, "README_qsnp_ldclump_manifest.tsv"), sep = "\t")

print(summary_dt)
