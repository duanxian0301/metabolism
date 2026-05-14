args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 55_export_ctwas_pair_summary.R <pair> <factor_trait> <disease_trait>")
}

pair <- args[[1]]
factor_trait <- args[[2]]
disease_trait <- args[[3]]

base_linux <- file.path("/home/shenjing/ctwas_paths", paste0("metabolic_mctwas_work_", pair), "single_trait")
out_dir <- file.path("/mnt/d/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/results", paste0("16_ctwas_", pair), "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_obj <- function(path) {
  obj <- readRDS(path)
  as.data.frame(obj, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

export_trait <- function(trait) {
  res_dir <- file.path(base_linux, trait, "results")
  if (!dir.exists(res_dir)) {
    stop(sprintf("Missing result dir: %s", res_dir))
  }

  z_gene <- read_obj(file.path(res_dir, sprintf("%s.z_gene.RDS", trait)))
  finemap <- read_obj(file.path(res_dir, sprintf("%s.finemap_res.RDS", trait)))
  susie_alpha <- read_obj(file.path(res_dir, sprintf("%s.susie_alpha_res.RDS", trait)))
  boundary <- read_obj(file.path(res_dir, sprintf("%s.boundary_genes.RDS", trait)))
  run_summary <- read.delim(file.path(res_dir, sprintf("%s_run_summary.tsv", trait)), sep = "\t", stringsAsFactors = FALSE)

  write_tsv(z_gene, file.path(out_dir, sprintf("%s_z_gene.tsv", trait)))
  write_tsv(finemap, file.path(out_dir, sprintf("%s_finemap_res.tsv", trait)))
  write_tsv(susie_alpha, file.path(out_dir, sprintf("%s_susie_alpha_res.tsv", trait)))
  write_tsv(boundary, file.path(out_dir, sprintf("%s_boundary_genes.tsv", trait)))

  if ("z" %in% names(z_gene)) {
    z_gene$abs_z <- abs(as.numeric(z_gene$z))
    z_top50 <- z_gene[order(-z_gene$abs_z), , drop = FALSE]
    z_top50 <- head(z_top50, 50)
    z_top50$abs_z <- NULL
    write_tsv(z_top50, file.path(out_dir, sprintf("%s_top50_z_gene.tsv", trait)))
  }

  if ("susie_pip" %in% names(finemap)) {
    finemap$susie_pip <- suppressWarnings(as.numeric(finemap$susie_pip))
    fm_top50 <- finemap[order(-finemap$susie_pip), , drop = FALSE]
    fm_top50 <- head(fm_top50, 50)
    write_tsv(fm_top50, file.path(out_dir, sprintf("%s_top50_finemap.tsv", trait)))
  }

  counts <- data.frame(
    trait = trait,
    n_z_gene = nrow(z_gene),
    n_finemap = nrow(finemap),
    n_susie_alpha = nrow(susie_alpha),
    n_boundary_genes = nrow(boundary),
    n_weights = if ("n_weights" %in% names(run_summary)) run_summary$n_weights[[1]] else NA,
    n_harmonized = if ("n_harmonized" %in% names(run_summary)) run_summary$n_harmonized[[1]] else NA,
    run_time = if ("run_time" %in% names(run_summary)) run_summary$run_time[[1]] else NA,
    stringsAsFactors = FALSE
  )
  write_tsv(counts, file.path(out_dir, sprintf("%s_counts.tsv", trait)))
}

export_trait(factor_trait)
export_trait(disease_trait)

message(sprintf("Wrote cTWAS summaries to %s", out_dir))
