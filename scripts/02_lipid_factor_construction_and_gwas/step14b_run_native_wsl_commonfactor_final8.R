options(stringsAsFactors = FALSE)
.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))

library(GenomicSEM)
library(data.table)

root_dir <- "/mnt/d/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
input_root <- file.path(root_dir, "step14_native_wsl_factor_gwas_inputs")
input_dir <- file.path(input_root, "input")
manifest_file <- file.path(input_root, "chunk_manifest.tsv")
ldsc_file <- file.path(root_dir, "step12_ldsc_lipid_module_final8", "lipid_module_final8_multivariate_ldsc.rds")
output_dir <- file.path(root_dir, "step14_native_wsl_commonfactor_results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

factor_defs <- list(
  F1_TG_rich_axis = c("M_HDL_TG", "VLDL_size", "MUFA", "S_VLDL_TG"),
  F2_HDL_core_axis = c("ApoA1", "HDL_CE"),
  F3_CE_structural_axis = c("XS_VLDL_FC", "VLDL_CE")
)

for (fac in names(factor_defs)) {
  dir.create(file.path(output_dir, fac), showWarnings = FALSE, recursive = TRUE)
}

manifest <- fread(manifest_file)
LDSCoutput <- readRDS(ldsc_file)

trait_order <- colnames(LDSCoutput$S)
trait_keep <- unique(unlist(factor_defs, use.names = FALSE))

if (!all(trait_keep %in% trait_order)) {
  stop("One or more factor traits are missing from LDSC trait order.")
}

start_chunk <- as.integer(Sys.getenv("START_CHUNK", unset = "1"))
end_chunk <- as.integer(Sys.getenv("END_CHUNK", unset = as.character(nrow(manifest))))
use_cores <- as.integer(Sys.getenv("USE_CORES", unset = "4"))

start_chunk <- max(1L, start_chunk)
end_chunk <- min(nrow(manifest), end_chunk)
use_cores <- max(1L, use_cores)

subset_covstruc <- function(covstruc, keep_traits) {
  idx <- match(keep_traits, colnames(covstruc$S))
  full_pairs <- which(lower.tri(covstruc$S, diag = TRUE), arr.ind = TRUE)
  keep_pairs <- which(
    full_pairs[, 1] %in% idx & full_pairs[, 2] %in% idx
  )

  list(
    V = covstruc$V[keep_pairs, keep_pairs, drop = FALSE],
    S = covstruc$S[idx, idx, drop = FALSE],
    I = covstruc$I[idx, idx, drop = FALSE],
    N = covstruc$N[, keep_pairs, drop = FALSE],
    m = covstruc$m
  )
}

for (i in seq.int(start_chunk, end_chunk)) {
  chunk_file <- file.path(input_dir, manifest$file_name[i])
  if (!file.exists(chunk_file)) {
    message("Skipping missing chunk: ", chunk_file)
    next
  }

  message("Running chunk ", i, " / ", nrow(manifest), ": ", manifest$file_name[i])
  snp_chunk <- fread(chunk_file)

  for (fac in names(factor_defs)) {
    traits <- factor_defs[[fac]]
    cov_sub <- subset_covstruc(LDSCoutput, traits)

    beta_cols <- paste0("beta.", traits)
    se_cols <- paste0("se.", traits)
    missing_cols <- setdiff(c(beta_cols, se_cols), names(snp_chunk))
    if (length(missing_cols)) {
      stop("Missing SNP chunk columns for ", fac, ": ", paste(missing_cols, collapse = ", "))
    }

    snp_sub <- snp_chunk[, c("SNP", "A1", "A2", "MAF", beta_cols, se_cols), with = FALSE]
    out_file <- file.path(output_dir, fac, paste0(fac, "_", i, ".tsv"))
    if (file.exists(out_file)) {
      message("Skipping completed factor chunk: ", fac, " chunk ", i)
      next
    }

    result <- tryCatch(
      commonfactorGWAS(
        covstruc = cov_sub,
        SNPs = snp_sub,
        estimation = "DWLS",
        cores = use_cores,
        parallel = use_cores > 1,
        GC = "none",
        MPI = FALSE,
        smooth_check = TRUE
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      message("Factor ", fac, " failed on chunk ", i, ": ", conditionMessage(result))
      next
    }

    fwrite(as.data.frame(result), out_file, sep = "\t")
  }
}

message("Native WSL commonfactorGWAS loop finished.")
