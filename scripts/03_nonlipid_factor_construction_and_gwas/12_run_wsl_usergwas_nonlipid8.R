options(stringsAsFactors = FALSE)
.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))

library(GenomicSEM)
library(data.table)

root_dir <- "/mnt/d/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
input_root <- file.path(root_dir, "step21_native_wsl_factor_gwas_inputs_nonlipid8")
input_dir <- file.path(input_root, "input")
manifest_file <- file.path(input_root, "chunk_manifest.tsv")
ldsc_file <- file.path(root_dir, "step19_ldsc_nonlipid_module_ultrapure8", "nonlipid_module_ultrapure8_multivariate_ldsc.rds")
output_dir <- file.path(root_dir, "step22_native_wsl_usergwas_nonlipid8_results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "F1"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "F2"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "F3"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Q_SNP"), showWarnings = FALSE, recursive = TRUE)

manifest <- fread(manifest_file)
LDSCoutput <- readRDS(ldsc_file)

start_chunk <- as.integer(Sys.getenv("START_CHUNK", unset = "1"))
end_chunk <- as.integer(Sys.getenv("END_CHUNK", unset = as.character(nrow(manifest))))
use_cores <- as.integer(Sys.getenv("USE_CORES", unset = "8"))

model <- "
F1_ketone_axis =~ Acetoacetate + bOHbutyrate
F2_amino_acid_axis =~ Val + Leu + Phe
F3_energy_bridge_axis =~ Acetate + Lactate + Glucose
F1_ketone_axis ~~ F2_amino_acid_axis
F1_ketone_axis ~~ F3_energy_bridge_axis
F2_amino_acid_axis ~~ F3_energy_bridge_axis
F1_ketone_axis ~ SNP
F2_amino_acid_axis ~ SNP
F3_energy_bridge_axis ~ SNP
"

writeLines(model, file.path(output_dir, "nonlipid8_3factor_userGWAS_model.txt"))

sub_targets <- c(
  "F1_ketone_axis~SNP",
  "F2_amino_acid_axis~SNP",
  "F3_energy_bridge_axis~SNP"
)

for (i in seq.int(start_chunk, end_chunk)) {
  chunk_file <- file.path(input_dir, manifest$file_name[i])
  if (!file.exists(chunk_file)) next

  f1_out <- file.path(output_dir, "F1", paste0("F1_", i, ".tsv"))
  f2_out <- file.path(output_dir, "F2", paste0("F2_", i, ".tsv"))
  f3_out <- file.path(output_dir, "F3", paste0("F3_", i, ".tsv"))
  q_out <- file.path(output_dir, "Q_SNP", paste0("Q_SNP_", i, ".tsv"))

  if (file.exists(f1_out) && file.exists(f2_out) && file.exists(f3_out)) next

  snp_chunk <- fread(chunk_file)
  result <- tryCatch(
    userGWAS(
      covstruc = LDSCoutput,
      SNPs = snp_chunk,
      estimation = "DWLS",
      model = model,
      printwarn = TRUE,
      sub = sub_targets,
      cores = use_cores,
      parallel = use_cores > 1,
      GC = "none",
      MPI = FALSE,
      smooth_check = TRUE,
      fix_measurement = TRUE,
      std.lv = TRUE,
      Q_SNP = TRUE
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) next

  fwrite(as.data.frame(result[[1]]), f1_out, sep = "\t")
  fwrite(as.data.frame(result[[2]]), f2_out, sep = "\t")
  fwrite(as.data.frame(result[[3]]), f3_out, sep = "\t")
  if (length(result) >= 4 && !is.null(result[[4]])) {
    fwrite(as.data.frame(result[[4]]), q_out, sep = "\t")
  }
}
