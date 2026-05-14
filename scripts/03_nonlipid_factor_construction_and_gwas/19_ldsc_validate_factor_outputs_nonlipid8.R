library(GenomicSEM)
library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
std_dir <- file.path(
  root_dir,
  "step22_native_wsl_usergwas_nonlipid8_results",
  "merged_nonlipid_final8",
  "standard_txt"
)
ascii_work_dir <- "D:/codex/GenomicSEM/metabolic/step23_nonlipid_final8_ldsc_work"
out_dir <- file.path(root_dir, "step23_ldsc_validation_nonlipid_final8")

dir.create(ascii_work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(ascii_work_dir)

source_files <- c(
  file.path(std_dir, "nonlipid_final8_F1_standard.txt"),
  file.path(std_dir, "nonlipid_final8_F2_standard.txt"),
  file.path(std_dir, "nonlipid_final8_F3_standard.txt")
)
trait_names <- c("nonlipid8_F1", "nonlipid8_F2", "nonlipid8_F3")
N <- c(599249, 599249, 599249)

if (!all(file.exists(source_files))) {
  stop("Missing standard factor files in: ", std_dir)
}

files <- file.path(ascii_work_dir, paste0(trait_names, ".txt"))
for (i in seq_along(source_files)) {
  file.copy(source_files[i], files[i], overwrite = TRUE)
}

hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

munge(
  files = files,
  hm3 = hm3,
  trait.names = trait_names,
  N = N
)

traits <- file.path(ascii_work_dir, paste0(trait_names, ".sumstats.gz"))

ldsc_out <- ldsc(
  traits = traits,
  sample.prev = rep(NA, length(traits)),
  population.prev = rep(NA, length(traits)),
  ld = ld,
  wld = wld,
  trait.names = trait_names,
  ldsc.log = file.path(out_dir, "nonlipid_final8_factor_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "nonlipid_final8_factor_ldsc.rds"))

summary_dt <- data.table(
  trait = colnames(ldsc_out$S),
  h2 = diag(ldsc_out$S),
  intercept = diag(ldsc_out$I)
)

rg_mat <- cov2cor(ldsc_out$S)
rg_dt <- data.table(trait = rownames(rg_mat))
rg_dt <- cbind(rg_dt, as.data.table(rg_mat))

fwrite(summary_dt, file.path(out_dir, "nonlipid_final8_factor_ldsc_summary.tsv"), sep = "\t")
fwrite(rg_dt, file.path(out_dir, "nonlipid_final8_factor_rg_matrix.tsv"), sep = "\t")
