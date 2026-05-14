# Metabolic factor-NDD genomic workflow

This repository contains the analysis code used for the manuscript:

**Metabolic Genetic Axes Linked to Alzheimer and Parkinson Disease Through Distinct Polygenic, Cellular, and Target-Annotation Programs**

The code is organized in manuscript workflow order and intentionally excludes raw GWAS summary statistics, intermediate work directories, generated figures, supplementary tables, Word files, and temporary/debug scripts.

## Repository structure

- `scripts/01_metabolic_trait_qc_and_selection/`: metabolic-trait QC, univariate LDSC entry filtering, and lipid/non-lipid candidate panel selection.
- `scripts/02_lipid_factor_construction_and_gwas/`: lipid factor construction, Genomic SEM factor GWAS generation, factor QC, and internal validation.
- `scripts/03_nonlipid_factor_construction_and_gwas/`: non-lipid factor construction, Genomic SEM factor GWAS generation, factor QC, and internal validation.
- `scripts/04_ndd_overlap_locus_gene_cell_target_workflow/`: factor-NDD LDSC, MiXeR, pleiotropy-informed shared loci, coloc/PWCoCo, cTWAS/SMR, scPagwas, scTenifoldKnk, target annotation, and final supplementary table assembly.
- `scripts/05_figure_generation/`: final figure-generation scripts for manuscript figures.
- `scripts/06_manuscript_and_supplement_builders/`: scripts used to assemble final manuscript, supplementary methods/information, and supplementary tables from processed outputs.
- `docs/pipeline_manifest.csv`: ordered mapping from release file to original local file.

## Data availability

This code release does not redistribute third-party GWAS summary statistics, single-cell data, eQTL reference panels, LDSC reference files, MiXeR reference files, or generated manuscript result files. Users should download these resources from the original providers and update local path variables inside scripts before execution.

## Software dependencies

The workflow uses R, Python, WSL/Linux command-line tools, LDSC, GenomicSEM, MiXeR, coloc, scPagwas, scTenifoldKnk, SMR, and standard scientific Python/R packages. A non-exhaustive package list is provided in `requirements/`.

## Notes on paths

Many scripts preserve the project-local path variables used for the manuscript analysis. Before rerunning the workflow on another machine, update root directories and external resource paths at the top of each script.
