from __future__ import annotations

import argparse
from pathlib import Path


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
SCRIPTS = ROOT / "scripts"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait1", required=True)
    parser.add_argument("--trait2", required=True)
    parser.add_argument("--trait2-label", required=True)
    args = parser.parse_args()

    pair = f"{args.trait1}_{args.trait2}"
    out_path = SCRIPTS / f"39_config_pleiofdr_{pair}.txt"
    traitfiles_cell = "{'" + f"{args.trait2_label}_fdr.mat" + "'}"
    traitnames_cell = "{'" + f"{args.trait2}" + "'}"

    text = f"""# Metabolic {args.trait1} x {args.trait2} conjFDR config.\n# Copy this file to D:\\pleioFDR\\pleiofdr-master\\config.txt before running runme.m.\n\nreffile=ref9545380_1kgPhase3eur_LDr2p1.mat\ntraitfolder=data\n\ntraitfile1={args.trait1}_fdr.mat\ntraitname1={args.trait1}\ntraitfiles={traitfiles_cell}\ntraitnames={traitnames_cell}\n\noutputdir=output_{pair}\n\nstattype=conjfdr\nfdrthresh=0.05\n\nrandprune=true\nrandprune_n=500\n\nexclude_chr_pos = [6 25119106 33854733; 8 7200000 12500000]\n\nmanh_fontsize_genenames=12\nmanh_yspace=0.75\nmanh_ymargin=0.25\nmanh_colorlist=[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]\n\nrefinfo=9545380.ref\n\nreset_pruneidx=true\nrandprune_repeats=default\npthresh=1\nperform_gc=true\nuse_standard_gc=false\nrandprune_gc=true\nexclude_from_discovery=false\nmafthresh = 0.005\nexclude_ambiguous_snps = true\nonscreen = false\ndummy_zscore = false\nexit_matlab_upon_completion = true\n"""
    out_path.write_text(text, encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
