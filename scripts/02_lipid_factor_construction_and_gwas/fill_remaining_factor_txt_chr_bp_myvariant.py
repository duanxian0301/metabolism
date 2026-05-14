from pathlib import Path
import time

import pandas as pd
import requests


TARGETS = [
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F1_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F2_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F3_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F1_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F2_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F3_standard.txt"),
]

CACHE = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\factor_txt_missing_chr_bp_myvariant_cache.tsv")
SUMMARY = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\factor_standard_txt_chr_bp_finalfill_summary.tsv")
CHUNK_SIZE = 1000
SLEEP_SEC = 0.2
URL = "https://myvariant.info/v1/query"


def load_missing_snps():
    dt = pd.read_csv(TARGETS[0], sep="\t", usecols=["SNP", "CHR", "BP"])
    miss = dt[dt["CHR"].isna() | dt["BP"].isna()]["SNP"].drop_duplicates().tolist()
    return miss


def query_chunk(chunk, session):
    for attempt in range(4):
        try:
            response = session.post(
                URL,
                data={
                    "q": ",".join(chunk),
                    "scopes": "dbsnp.rsid",
                    "fields": "dbsnp.rsid,dbsnp.chrom,dbsnp.hg19.start,dbsnp.hg19.end",
                    "species": "human",
                },
                timeout=120,
            )
            response.raise_for_status()
            records = response.json()
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))
    rows = []
    seen = set()
    for rec in records:
        query = rec.get("query")
        if not query or query in seen:
            continue
        dbsnp = rec.get("dbsnp", {})
        hg19 = dbsnp.get("hg19", {})
        chrom = dbsnp.get("chrom")
        start = hg19.get("start")
        if chrom is None or start is None:
            continue
        try:
            chr_int = int(chrom)
            bp_int = int(start)
        except Exception:
            continue
        rows.append({"SNP": query, "CHR": chr_int, "BP": bp_int})
        seen.add(query)
    return rows


def build_lookup(missing):
    existing = None
    if CACHE.exists():
        existing = pd.read_csv(CACHE, sep="\t")
        existing = existing.drop_duplicates("SNP")
    else:
        existing = pd.DataFrame(columns=["SNP", "CHR", "BP"])

    done = set(existing["SNP"].astype(str))
    pending = [s for s in missing if s not in done]
    if not pending:
        return existing

    session = requests.Session()
    collected = [existing]
    total = len(pending)
    for i in range(0, total, CHUNK_SIZE):
        chunk = pending[i : i + CHUNK_SIZE]
        rows = query_chunk(chunk, session)
        chunk_df = pd.DataFrame(rows)
        if not chunk_df.empty:
            collected.append(chunk_df)
            merged = pd.concat(collected, ignore_index=True).drop_duplicates("SNP", keep="first")
            merged.to_csv(CACHE, sep="\t", index=False)
        print(f"queried {min(i + CHUNK_SIZE, total)}/{total}; found {len(rows)} in chunk", flush=True)
        time.sleep(SLEEP_SEC)
    return pd.concat(collected, ignore_index=True).drop_duplicates("SNP", keep="first")


def apply_lookup(lookup):
    summary_rows = []
    for path in TARGETS:
        dt = pd.read_csv(path, sep="\t")
        miss_before = ((dt["CHR"].isna()) | (dt["BP"].isna())).sum()
        fill_map = lookup.rename(columns={"CHR": "CHR_fill", "BP": "BP_fill"})
        dt = dt.merge(fill_map, on="SNP", how="left")
        dt["CHR"] = dt["CHR"].fillna(dt["CHR_fill"])
        dt["BP"] = dt["BP"].fillna(dt["BP_fill"])
        dt = dt.drop(columns=["CHR_fill", "BP_fill"])
        miss_after = ((dt["CHR"].isna()) | (dt["BP"].isna())).sum()
        dt.to_csv(path, sep="\t", index=False)
        summary_rows.append(
            {
                "file": str(path),
                "n_rows": len(dt),
                "n_missing_before": int(miss_before),
                "n_missing_after": int(miss_after),
                "n_newly_filled": int(miss_before - miss_after),
                "pct_filled_total": 1 - (miss_after / len(dt)),
                "lookup_cache": str(CACHE),
            }
        )
    pd.DataFrame(summary_rows).to_csv(SUMMARY, sep="\t", index=False)


def main():
    missing = load_missing_snps()
    lookup = build_lookup(missing)
    apply_lookup(lookup)
    print(SUMMARY)


if __name__ == "__main__":
    main()
