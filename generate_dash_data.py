# generate_dash_data.py
import csv
import re
import glob
import os
from pathlib import Path
import pandas as pd
import pycountry
from Bio import SeqIO

# ---------- Color maps ----------
HBV_COLOR_MAP = {
    "0000ff": "HBV-A",
    "996633": "HBV-B",
    "00ffff": "HBV-C",
    "00ff00": "HBV-D",
    "ff00ff": "HBV-E",
    "ff8000": "HBV-F",
    "800080": "HBV-G",
    "ff0000": "HBV-H",
    "ffff00": "HBV-I",
    "800040": "HBV-J",
}
HCV_COLOR_MAP = {
    "0000ff": "HCV-1",
    "996633": "HCV-2",
    "00ffff": "HCV-3",
    "00ff00": "HCV-4",
    "ff00ff": "HCV-5",
    "ff8000": "HCV-6",
    "800080": "HCV-7",
    "ff0000": "HCV-8",
}

ID_RE = re.compile(r"([A-Z]+\d+\.\d+)")  # e.g., OQ948327.1, PV555119.1

# ---------- Helpers ----------
def load_recombinant_ids(file_path: str) -> set:
    p = Path(file_path)
    if not p.is_file():
        return set()
    return {line.strip() for line in p.read_text().splitlines() if line.strip()}

def extract_country_from_location(location: str) -> str:
    if not isinstance(location, str):
        return "Unknown"

    region_overrides = {
        "St. Petersburg": "Russia",
        "Altai Republic": "Russia",
        "Yakutia (Sakha) Republic": "Russia",
        "Ural federal district": "Russia",
        "Tyva Republic": "Russia",
        "Novosibirsk region": "Russia",
        "West Bank": "Palestine",
        "Cape Verde": "Cabo Verde",
        "Democratic Republic of the Congo": "Democratic Republic of the Congo",
        "Cote d'Ivoire": "Ivory Coast",
    }
    for region, country in region_overrides.items():
        if region.lower() in location.lower():
            return country

    parts = re.split(r"[:/]", location)
    parts = [p.strip() for p in parts if p.strip()]
    for part in reversed(parts):
        try:
            country = pycountry.countries.lookup(part)
            return country.name
        except LookupError:
            continue
    return "Unknown"

def extract_color_genotypes(tree_file: str, color_map: dict) -> dict:
    tip_genotypes = {}
    with open(tree_file) as f:
        for line in f:
            line = line.strip()
            m = re.match(r"'([^']+)'\[\&!color=#([0-9a-fA-F]{6})\]", line)
            if m:
                taxon = m.group(1)
                color = m.group(2).lower()
                genotype = color_map.get(color, "Unknown")
                tip_genotypes[taxon] = genotype
    return tip_genotypes

def parse_fasta_ids(fasta_path: str) -> dict:
    ids = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = record.id  # keep accession exactly
        ids[acc] = acc
    return ids

def parse_metadata(metadata_csv: str) -> dict:
    meta = {}
    with open(metadata_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row.get("accession_id")
            location = row.get("location")
            date = row.get("isolate_collection_date") or row.get("release_date")
            country = extract_country_from_location(location)
            meta[acc] = {"Country": country, "Date": date}
    return meta

def generate_summary_tsv(
    fasta_file: str,
    tree_file: str,
    metadata_file: str,
    color_map: dict,
    output_file: str,
    recombs_file: str | None = None,
):
    fasta_ids = parse_fasta_ids(fasta_file)
    metadata = parse_metadata(metadata_file)
    genotype_map = extract_color_genotypes(tree_file, color_map)
    recombinant_ids = load_recombinant_ids(recombs_file) if recombs_file else set()

    with open(output_file, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["Taxa", "Genotype", "Country", "Date"])
        for acc in fasta_ids:
            genotype = "Recombinant" if acc in recombinant_ids else genotype_map.get(acc, "Unknown")
            country = metadata.get(acc, {}).get("Country", "Unknown")
            raw_date = metadata.get(acc, {}).get("Date", "Unknown")

            # Extract YYYY
            if raw_date and raw_date != "Unknown" and "-" in raw_date:
                year = raw_date.split("-")[0]
            else:
                year = raw_date

            writer.writerow([acc, genotype, country, year])

# ---------- Mutation extractors ----------
def extract_hbv_mutations(input_folder: str, output_file: str) -> pd.DataFrame:
    """
    Reads g2p HBV CSVs and writes a tall table:
    columns -> ID, drug, mutation
    """
    pattern = os.path.join(input_folder, "g2p_hbv_data_hbv_*.csv")
    paths = sorted(glob.glob(pattern))
    if not paths:
        print(f"⚠️ No HBV G2P files found at {pattern}")
        pd.DataFrame(columns=["ID", "drug", "mutation"]).to_csv(output_file, index=False)
        return pd.DataFrame(columns=["ID", "drug", "mutation"])

    drug_map = {
        "Lamivudine": ("Lamivudine_class", "Lamivudine_muts"),
        "Adefovir": ("Adefovir_class", "Adefovir_muts"),
        "Entecavir": ("Entecavir_class", "Entecavir_muts"),
        "Tenofovir": ("Tenofovir_class", "Tenofovir_muts"),
        "Telbivudine": ("Telbivudine_class", "Telbivudine_muts"),
    }

    rows = []
    for p in paths:
        df = pd.read_csv(p)
        # extract GenBank-like ID from the 'id' column
        df["ID"] = (
            df.get("id", "")
              .astype(str)
              .apply(lambda s: (m := ID_RE.search(s)) and m.group(1))
        )

        for _, r in df.iterrows():
            rid = r.get("ID")
            if not rid:
                continue
            for drug, (class_col, muts_col) in drug_map.items():
                cls = str(r.get(class_col, "")).strip().lower()
                if cls == "resistant":
                    muts = str(r.get(muts_col, "")).split(",")
                    for mut in muts:
                        mut = mut.strip()
                        if mut and mut.lower() != "none":
                            rows.append({"ID": rid, "drug": drug, "mutation": mut})

    out = pd.DataFrame(rows)
    out.to_csv(output_file, index=False)
    print(f"✅ HBV mutations: wrote {len(out):,} rows to {output_file}")
    return out

def extract_hcv_mutations(input_folder: str, output_file: str) -> pd.DataFrame:
    """
    Reads g2p HCV CSVs and writes a tall table:
    columns -> ID, gene, mutation
    """
    pattern = os.path.join(input_folder, "g2p_hcv_data_hcv_*.csv")
    paths = sorted(glob.glob(pattern))
    if not paths:
        print(f"⚠️ No HCV G2P files found at {pattern}")
        pd.DataFrame(columns=["ID", "gene", "mutation"]).to_csv(output_file, index=False)
        return pd.DataFrame(columns=["ID", "gene", "mutation"])

    resistance_cols = {
        "all_resistance_muts_ns5b": "ns5b",
        "all_resistance_muts_ns5a": "ns5a",
        "all_resistance_muts_ns3": "ns3",
    }

    rows = []
    for p in paths:
        df = pd.read_csv(p, quotechar='"')
        df["ID"] = df.get("id", "").astype(str).apply(lambda s: (m := ID_RE.search(s)) and m.group(1))

        for col, gene in resistance_cols.items():
            if col not in df.columns:
                continue
            for _, r in df.iterrows():
                rid = r.get("ID")
                if not rid:
                    continue
                val = r.get(col, "")
                if pd.notna(val):
                    muts = [m.strip() for m in str(val).split(",") if m.strip() and m.strip().lower() != "not available"]
                    for mut in muts:
                        rows.append({"ID": rid, "gene": gene, "mutation": mut})

    out = pd.DataFrame(rows)
    out.to_csv(output_file, index=False)
    print(f"✅ HCV mutations: wrote {len(out):,} rows to {output_file}")
    return out

# ---------- File paths ----------
hbv_fasta = "hbv_sequences.fasta"
hbv_tree = "hbv_sequences_wrefs_aligned_edit_assigned.tree"
hbv_metadata = "hbv_metadata.csv"
hbv_output = "HBV_summary.tsv"  # keep case used by Dash

hcv_fasta = "hcv_sequences.fasta"
hcv_tree = "hcv_sequences_wrefs_aligned_edit_assigned.tree"
hcv_metadata = "hcv_metadata.csv"
hcv_output = "HCV_summary.tsv"

# mutation folders + outputs
hbv_g2p_dir = "hbv_g2p_dash"
hbv_mut_out = "hbv_resistant_mutations.csv"

hcv_g2p_dir = "hcv_g2p_dash"
hcv_mut_out = "hcv_resistance_mutations.csv"

# ---------- Run ----------
if __name__ == "__main__":
    # 1) Summary TSVs (genotypes, country, year)
    generate_summary_tsv(hbv_fasta, hbv_tree, hbv_metadata, HBV_COLOR_MAP, hbv_output, recombs_file="hbv_recombinants.txt")
    generate_summary_tsv(hcv_fasta, hcv_tree, hcv_metadata, HCV_COLOR_MAP, hcv_output, recombs_file="hcv_recombinants.txt")
    print("✅ TSV summary files created successfully.")

    # 2) Mutation tall tables for the dashboard
    extract_hbv_mutations(hbv_g2p_dir, hbv_mut_out)
    extract_hcv_mutations(hcv_g2p_dir, hcv_mut_out)
    print("✅ Mutation CSVs created successfully.")
    
import pandas as pd
import numpy as np
from pathlib import Path

# Optional: import your cached standardizer if it exists
try:
    from Hepatitis_dash_vagner_derek_darren_muts import standardize_country_cached
except Exception:
    def standardize_country_cached(x):  # fallback noop
        return x

# ------------- CONFIG -------------
HBV_SUMMARY = "HBV_summary.tsv"
HCV_SUMMARY = "HCV_summary.tsv"
IHME_FILE   = "IHME-GBD_2021_DATA-9e7ec2c0-1.csv"  # <-- adjust to your actual filename
OUTPUT_HBV  = "hbv_burden_adjusted_coverage.csv"
OUTPUT_HCV  = "hcv_burden_adjusted_coverage.csv"

CAUSE_LOOKUP = {
    "HBV": "Total burden related to hepatitis B",
    "HCV": "Total burden related to hepatitis C",
}

# ------------- CORE -------------
def _load_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # Expect columns: Taxa (Accession), Genotype, Country, Year (or Date -> Year already done in your pipeline)
    if "Country_standard" not in df.columns:
        df["Country_standard"] = df["Country"].apply(standardize_country_cached)
    # Ensure Year numeric if present
    if "Year" in df.columns:
        df["Year"] = pd.to_numeric(df["Year"], errors="coerce").astype("Int64")
    return df

def _load_ihme(path: str) -> pd.DataFrame:
    ihme = pd.read_csv(path)
    # Normalize columns
    rename_map = {}
    if "location" in ihme.columns:
        rename_map["location"] = "Country_standard"
    ihme = ihme.rename(columns=rename_map)
    ihme["Country_standard"] = ihme["Country_standard"].apply(standardize_country_cached)
    return ihme

def compute_burden_adjusted_coverage(
    virus: str,
    summary_df: pd.DataFrame,
    ihme_df: pd.DataFrame,
    by_year: bool = False
) -> pd.DataFrame:
    """
    virus: "HBV" or "HCV"
    summary_df: HBV/HCV summary with columns: Country_standard, Genotype [, Year]
    ihme_df: IHME with columns: measure, metric, sex, age, cause, year, val, Country_standard
    by_year: if True, compute per (country, year, genotype); else aggregate across years to country/genotype
    """
    virus = virus.upper()
    cause_filter = CAUSE_LOOKUP[virus]

    # ---- IHME: prevalence numbers for Both, All ages, selected virus ----
    ihme_prev = ihme_df[
        (ihme_df["cause"] == cause_filter) &
        (ihme_df["measure"] == "Prevalence") &
        (ihme_df["metric"] == "Number") &
        (ihme_df["sex"] == "Both") &
        (ihme_df["age"] == "All ages")
    ].copy()

    if ihme_prev.empty:
        print(f"⚠️ IHME prevalence (Number) is empty for {virus}. Check IHME file and CAUSE values.")
        # build an empty shell with expected columns
        cols = ["Country_standard", "Genotype", "Seq_count", "Total_seq",
                "Genotype_prop", "HBV_HCV_prevalence_abs", "Est_infections_genotype",
                "Coverage_ratio"]
        if by_year: cols.insert(1, "Year")
        return pd.DataFrame(columns=cols)

    # ---- Sequence counts & genotype proportions ----
    if by_year and "Year" in summary_df.columns and summary_df["Year"].notna().any():
        seq = (
            summary_df.dropna(subset=["Year"])
                      .groupby(["Country_standard", "Year", "Genotype"])
                      .size().reset_index(name="Seq_count")
        )
        totals = (
            seq.groupby(["Country_standard", "Year"])["Seq_count"]
               .sum().reset_index(name="Total_seq")
        )
        seq = seq.merge(totals, on=["Country_standard", "Year"], how="left")
        seq["Genotype_prop"] = np.where(
            seq["Total_seq"] > 0,
            seq["Seq_count"] / seq["Total_seq"],
            np.nan
        )
        # merge IHME by (country, year)
        prev = ihme_prev[["Country_standard", "year", "val"]].rename(columns={"year": "Year", "val": "Prev_Number"})
        out = seq.merge(prev, on=["Country_standard", "Year"], how="left")
    else:
        # aggregate across years to per-country
        seq = (
            summary_df.groupby(["Country_standard", "Genotype"])
                      .size().reset_index(name="Seq_count")
        )
        totals = (
            seq.groupby("Country_standard")["Seq_count"]
               .sum().reset_index(name="Total_seq")
        )
        seq = seq.merge(totals, on="Country_standard", how="left")
        seq["Genotype_prop"] = np.where(
            seq["Total_seq"] > 0,
            seq["Seq_count"] / seq["Total_seq"],
            np.nan
        )
        # For IHME, average prevalence across years (or pick a year). Here we use the latest year per country.
        prev_latest = (
            ihme_prev.sort_values("year")
                     .groupby("Country_standard", as_index=False)
                     .tail(1)  # latest row
        )
        prev_latest = prev_latest[["Country_standard", "val"]].rename(columns={"val": "Prev_Number"})
        out = seq.merge(prev_latest, on="Country_standard", how="left")

    # ---- Estimated infections per genotype & coverage ----
    out["Est_infections_genotype"] = out["Prev_Number"] * out["Genotype_prop"]
    # Avoid division by zero
    out["Coverage_ratio"] = np.where(
        (out["Est_infections_genotype"] > 0) & out["Seq_count"].notna(),
        out["Seq_count"] / out["Est_infections_genotype"],
        np.nan
    )

    # Tidy columns
    label_prev = f"{virus}_prevalence_abs"
    out = out.rename(columns={"Prev_Number": label_prev})
    ordered_cols = ["Country_standard"]
    if by_year: ordered_cols += ["Year"]
    ordered_cols += [
        "Genotype", "Seq_count", "Total_seq", "Genotype_prop",
        label_prev, "Est_infections_genotype", "Coverage_ratio"
    ]
    out = out[ordered_cols].sort_values(ordered_cols[: (2 if by_year else 1)] + ["Genotype"]).reset_index(drop=True)
    return out

def main(by_year: bool = False):
    hbv = _load_summary(HBV_SUMMARY)
    hcv = _load_summary(HCV_SUMMARY)
    ihme = _load_ihme(IHME_FILE)

    hbv_cov = compute_burden_adjusted_coverage("HBV", hbv, ihme, by_year=by_year)
    hcv_cov = compute_burden_adjusted_coverage("HCV", hcv, ihme, by_year=by_year)

    hbv_cov.to_csv(OUTPUT_HBV, index=False)
    hcv_cov.to_csv(OUTPUT_HCV, index=False)

    print(f"✅ Saved: {OUTPUT_HBV} ({len(hbv_cov):,} rows)")
    print(f"✅ Saved: {OUTPUT_HCV} ({len(hcv_cov):,} rows)")

if __name__ == "__main__":
    # Set by_year=True if you want (country, year, genotype) granularity
    main(by_year=False)

