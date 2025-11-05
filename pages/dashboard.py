import dash
from dash import html, register_page, dcc, callback, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from datetime import datetime
import pycountry
from functools import lru_cache
from dash import callback_context
from concurrent.futures import ThreadPoolExecutor
import re
import json
from dash.exceptions import PreventUpdate
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly_white"


register_page(__name__, path="/", name="Dashboard", order=0)
register_page(__name__, path="/dashboard", name="Dashboard")

# Import data loading functions
from data_loader import load_and_preprocess_data

# Initialize data store (will be loaded on first access)
data_store = None

def get_data_store():
    """Helper function to access the global data store"""
    global data_store
    if data_store is None:
        try:
            from data_loader import load_and_preprocess_data
            data_store = load_and_preprocess_data()
            print("Data store loaded successfully")
        except Exception as e:
            print(f"Error loading data: {e}")
            # Create empty data store structure to prevent further errors
            data_store = {
                'hbv_data': pd.DataFrame(),
                'hcv_data': pd.DataFrame(), 
                'ihme_df': pd.DataFrame(),
                'population_df': pd.DataFrame(),
                'coord_lookup': {},
                'hbv_mut': pd.DataFrame(),
                'hcv_mut': pd.DataFrame()
            }
    return data_store

# === CONFIG & CONSTANTS ======================================================
HBV_GENOTYPE_COLORS = {
    "HBV-A": "#cdac02",
    "HBV-B": "#951de0",
    "HBV-C": "#016301",
    "HBV-D": "#1f23bb",
    "HBV-E": "#770104",
    "HBV-F": "#7f7340",
    "HBV-G": "#da760c",
    "HBV-H": "#dea8b1",
    "HBV-I": "#e9e905",
    "HBV-J": "#07e705",
    "Recombinant": "#e00603"
}

HCV_GENOTYPE_COLORS = {
    "HCV-1": "#cdac02",
    "HCV-2": "#951de0",
    "HCV-3": "#016301",
    "HCV-4": "#1f23bb",
    "HCV-5": "#770104",
    "HCV-6": "#7f7340",
    "HCV-7": "#da760c",
    "HCV-8": "#FF69B4",
    "Recombinant": "#e00603"
}

BURDEN_MEASURE_FALLBACK = "Prevalence|Number"  # used if dropdown missing
        
# === HELPERS & FIGURE BUILDERS ==============================================
def compute_gap_df(
    virus: str,
    filtered_seq_df: pd.DataFrame,
    ihme_df: pd.DataFrame,
    selected_years: list | tuple,
    who_regions: list | None,
    countries: list | None,
    ihme_metric_choice: str | None,
    target_per_10k: float = 5.0
) -> pd.DataFrame:
    # Example placeholder logic
    if filtered_seq_df is None or filtered_seq_df.empty:
        return pd.DataFrame(columns=["Country_standard", "observed_sequences", "burden", "expected_sequences", "coverage_gap", "coverage_ratio"])

    # --- count observed sequences ---
    obs = (
        filtered_seq_df.groupby("Country_standard")
        .size()
        .reset_index(name="observed_sequences")
    )

    # --- get IHME burden for same filters ---
    try:
        measure, metric = (ihme_metric_choice or BURDEN_MEASURE_FALLBACK).split("|")
    except ValueError:
        measure, metric = "Prevalence", "Number"

    cause_lookup = {"HBV": "Total burden related to hepatitis B", "HCV": "Total burden related to hepatitis C"}
    cause = cause_lookup.get(virus.upper(), "")

    burden = ihme_df[
        (ihme_df["cause"] == cause)
        & (ihme_df["measure"] == measure)
        & (ihme_df["metric"] == metric)
        & (ihme_df["sex"] == "Both")
        & (ihme_df["age"] == "All ages")
    ].copy()

    # Filter by years
    if selected_years:
        y0, y1 = selected_years
        burden = burden[(burden["year"] >= y0) & (burden["year"] <= y1)]

    # Aggregate per country (latest year in range)
    if not burden.empty:
        burden = burden[burden["year"] == burden["year"].max()]
    burden = burden.groupby("Country_standard", as_index=False)["val"].sum().rename(columns={"val": "burden"})

    # --- merge and compute expected ---
    df = pd.merge(obs, burden, on="Country_standard", how="outer").fillna(0)
    
    # FIXED: Handle zero burden cases properly
    df["expected_sequences"] = np.where(
        df["burden"] > 0,
        (df["burden"] / 10000.0) * target_per_10k,
        0
    )
    
    # FIXED: Calculate coverage ratio safely
    df["coverage_ratio"] = np.where(
        df["expected_sequences"] > 0,
        df["observed_sequences"] / df["expected_sequences"],
        np.where(df["observed_sequences"] > 0, np.inf, 0)  # Infinite if we have sequences but no expected
    )
    
    df["coverage_gap"] = df["expected_sequences"] - df["observed_sequences"]
    df.loc[df["coverage_gap"] < 0, "coverage_gap"] = 0  # clip negative
    
    who_map = (
        ihme_df[["Country_standard", "WHO_Regions"]]
        .dropna(subset=["Country_standard"])
        .drop_duplicates("Country_standard")
    )
    
    df = df.merge(who_map, on="Country_standard", how="left")
    return df
    
def ihme_latest_by_country(ihme_df, virus, measure_metric, regions=None, countries=None, years=None):
    """
    Returns a DF with columns: Country_standard, Metric_raw, Metric, year
    filtered to Both sexes, All ages, selected cause, selected measure|metric,
    optional WHO regions / countries, and the latest year in the selected range.
    """
    try:
        ihme_measure, ihme_metric = (measure_metric or "Prevalence|Number").split("|")
    except Exception:
        ihme_measure, ihme_metric = "Prevalence", "Number"

    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
    }
    cause_filter = cause_lookup.get((virus or "HBV").upper(), "")

    base = ihme_df[
        (ihme_df["sex"] == "Both") &
        (ihme_df["age"] == "All ages") &
        (ihme_df["cause"] == cause_filter) &
        (ihme_df["measure"] == ihme_measure) &
        (ihme_df["metric"] == ihme_metric)
    ].copy()

    if regions:
        base = base[base["WHO_Regions"].isin(regions)]
    if countries:
        base = base[base["Country_standard"].isin(countries)]

    # Limit to selected range then pick latest year within that range
    if years and len(years) == 2:
        y0, y1 = int(years[0]), int(years[1])
        base = base[(base["year"] >= y0) & (base["year"] <= y1)]
    if base.empty:
        return pd.DataFrame(columns=["Country_standard", "Metric_raw", "Metric", "year"])

    latest_year = int(base["year"].max())
    latest = base[base["year"] == latest_year].copy()

    latest["Metric_raw"] = pd.to_numeric(latest["val"], errors="coerce")
    latest["Metric"] = latest["Metric_raw"].apply(lambda x: np.log10(x) if (np.isfinite(x) and x > 0) else np.nan)

    return latest[["Country_standard", "Metric_raw", "Metric", "year"]]

_ACCESSION_RE = re.compile(r"([A-Z0-9]+\.\d+)")

def _to_df(obj) -> pd.DataFrame:
    """Coerce list/dict/DF/None into a DataFrame copy."""
    if obj is None:
        return pd.DataFrame()
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, list):
        return pd.DataFrame(obj)
    if isinstance(obj, dict):
        # If it's a single record dict, wrap it
        if all(not isinstance(v, (list, tuple)) for v in obj.values()):
            return pd.DataFrame([obj])
        return pd.DataFrame(obj)
    return pd.DataFrame(obj)

def _rename_flex(df: pd.DataFrame) -> pd.DataFrame:
    """Case/alias-insensitive renames into canonical names used in plots."""
    if df.empty:
        return df
    col_map = {}
    for c in df.columns:
        lc = c.lower().strip()
        if lc in {"country_standard", "country_std", "country name", "country"}:
            col_map[c] = "Country_standard"
        elif lc in {"who_regions", "who region", "region", "who"}:
            col_map[c] = "WHO_Regions"
        elif lc in {"genotype", "genotypes", "geno"}:
            col_map[c] = "Genotype"
        elif lc in {"year", "yr"}:
            col_map[c] = "Year"
        elif lc in {"date"}:
            col_map[c] = "Date"
        elif lc in {"taxa"}:
            col_map[c] = "Taxa"
        elif lc in {"id", "accession"}:
            col_map[c] = "ID"
        elif lc in {"population"}:
            col_map[c] = "Population"
    return df.rename(columns=col_map)

def _ensure_year(df: pd.DataFrame) -> pd.DataFrame:
    """Create numeric Year from Year or Date, drop rows without it."""
    if "Year" not in df.columns and "Date" in df.columns:
        df["Year"] = pd.to_numeric(df["Date"], errors="coerce")
    if "Year" not in df.columns:
        df["Year"] = pd.NA
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["Year"]).copy()
    df["Year"] = df["Year"].astype(int)
    return df

def _normalize_seq_df(df: pd.DataFrame,
                      required=("Country_standard", "Year"),
                      fill_region=True) -> pd.DataFrame:
    """Normalize an input (list/dict/df) to the columns our plots expect."""
    df = _to_df(df)
    df = _rename_flex(df)
    df = _ensure_year(df)
    if fill_region and "WHO_Regions" not in df.columns:
        df["WHO_Regions"] = "Unknown"
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"normalize: missing {missing}; present columns: {list(df.columns)}"
        )
    return df

def _extract_id(series: pd.Series) -> pd.Series:
    """Extract accession-like ID (e.g., AB123456.1) from Taxa strings."""
    return series.astype(str).str.extract(_ACCESSION_RE, expand=False)

def _enrich_mutation_df(seq_df, seq_source_df) -> pd.DataFrame:
    """
    Ensure mutation DF has Country_standard/WHO_Regions/Year/Genotype by
    merging from the sequence source DF on ID (or Taxa-derived ID).
    """
    seq_df = _rename_flex(_to_df(seq_df))
    src = _normalize_seq_df(_rename_flex(_to_df(seq_source_df)),
                            required=("Country_standard", "Year", "Genotype"))

    # Make sure both frames have an ID for merging
    if "ID" not in seq_df.columns:
        if "Taxa" in seq_df.columns:
            seq_df["ID"] = _extract_id(seq_df["Taxa"])
    if "ID" not in src.columns:
        if "Taxa" in src.columns:
            src["ID"] = _extract_id(src["Taxa"])

    # If we already have all columns, return early
    needed = {"Country_standard", "WHO_Regions", "Year", "Genotype"}
    if needed.issubset(set(seq_df.columns)):
        return _normalize_seq_df(seq_df, required=tuple(needed))

    # Merge the missing fields from source
    attach_cols = ["ID", "Country_standard", "WHO_Regions", "Year", "Genotype"]
    attach_cols = [c for c in attach_cols if c in src.columns]
    if "ID" not in seq_df.columns or "ID" not in attach_cols:
        # Can't merge; raise a descriptive error
        raise KeyError("Mutation DF has no ID/Taxa to join on. Provide ID or Taxa.")
    merged = seq_df.merge(src[attach_cols].drop_duplicates("ID"),
                          on="ID", how="left", suffixes=("", "_src"))
    return _normalize_seq_df(merged, required=tuple(needed))


def merge_population_nearest(counts_df: pd.DataFrame,
                             pop_df: pd.DataFrame,
                             tol_years: int = 3) -> pd.DataFrame:
    """
    Per-country nearest-year merge within ±tol_years.
    Ensures Year is the same dtype (int64) on both sides and sorted.
    """
    if counts_df is None or len(counts_df) == 0:
        return counts_df.assign(Population=np.nan)

    # Copy and select req cols on population
    left = counts_df.copy()
    pop     = pop_df[["Country_standard", "Year", "Population"]].copy()

    # Coerce types consistently
    left["Country_standard"] = left["Country_standard"].astype(str)
    pop["Country_standard"]     = pop["Country_standard"].astype(str)

    # Coerce Year to numeric, drop NaNs, cast BOTH to EXACT SAME dtype (int64)
    left["Year"] = pd.to_numeric(left["Year"], errors="coerce")
    pop["Year"]     = pd.to_numeric(pop["Year"],  errors="coerce")
    left = left.dropna(subset=["Year"])
    pop     = pop.dropna(subset=["Year"])
    left["Year"] = left["Year"].astype("int64")
    pop["Year"]     = pop["Year"].astype("int64")

    # Ensure Population numeric
    pop["Population"] = pd.to_numeric(pop["Population"], errors="coerce")

    # Merge per country so keys are sorted within each group
    out = []
    for ctry, lgrp in left.groupby("Country_standard", sort=False):
        rgrp = pop[pop["Country_standard"] == ctry]
        if rgrp.empty:
            out.append(lgrp.assign(Population=np.nan))
            continue

        lgrp = lgrp.sort_values("Year", kind="mergesort")
        rgrp = rgrp.sort_values("Year", kind="mergesort")

        merged = pd.merge_asof(
            lgrp,
            rgrp,
            on="Year",                      # SAME dtype on both sides (int64)
            tolerance=int(tol_years),      # tolerance in years
            direction="nearest",
            allow_exact_matches=True,
        )
        out.append(merged)

    return pd.concat(out, ignore_index=True)



def merge_population_nearest_two_pass(counts_df: pd.DataFrame,
                                      pop_df: pd.DataFrame,
                                      tol_years_first: int = 3,
                                      tol_years_wide: int = 50) -> pd.DataFrame:
    """
    Tight nearest-year merge first; where Population is NaN, widen the window and fill.
    """
    first = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_first)

    # If Population missing anywhere, try a wider tolerance and fill only those rows
    if "Population" in first.columns and first["Population"].isna().any():
        widened = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_wide)
        first["Population"] = first["Population"].fillna(widened["Population"])

    return first

def _country_pie_heading(virus: str, years_text: str, top_n: int, has_filters: bool) -> str:
    scope = "current selection" if has_filters else "all data"
    return f"{virus} sequences by country · {years_text} — Top {top_n} ({scope})"
    
def _fmt_list(values, max_items=3, *, empty_label="All"):
    vals = [v for v in (values or []) if v]
    if not vals:
        return empty_label
    if len(vals) <= max_items:
        return ", ".join(vals)
    return f"{', '.join(vals[:max_items])} +{len(vals)-max_items} more"

def _mutations_heading(virus, years_text, filters_text, has_filters):
    scope = "current selection" if has_filters else "all data"
    v = (virus or "HBV").upper()
    prefix = "HBV resistance mutations" if v == "HBV" else "HCV mutations"
    facet  = "drug" if v == "HBV" else "gene"
    return f"{prefix} by {facet} — {filters_text} · {years_text} ({scope})"

def calculate_global_mutation_maximum():
    """Calculate the actual maximum mutation percentage from both HBV and HCV data"""
    data = get_data_store()
    
    max_percentage = 0
    
    # Check HBV mutations
    if not data['hbv_mut'].empty and not data['hbv_data'].empty:
        hbv_mutations = _enrich_mutation_df(data['hbv_mut'], data['hbv_data'])
        if not hbv_mutations.empty:
            hbv_counts = hbv_mutations.groupby("mutation")["ID"].nunique()
            hbv_total = len(data['hbv_data'])
            hbv_max = (hbv_counts.max() / hbv_total * 100) if hbv_total > 0 else 0
            max_percentage = max(max_percentage, hbv_max)
    
    # Check HCV mutations  
    if not data['hcv_mut'].empty and not data['hcv_data'].empty:
        hcv_mutations = _enrich_mutation_df(data['hcv_mut'], data['hcv_data'])
        if not hcv_mutations.empty:
            hcv_counts = hcv_mutations.groupby("mutation")["ID"].nunique()
            hcv_total = len(data['hcv_data'])
            hcv_max = (hcv_counts.max() / hcv_total * 100) if hcv_total > 0 else 0
            max_percentage = max(max_percentage, hcv_max)
    
    # Round up to nearest 10 and add some padding
    global_max = min(100, np.ceil(max_percentage / 10) * 10 + 10) if max_percentage > 0 else 50
    return global_max

def create_world_map(
    country_data: pd.DataFrame, 
    country_genotype_counts: pd.DataFrame, 
    coord_lookup: dict[str, dict[str, float]], 
    virus_type: str = "HBV", 
    display_mode: str = "raw"  # "raw", "PerMillion", or "ihme"
) -> go.Figure:
    
    genotype_colors = HBV_GENOTYPE_COLORS if virus_type == "HBV" else HCV_GENOTYPE_COLORS
    
    df = country_data.copy()
    if "Metric_raw" not in df.columns:
        df["Metric_raw"] = np.nan
    valid = df.dropna(subset=["Metric_raw"]).copy()

    fig = go.Figure()
    if valid.empty:
        return _empty_world("No country data available for current filters")

    # EPIDEMIOLOGY MODE (IHME data)
    if display_mode == "ihme":
        # Pre-defined global log limits for this specific metric and virus
        GLOBAL_LOG_MIN = 3.0  # 10^3 = 1,000
        GLOBAL_LOG_MAX = 8.0  # 10^8 = 100,000,000
    
        # Filter out zero values for log scaling
        valid_nonzero = valid[valid["Metric_raw"] > 0].copy()
        
        if valid_nonzero.empty:
            fig = _empty_world("No IHME data available for current filters")
            fig.add_annotation(text="Try adjusting year, age, or sex filters", x=0.5, y=0.4, showarrow=False)
            return fig
        
        # Apply logarithmic transformation
        valid_nonzero["log_value"] = np.log10(valid_nonzero["Metric_raw"])
        
        # Use appropriate color scale based on virus type
        HCV_COLOR_SCALE = [[0.0, "#FFF7BC"], [0.25, "#FEC44F"], [0.5, "#EC7014"], [0.75, "#993404"], [1.0, "#662506"]]
        HBV_COLOR_SCALE = [[0.0, "#F7FCF0"], [0.25, "#A8DDB5"], [0.5, "#2B8CBE"], [0.75, "#084081"], [1.0, "#06214D"]]
        
        colorscale = HBV_COLOR_SCALE if virus_type == "HBV" else HCV_COLOR_SCALE
        
        fig.add_trace(
            go.Choropleth(
                locations=valid_nonzero["Country_standard"],
                locationmode="country names",
                z=valid_nonzero["log_value"],
                # Use the GLOBAL min and max, not the data min/max
                zmin=GLOBAL_LOG_MIN,
                zmax=GLOBAL_LOG_MAX,
                colorscale=colorscale,
                colorbar=dict(
                    title="Prevalence (Log10 Scale)",
                    len=0.6,
                    thickness=20,
                    # Ticks are based on the global scale
                    tickvals=np.arange(GLOBAL_LOG_MIN, GLOBAL_LOG_MAX + 1),
                    ticktext=[f"10^{int(x)}" for x in np.arange(GLOBAL_LOG_MIN, GLOBAL_LOG_MAX + 1)],
                ),
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertemplate=(
                    "<b>%{location}</b><br>" +
                    "Prevalence: %{customdata:,.0f}<br>" +
                    "Log10: %{z:.2f}<extra></extra>"
                ),
                customdata=valid_nonzero["Metric_raw"].values,
            )
        )
    
    # PER MILLION MODE
    elif display_mode == "PerMillion":
        # log10 per‑million; keep as float ndarray
        z_vals = (
            valid["Metric_raw"].apply(lambda x: np.log10(x) if (pd.notna(x) and x > 0) else np.nan)
            .astype(float)
            .to_numpy()
        )
        if np.all(np.isnan(z_vals)):
            return _empty_world("No per‑million values available for current filters")
        vmin = float(np.nanmin(z_vals)) if np.isfinite(np.nanmin(z_vals)) else -5.0
        vmax = float(np.nanmax(z_vals)) if np.isfinite(np.nanmax(z_vals)) else 0.0

        HBV_COLOR_SCALE = [[0.0, "#deebf7"], [0.25, "#9ecae1"], [0.5, "#6baed6"], [0.75, "#3182bd"], [1.0, "#08519c"]]
        HCV_COLOR_SCALE = [[0.0, "#feedde"], [0.25, "#fdbe85"], [0.5, "#fd8d3c"], [0.75, "#e6550d"], [1.0, "#a63603"]]
        colorscale = HBV_COLOR_SCALE if virus_type == "HBV" else HCV_COLOR_SCALE

        fig.add_trace(
            go.Choropleth(
                locations=valid["Country_standard"],
                locationmode="country names",
                z=z_vals,
                zmin=vmin,
                zmax=vmax,
                colorscale=colorscale,
                colorbar_title="Log10 per million",
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertemplate="<b>%{location}</b><br>Per million: %{customdata:.2f}<extra></extra>",
                customdata=valid["Metric_raw"].astype(float),
            )
        )
    
    # RAW COUNTS MODE (default)
    else:
        # Discrete bins for raw counts; drive bounds with non‑zero countries
        bins = [0, 1, 5, 20, 100, 500, 2000, 4000, float("inf")]
        labels = ["0", "1–4", "5–19", "20–99", "100–499", "500–1,999", "2,000–3,999", "4,000+"]
        valid["bin"] = pd.cut(valid["Metric_raw"], bins=bins, labels=labels, include_lowest=True, right=False)
        bin_to_idx = {lab: i for i, lab in enumerate(labels)}
        valid["z_value"] = valid["bin"].map(bin_to_idx)

        nonzero = valid[valid["Metric_raw"] > 0].copy()
        driving = nonzero if not nonzero.empty else valid
        z_numeric = driving["z_value"].astype(float).to_numpy()
        locations = driving["Country_standard"]

        HCV_COLORS = ["#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506"]
        HBV_COLORS = ["#F7FCF0", "#E0F3DB", "#A8DDB5", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081", "#06214D"]
        colors = HBV_COLORS if virus_type == "HBV" else HCV_COLORS
        colorscale = [[i / (len(labels) - 1), c] for i, c in enumerate(colors)]

        fig.add_trace(
            go.Choropleth(
                locations=locations,
                locationmode="country names",
                z=z_numeric,
                zmin=0,
                zmax=len(labels) - 1,
                colorscale=colorscale,
                showscale=True,
                colorbar=dict(
                    tickvals=list(range(len(labels))),
                    ticktext=labels,
                    title="Sequence count",
                    len=0.6,
                    thickness=20,
                ),
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertext=driving.apply(
                    lambda r: f"<b>{r['Country_standard']}</b><br>Exact count: {float(r['Metric_raw']):.0f}<br>Range: {r['bin']}",
                    axis=1,
                ),
                hoverinfo="text",
            )
        )

    # Genotype overlay markers (only show for sequence data, not IHME data)
    if display_mode != "ihme":
        lons, lats, texts, sizes = [], [], [], []
        for country in df["Country_standard"].dropna().unique():
            subset = country_genotype_counts[country_genotype_counts["Country_standard"] == country]
            subset = subset[subset["Count"] > 0]
            total = int(subset["Count"].sum()) if not subset.empty else 0
            coords = coord_lookup.get(country)
            if total <= 0 or not coords:
                continue
            lat = float(coords["latitude"])
            lon = float(coords["longitude"])
            genotype_text = "<br>".join(
                f"{row.Genotype}: {row.Count} ({row.Count/total:.1%})" for _, row in subset.sort_values("Count", ascending=False).iterrows()
            )
            texts.append(f"<b>{country}</b><br>Total: {total}<br>{genotype_text}")
            lats.append(lat)
            lons.append(lon)
            sizes.append(10 + min(20, total ** 0.2))

        if lons:
            fig.add_trace(
                go.Scattergeo(
                    lon=lons,
                    lat=lats,
                    text=texts,
                    hoverinfo="text",
                    mode="markers",
                    marker=dict(size=sizes, color="lightgrey", opacity=0.7, line=dict(width=1.5, color="black")),
                    showlegend=False,
                )
            )

    # Apply to all modes after traces so fitbounds works
    fig.update_geos(
        projection_type="natural earth",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds="locations",
        domain=dict(x=[0, 1], y=[0.2, 1]),
    )
    fig.update_layout(height=700, margin=dict(t=30, b=30, l=10, r=10))
    return fig
    
def create_coverage_map(
    cov_df: pd.DataFrame,
    coord_lookup: dict[str, any],
    coords_df: pd.DataFrame | None = None,
    virus_type: str = "HBV",
    who_regions: list[str] | None = None,
    countries: list[str] | None = None,
) -> go.Figure:
    
    valid = cov_df.copy()

    # --- normalize basics ---
    if "Country_standard" in valid.columns:
        valid["Country_standard"] = valid["Country_standard"].astype(str).str.strip()
    
    # FIXED: Safe coverage ratio calculation
    if "coverage_ratio" not in valid.columns:
        # Calculate coverage ratio safely
        valid["coverage_ratio"] = np.where(
            (valid.get("expected_sequences", 0) > 0) & (valid.get("observed_sequences", 0) >= 0),
            valid.get("observed_sequences", 0) / valid.get("expected_sequences", 1),
            0
        )
    
    valid["coverage_ratio"] = pd.to_numeric(valid.get("coverage_ratio"), errors="coerce")
    
    # FIXED: Replace inf and handle zeros
    valid["coverage_ratio"] = valid["coverage_ratio"].replace([np.inf, -np.inf], 10)  # Cap infinite values
    valid["coverage_ratio"] = valid["coverage_ratio"].fillna(0)
    
    # --- resolve WHO region column robustly ---
    region_col = None
    for cand in ("WHO_Regions", "WHO_region", "who_regions", "who_region"):
        if cand in valid.columns:
            region_col = cand
            break
    
    # If we need WHO region filtering but the column is missing, attach from coords_df
    if who_regions and not region_col and coords_df is not None:
        add = (
            coords_df[["Country_standard", "WHO"]]
            .drop_duplicates("Country_standard")
            .rename(columns={"WHO": "WHO_Regions"})
        )
        valid = valid.merge(add, on="Country_standard", how="left")
        region_col = "WHO_Regions"
    
    # --- apply filters (only if the column exists) ---
    if who_regions and region_col:
        valid = valid[valid[region_col].isin(who_regions)]
    if countries and "Country_standard" in valid.columns:
        valid = valid[valid["Country_standard"].isin(countries)]
    
    # keep only rows with valid country + numeric coverage
    valid = valid[valid["Country_standard"].notna()]
    valid = valid[valid["coverage_ratio"].notna()]
    
    # FIXED: Better clipping for visualization
    valid["coverage_capped"] = valid["coverage_ratio"].clip(upper=5.0)  # Cap at 5x coverage

    # Discrete bins (left-closed, right-open) for coverage
    bins = [0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 5.01]
    labels = ["<0.10", "0.10–0.25", "0.25–0.50", "0.50–0.75", "0.75–1.0", "1.0–1.5", "1.5–2.0", "≥2.0"]
    valid["bin"] = pd.cut(valid["coverage_capped"], bins=bins, labels=labels, include_lowest=True, right=False)
    bin_to_idx = {lab: i for i, lab in enumerate(labels)}
    valid["z"] = valid["bin"].map(bin_to_idx).astype(float)

    # Color scale (cool→warm)
    HBV_COLOR_SCALE = [[0.0, "#deebf7"], [0.25, "#9ecae1"], [0.5, "#6baed6"], [0.75, "#3182bd"], [1.0, "#08519c"]]
    HCV_COLOR_SCALE = [[0.0, "#feedde"], [0.25, "#fdbe85"], [0.5, "#fd8d3c"], [0.75, "#e6550d"], [1.0, "#a63603"]]
    colorscale = HBV_COLOR_SCALE if virus_type == "HBV" else HCV_COLOR_SCALE

    fig = go.Figure()

    if valid.empty:
        fig.add_annotation(
            text="No coverage data for current filters",
            x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False
        )
        fitbounds_setting = False
    else:
        # --- ADD TRACE FIRST ---
        fig.add_trace(
            go.Choropleth(
                locations=valid["Country_standard"],
                locationmode="country names",
                z=valid["z"].astype(float).to_numpy(),
                zmin=0, zmax=len(labels) - 1,
                colorscale=colorscale,
                showscale=True,
                colorbar=dict(
                    tickvals=list(range(len(labels))),
                    ticktext=labels,
                    title="Coverage ratio",
                    len=0.6,
                    thickness=18,
                ),
                hovertext=valid.apply(
                    lambda r: (
                        f"<b>{r['Country_standard']}</b><br>"
                        f"Observed: {r.get('observed_sequences', 0):,.0f}<br>"
                        f"Expected: {r.get('expected_sequences', 0):,.0f}<br>"
                        f"Coverage: {r['coverage_ratio']:.2f}×"
                    ),
                    axis=1,
                ),
                hoverinfo="text",
                marker_line_color="rgba(0,0,0,0.25)",
                marker_line_width=0.5,
                name=""
            )
        )
        fitbounds_setting = "locations"

    # --- THEN FIT TO CURRENT LOCATIONS ---
    fig.update_layout(
        geo=dict(
            center=None,
            projection=dict(type="natural earth", scale=1),
        )
    )
    fig.update_geos(
        showcountries=True, countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds=fitbounds_setting,
        domain=dict(x=[0, 1], y=[0.2, 1]),
    )

    # Make the camera reset when filters change
    rev_parts = [
        virus_type,
        "WR:" + ",".join(sorted(who_regions or [])),
        "CT:" + ",".join(sorted(countries or [])),
        f"N:{len(valid)}"
    ]
    fig.update_layout(
        height=700,
        margin=dict(t=30, b=30, l=10, r=10),
        uirevision="|".join(rev_parts)
    )

    return fig
    
def make_line_trend(
    filtered_df: pd.DataFrame, 
    selected_virus: str, 
) -> go.Figure:
    
    genotype_colors = HBV_GENOTYPE_COLORS if selected_virus == "HBV" else HCV_GENOTYPE_COLORS
    
    # Group by year and genotype first
    line_data = (
        filtered_df.groupby(["Year", "Genotype"])
        .size().reset_index(name="Genome Sequences")
    )
    
    # Apply rolling average for each genotype
    smoothed_data = []
    for genotype in line_data["Genotype"].unique():
        genotype_df = line_data[line_data["Genotype"] == genotype].copy()
        genotype_df = genotype_df.sort_values("Year")
        
        # Apply 3-year rolling average
        genotype_df["Smoothed_Sequences"] = (
            genotype_df["Genome Sequences"]
            .rolling(window=3, min_periods=1, center=True)
            .mean()
        )
        smoothed_data.append(genotype_df)
    
    smoothed_df = pd.concat(smoothed_data, ignore_index=True)
    
    # FIXED: Handle zeros for log scale
    smoothed_df["Smoothed_Sequences"] = smoothed_df["Smoothed_Sequences"].replace(0, 0.1)  # Avoid log(0)
    
    fig = px.line(
        smoothed_df,
        x="Year",
        y="Smoothed_Sequences",
        color="Genotype",
        color_discrete_map=genotype_colors,
        markers=True,
        line_shape="spline"
    )
    
    fig.update_traces(mode='lines+markers', marker=dict(size=4))
    
    if not smoothed_df.empty:
        year_min = int(smoothed_df["Year"].min())
        year_max = int(smoothed_df["Year"].max())
        fig.update_layout(
            xaxis=dict(
                tickmode="array",
                tickvals=list(range(year_min, year_max+1, 4))
            )
        )
    
    # FIXED: Safe log range calculation
    y_pos = line_data["Genome Sequences"].replace(0, 0.1).dropna()  # Avoid zeros
    if len(y_pos):
        y_min = float(np.nanmax([0.1, np.nanmin(y_pos)]))  # Ensure positive
        y_max = float(np.nanmax(y_pos))
        lo = np.log10(max(0.1, y_min/1.5))
        hi = np.log10(y_max*1.5)
    else:
        lo, hi = 0, 2

    fig.update_layout(
        xaxis_title="Year",
        yaxis=dict(
            title="Number of Sequences (log scale)",
            type="log",
            autorange=False,
            range=[lo, hi],
            tickformat="~s",
            gridcolor="rgba(0,0,0,0.08)",
            linecolor="rgba(0,0,0,0.25)",
            zeroline=False
        ),
        height=400,
        margin=dict(t=100, b=0, l=0, r=0),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        plot_bgcolor="white",
        paper_bgcolor="white",
        meta={"y_full_log_range": [lo, hi]}
    )
    
    return fig

def make_genotype_bar(
    filtered_df: pd.DataFrame,
    population_df: pd.DataFrame,
    selected_virus: str,
    display_mode: str,
) -> go.Figure:
    """
    Totals-by-genotype bar (one bar per genotype), log y-axis.
    - Uses HBV_GENOTYPE_COLORS / HCV_GENOTYPE_COLORS.
    - Supports 'raw' (Total) and 'PerMillion'.
    - Keeps a legend (useful if you reuse the color map elsewhere).
    """
    # Pick palette
    selected_virus = selected_virus or "HBV"
    genotype_colors = HBV_GENOTYPE_COLORS if selected_virus.upper() == "HBV" else HCV_GENOTYPE_COLORS

    # 1) Aggregate counts at (Country, Year, Genotype)
    cyg = (
        filtered_df.groupby(["Country_standard", "Year", "Genotype"], as_index=False)
                   .size().rename(columns={"size": "Count"})
    )

    # 2) Merge population (same helper you used before)
    cyg = merge_population_nearest_two_pass(
        cyg, population_df, tol_years_first=3, tol_years_wide=50
    )

    # 3) Aggregate to genotype totals and population denominators
    agg = (
        cyg.groupby("Genotype", as_index=False)
           .agg(Total=("Count", "sum"), Pop=("Population", "sum"))
    )

    # Normalise genotype labels for consistent ordering and formatting
    def _norm(g):
        g = str(g).strip()
        if selected_virus.upper() == "HBV":
            if len(g) == 1 and g.isalpha():
                return f"HBV-{g.upper()}"
            elif g.startswith("HBV-"):
                return g
            elif g == "Recombinant":
                return "Recombinant"  # Keep as is for color mapping
            else:
                return f"HBV-{g}"
        else:  # HCV
            if g.isdigit() or (g.replace('.', '').isdigit() and g.count('.') <= 1):
                return f"HCV-{g}"
            elif g.startswith("HCV-"):
                return g
            elif g == "Recombinant":
                return "Recombinant"  # Keep as is for color mapping
            else:
                return f"HCV-{g}"

    agg["Genotype"] = agg["Genotype"].apply(_norm)

    # 4) Compute per-million when requested
    agg["PerMillion"] = np.where(
        (agg["Pop"].notna()) & (agg["Pop"] > 0),
        (agg["Total"] / agg["Pop"]) * 1_000_000.0,
        np.nan
    )
    y_col = "PerMillion" if (display_mode or "raw") == "PerMillion" else "Total"
    y_title = "Sequences per Million" if y_col == "PerMillion" else "Number of Sequences"

    # 5) Ensure stable order
    if selected_virus.upper() == "HBV":
        target_order = [f"HBV-{c}" for c in list("ABCDEFGHIJ")] + ["Recombinant"]
    else:
        # HCV genotypes: 1-7 plus recombinant
        target_order = [f"HCV-{i}" for i in range(1, 9)] + ["Recombinant"]

    # Add missing categories as zero so the axis is complete
    base = pd.DataFrame({"Genotype": target_order})
    agg = base.merge(agg, on="Genotype", how="left").fillna({y_col: 0, "Total": 0})

    # Remove any genotypes that don't exist in our data (all zeros)
    agg = agg[agg[y_col] > 0] if len(agg[agg[y_col] > 0]) > 0 else agg

    # 6) Prepare log range (avoid log(0))
    vals = pd.to_numeric(agg[y_col], errors="coerce").dropna()
    if len(vals) == 0 or vals.sum() == 0:
        # Empty data case
        ymin, ymax = 0.1, 10
    else:
        ymin = float(vals[vals > 0].min()) if (vals > 0).any() else 0.1
        ymax = float(vals.max()) if len(vals) else 10
    
    low_pad = max(0.1, ymin / 1.2)
    high_pad = max(1.0, ymax * 1.3)
    log_range = [np.log10(low_pad), np.log10(high_pad)]

    # 7) Create custom color mapping for the actual genotypes present
    present_genotypes = agg["Genotype"].tolist()
    custom_color_map = {}
    
    for genotype in present_genotypes:
        # Look for exact match first
        if genotype in genotype_colors:
            custom_color_map[genotype] = genotype_colors[genotype]
        # For recombinant - check if it's just "Recombinant" and use the color from dictionary
        elif genotype == "Recombinant" and "Recombinant" in genotype_colors:
            custom_color_map[genotype] = genotype_colors["Recombinant"]
        # For HBV genotypes without prefix
        elif selected_virus.upper() == "HBV" and len(genotype) == 1 and genotype.isalpha():
            hbv_key = f"HBV-{genotype}"
            if hbv_key in genotype_colors:
                custom_color_map[genotype] = genotype_colors[hbv_key]
        # For HCV genotypes without prefix
        elif selected_virus.upper() == "HCV" and genotype.isdigit():
            hcv_key = f"HCV-{genotype}"
            if hcv_key in genotype_colors:
                custom_color_map[genotype] = genotype_colors[hcv_key]
        else:
            # Fallback color if not found
            custom_color_map[genotype] = "#CCCCCC"

    # 8) Plot with custom color mapping
    bar_fig = px.bar(
        agg,
        x="Genotype",
        y=y_col,
        color="Genotype",
        category_orders={"Genotype": [g for g in target_order if g in agg["Genotype"].values]},
        color_discrete_map=custom_color_map,  # Use the custom mapping
        template="plotly_white",
    )

    # Hover text
    unit = " per million" if y_col == "PerMillion" else ""
    bar_fig.update_traces(
        hovertemplate="<b>%{x}</b><br>%{y:,}" + unit + " sequences<extra></extra>"
    )

    # Layout & legend
    bar_fig.update_layout(
        title=f"Total {selected_virus.upper()} Sequences by Genotype",
        xaxis_title="Genotype",
        yaxis=dict(
            title=y_title + " (log scale)",
            type="log",
            autorange=False,
            range=log_range,
            tickformat="~s",
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False,
            linecolor="rgba(0,0,0,0.25)",
        ),
        bargap=0.25,
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=420,
        margin=dict(t=70, b=40, l=40, r=20),
        legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center", yanchor="bottom"),
    )

    return bar_fig
    
def make_country_pie(df, selected_regions=None, selected_countries=None, selected_years=None,
                     virus="HBV", top_n=10):
    df = _normalize_seq_df(df, required=("Country_standard", "Year"))
    
    """
    Donut showing Top-N countries by sequence count.
    Does NOT include 'Other' in the pie chart - only shows actual countries.
    Adds a footnote like '+ 79 more countries (2,988 sequences not shown)'.
    Returns (fig, title).
    """

    # Apply region/country filters first
    if selected_regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(selected_regions)]
    if selected_countries:
        df = df[df["Country_standard"].isin(selected_countries)]

    # Compute the full year span *after* region/country filters
    year_series = pd.to_numeric(df["Year"], errors="coerce")
    data_ymin = int(year_series.min()) if year_series.notna().any() else None
    data_ymax = int(year_series.max()) if year_series.notna().any() else None

    # Apply year filter (if any) and build the label
    if selected_years and len(selected_years) == 2 and all(v is not None for v in selected_years):
        y0, y1 = int(selected_years[0]), int(selected_years[1])
        years_text = f"{y0}–{y1}"
        df = df[(year_series >= y0) & (year_series <= y1)]
        # A years filter counts only if it narrows the full span
        years_filter_active = (data_ymin is not None and data_ymax is not None) and (y0 > data_ymin or y1 < data_ymax)
    else:
        years_text = "All years"
        years_filter_active = False

    # Determine if *any* filters are active
    has_filters = bool(
        (selected_regions and len(selected_regions) > 0) or
        (selected_countries and len(selected_countries) > 0) or
        years_filter_active
    )

    if df.empty:
        fig = px.pie(pd.DataFrame({"Country": ["No data"], "Count": [1]}),
                     names="Country", values="Count", hole=0.6)
        fig.update_traces(textinfo="none", hoverinfo="skip", showlegend=False)
        fig.update_layout(margin=dict(l=10, r=10, t=40, b=40))
        title = _country_pie_heading(virus, years_text, top_n, has_filters)
        return fig, title

    # --- counts & figure (unchanged from your working version) ---
    vc = (df["Country_standard"].fillna("Unknown")
          .value_counts(dropna=False).rename_axis("Country").reset_index(name="Count")
          .sort_values("Count", ascending=False))
    top = vc.head(top_n).copy()
    other = vc.iloc[top_n:]
    other_countries, other_count = int(other.shape[0]), int(other["Count"].sum())
    total = int(vc["Count"].sum())
    top["SharePct"] = (top["Count"] / total) * 100

    fig = px.pie(top, names="Country", values="Count", hole=0.55,
                 category_orders={"Country": list(top["Country"])})
    fig.update_traces(
        text=[f"{p:.1f}%" for p in top["SharePct"]],
        texttemplate="%{text}", textinfo="text", textposition="inside", sort=False,
        marker=dict(line=dict(color="#fff", width=1)),
        hovertemplate="<b>%{label}</b><br>Sequences: %{value:,}<br>"
                      "Share of total: %{customdata:.1f}%<extra></extra>",
        customdata=top["SharePct"], showlegend=True,
    )
    fig.update_layout(
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="center", x=0.5, bgcolor="rgba(0,0,0,0)",
                    font=dict(size=11), itemwidth=30, title_text=""),
        margin=dict(l=20, r=20, t=90, b=80), autosize=True
    )
    footnote = (f"+ {other_countries} more countries ({other_count:,} sequences not shown)"
                if other_count > 0 and other_countries > 0
                else "All countries shown (complete data)")
    fig.add_annotation(text=footnote, x=0.5, y=-0.12, xref="paper", yref="paper",
                       showarrow=False, font=dict(size=14), align="center")

    title = _country_pie_heading(virus, years_text, min(top_n, len(vc)), has_filters)
    return fig, title

def make_mutation_bar(
    mutation_df: pd.DataFrame,
    total_sequences: int,
    selected_virus: str,
    selected_filter: list | None = None,
    years_range: list | tuple | None = None,
    data_span: list | tuple | None = None,
    other_filters_active: bool = False,
) -> tuple[go.Figure, str]:
    import numpy as np
    v = (selected_virus or "HBV").upper()
    is_hbv = (v == "HBV")
    color = "#3182bd" if is_hbv else "#d94801"
    col_name = "drug" if is_hbv else "gene"

    # years text + active
    if years_range and len(years_range) == 2 and None not in years_range:
        y0, y1 = int(years_range[0]), int(years_range[1])
        years_text = f"{y0}–{y1}"
        if data_span and len(data_span) == 2 and None not in data_span:
            ymin, ymax = int(data_span[0]), int(data_span[1])
            years_active = (y0 > ymin) or (y1 < ymax)
        else:
            years_active = True
    else:
        years_text = "All years"
        years_active = False

    # filter by drug/gene
    filters_text = "All drugs" if is_hbv else "All genes"
    if selected_filter:
        mutation_df = mutation_df[mutation_df[col_name].notna()]
        selected_filter_lower = [str(s).strip().lower() for s in selected_filter]
        mutation_df = mutation_df[
            mutation_df[col_name].astype(str).str.strip().str.lower().isin(selected_filter_lower)
        ]
        filters_text = _fmt_list(
            selected_filter, max_items=3,
            empty_label=("All drugs" if is_hbv else "All genes")
        )

    has_filters = bool(other_filters_active or years_active or (selected_filter and len(selected_filter) > 0))

    fig = go.Figure()
    if mutation_df.empty or not total_sequences:
        fig.update_layout(
            xaxis={"visible": False}, yaxis={"visible": False},
            annotations=[{"text": "No mutations found for current selection",
                          "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5,
                          "showarrow": False, "font": {"size": 16}}],
            height=450, plot_bgcolor="white", paper_bgcolor="white",
            margin=dict(t=40, b=0, l=0, r=0)
        )
        title = _mutations_heading(v, years_text, filters_text, has_filters)
        return fig, title

    # unique sequences per mutation
    mutation_counts = (
        mutation_df.groupby("mutation")["ID"].nunique().reset_index()
        .rename(columns={"mutation": "Mutation", "ID": "Unique_Sequences"})
    )
    mutation_counts["Proportion"] = (mutation_counts["Unique_Sequences"] / total_sequences * 100.0).round(2)
    
    # FILTER OUT MUTATIONS WITH 0% - KEY FIX
    mutation_counts = mutation_counts[mutation_counts["Proportion"] > 0.1]
    
    # If no mutations left after filtering, return empty plot
    if mutation_counts.empty:
        fig.update_layout(
            xaxis={"visible": False}, yaxis={"visible": False},
            annotations=[{"text": "No mutations with >0% frequency found",
                          "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5,
                          "showarrow": False, "font": {"size": 16}}],
            height=450, plot_bgcolor="white", paper_bgcolor="white",
            margin=dict(t=40, b=0, l=0, r=0)
        )
        title = _mutations_heading(v, years_text, filters_text, has_filters)
        return fig, title
    
    mutation_counts = mutation_counts.sort_values("Proportion", ascending=False).head(20)
    
    # Calculate global y-axis maximum
    GLOBAL_Y_MAX = calculate_global_mutation_maximum()
    
    # Get current maximum from the data
    current_max = mutation_counts["Proportion"].max() if not mutation_counts.empty else 0
    
    # Determine final y_max - use whichever is larger: global max or current max + padding
    if current_max > GLOBAL_Y_MAX:
        # If current data exceeds global max, round up to nearest 10 above current max
        y_max = min(100, np.ceil(current_max / 10) * 10 + 10)
    else:
        # Use the global maximum for consistent scaling
        y_max = GLOBAL_Y_MAX

    fig = px.bar(
        mutation_counts,
        x="Mutation", y="Proportion",
        labels={"Proportion": "Sequences with Mutation (%)", "Mutation": "Mutation"},
        color_discrete_sequence=[color],
        template="plotly_white"
    )
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>Percentage: %{y:.1f}%<br>"
                      "Count: %{customdata[0]}<br>"
                      "Total sequences: %{customdata[1]}<extra></extra>",
        customdata=np.stack([
            mutation_counts["Unique_Sequences"].astype(int).to_numpy(),
            np.full(len(mutation_counts), int(total_sequences))
        ], axis=1),
        texttemplate="%{y:.1f}%", textposition="outside", cliponaxis=False
    )
    
    # Update layout with consistent y-axis range
    fig.update_layout(
        height=450, 
        margin=dict(t=16, b=0, l=0, r=0),
        plot_bgcolor="white", 
        paper_bgcolor="white",
        xaxis=dict(
            showgrid=False, 
            linecolor="rgba(0,0,0,0.25)", 
            tickangle=-45, 
            automargin=True
        ),
        yaxis=dict(
            title="Sequences with Mutation (%)", 
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False, 
            range=[-2, y_max],  # Consistent upper limit
            fixedrange=True,
            # Consistent tick marks for both viruses
            tickmode='linear',
            tick0=0,
            dtick=10
        )
    )

    title = _mutations_heading(v, years_text, filters_text, has_filters)
    return fig, title
    
def make_coverage_bar(
    gap_df: pd.DataFrame,
    selected_virus: str,
    order: str = "lowest",
    top_n: int = 20
) -> go.Figure:
    """
    Create a bar chart showing coverage gaps or best coverage ratios.
    """
    virus = (selected_virus or "HBV").upper()
    
    if gap_df.empty:
        return go.Figure().update_layout(
            height=450,
            annotations=[dict(text="No data available for current filters",
                              x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white"
        )

    # --- normalize expected columns & dtypes ---
    # allow for different casings or missing columns
    colmap = {c.lower(): c for c in gap_df.columns}
    def col(name):    # case-insensitive getter
        return colmap.get(name.lower(), name)

    for c in ["observed_sequences", "expected_sequences", "coverage_gap"]:
        if col(c) in gap_df.columns:
            gap_df[col(c)] = pd.to_numeric(gap_df[col(c)], errors="coerce")

    # compute Coverage_Ratio if missing (case-insensitively)
    if "coverage_ratio" in colmap:
        gap_df["Coverage_Ratio"] = pd.to_numeric(gap_df[col("coverage_ratio")], errors="coerce")
    else:
        # create it from observed/expected
        obs = gap_df[col("observed_sequences")] if col("observed_sequences") in gap_df.columns else np.nan
        exp = gap_df[col("expected_sequences")] if col("expected_sequences") in gap_df.columns else np.nan
        gap_df["Coverage_Ratio"] = np.where(pd.to_numeric(exp, errors="coerce") > 0,
                                            pd.to_numeric(obs, errors="coerce") / pd.to_numeric(exp, errors="coerce"),
                                            np.nan)

    # basic clean
    if col("coverage_gap") in gap_df.columns:
        gap_df = gap_df.dropna(subset=[col("coverage_gap")])
    gap_df["Coverage_Ratio"] = pd.to_numeric(gap_df["Coverage_Ratio"], errors="coerce")

    # --- select / order top N ---
    if order == "lowest":
        # "under-sequenced" = biggest additional genomes needed
        gap_df = gap_df.sort_values(col("coverage_gap"), ascending=False).head(top_n)
        gap_df = gap_df.sort_values(col("coverage_gap"), ascending=True)  # for horizontal bar (small→large)
        x_col = col("coverage_gap")
        x_label = "Estimated additional genomes needed"
        title = f"Top {len(gap_df)} under-sequenced countries ({virus})"
    else:
        # "best covered" = highest coverage ratio
        gap_df = gap_df.sort_values("Coverage_Ratio", ascending=False).head(top_n)
        gap_df = gap_df.sort_values("Coverage_Ratio", ascending=True)      # small→large left→right
        x_col = "Coverage_Ratio"
        x_label = "Coverage ratio (Observed / Expected)"
        title = f"Top {len(gap_df)} best-covered countries ({virus})"

    color = "#3182bd" if virus == "HBV" else "#d94801"
    fig = px.bar(
        gap_df,
        x=x_col,
        y="Country_standard",
        orientation="h",
        color_discrete_sequence=[color],
        labels={x_col: x_label, "Country_standard": ""},
        title=title,
        template="plotly_white"
    )
    fig.update_layout(
        height=450,
        margin=dict(t=60, r=10, b=30, l=10),
        xaxis=dict(title=x_label, gridcolor="rgba(0,0,0,0.08)", linecolor="rgba(0,0,0,0.25)"),
        yaxis=dict(title=""),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return fig
    
def _empty_world(message: str) -> go.Figure:
    """Create an empty world map with a message."""
    fig = go.Figure()
    fig.update_layout(
        geo=dict(
            showframe=False,
            showcoastlines=True,
            projection_type='equirectangular'
        ),
        annotations=[dict(
            text=message,
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )],
        height=500
    )
    return fig

#Time Series with Projections
def create_forecast_chart(ihme_df, selected_virus, selected_regions=None, selected_countries=None):
    """Show historical trends with proper statistical forecasting using only numpy/pandas"""
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
    }
    cause_filter = cause_lookup.get(selected_virus.upper())
    
    if cause_filter not in ihme_df["cause"].values:
        return _empty_plot(f"No {selected_virus} burden data available for forecasting")
    
    # Filter data
    base = ihme_df[
        (ihme_df["sex"] == "Both") &
        (ihme_df["age"] == "All ages") &
        (ihme_df["cause"] == cause_filter) &
        (ihme_df["metric"] == "Number")
    ].copy()
    
    if selected_regions:
        base = base[base["WHO_Regions"].isin(selected_regions)]
    if selected_countries:
        base = base[base["Country_standard"].isin(selected_countries)]
    
    # Get the actual latest year from your data
    latest_data_year = base["year"].max()
    forecast_years = 8  # Forecast to 2030
    
    fig = go.Figure()
    
    # Colors for different measures
    measure_colors = {
        "Prevalence": "#1f77b4",
        "Incidence": "#ff7f0e", 
        "Deaths": "#d62728"
    }
    
    for measure in ["Prevalence", "Incidence", "Deaths"]:
        measure_data = base[base["measure"] == measure].copy()
        
        if measure_data.empty:
            continue
            
        # Aggregate by year
        yearly_data = measure_data.groupby("year")["val"].sum().reset_index()
        yearly_data = yearly_data.sort_values("year")
        
        if len(yearly_data) < 3:  # Need at least 3 points for meaningful forecast
            # Just show historical data if insufficient
            fig.add_trace(go.Scatter(
                x=yearly_data["year"],
                y=yearly_data["val"],
                name=f"{measure} (Insufficient Data)",
                line=dict(color=measure_colors[measure], width=3),
                mode='lines+markers',
                marker=dict(size=6)
            ))
            continue
        
        # Simple linear regression using numpy (y = mx + b)
        years = yearly_data["year"].values
        values = yearly_data["val"].values
        
        # Calculate linear regression manually
        n = len(years)
        sum_x = np.sum(years)
        sum_y = np.sum(values)
        sum_xy = np.sum(years * values)
        sum_x2 = np.sum(years * years)
        
        # Slope (m) and intercept (b)
        denominator = n * sum_x2 - sum_x * sum_x
        if denominator != 0:
            m = (n * sum_xy - sum_x * sum_y) / denominator
            b = (sum_y - m * sum_x) / n
            
            # Create future years for forecast
            future_years = np.arange(latest_data_year + 1, latest_data_year + forecast_years + 1)
            all_years = np.concatenate([years, future_years])
            
            # Predict including future years
            predictions = m * all_years + b
            predictions = np.maximum(predictions, 0)  # Ensure non-negative
            
            # Separate historical and forecast
            historical_pred = predictions[:len(years)]
            forecast_pred = predictions[len(years):]
            
            # Add historical data
            fig.add_trace(go.Scatter(
                x=years,
                y=values,
                name=f"{measure} (Historical)",
                line=dict(color=measure_colors[measure], width=3),
                mode='lines+markers',
                marker=dict(size=6)
            ))
            
            # Add forecast
            fig.add_trace(go.Scatter(
                x=future_years,
                y=forecast_pred,
                name=f"{measure} (Forecast)",
                line=dict(color=measure_colors[measure], width=2, dash='dash'),
                mode='lines',
                hovertemplate=f"{measure} Forecast: %{{y:,.0f}}<extra></extra>"
            ))
            
            # Add simple confidence interval (±15%)
            upper_bound = forecast_pred * 1.15
            lower_bound = forecast_pred * 0.85
            
            fig.add_trace(go.Scatter(
                x=np.concatenate([future_years, future_years[::-1]]),
                y=np.concatenate([upper_bound, lower_bound[::-1]]),
                fill='toself',
                fillcolor=f'rgba{tuple(int(measure_colors[measure].lstrip("#")[i:i+2], 16) for i in (0, 2, 4)) + (0.2,)}',
                line=dict(color='rgba(255,255,255,0)'),
                name=f"{measure} Forecast Range",
                showlegend=False,
                hoverinfo='skip'
            ))
            
        else:
            # If linear regression fails, just show historical data
            fig.add_trace(go.Scatter(
                x=years,
                y=values,
                name=f"{measure} (Historical)",
                line=dict(color=measure_colors[measure], width=3),
                mode='lines+markers'
            ))
    
    if len(fig.data) == 0:
        return _empty_plot("Insufficient data for forecasting")
    
    # Add vertical line for current data boundary
    fig.add_vline(
        x=latest_data_year, 
        line_dash="dot", 
        line_color="red",
        line_width=2
    )
    
    # Add annotation for data boundary
    fig.add_annotation(
        x=latest_data_year,
        y=0.02,  # Position near bottom
        yref="paper",
        text="<b>Data Limit</b>",
        showarrow=False,
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="red",
        borderwidth=1,
        borderpad=4,
        font=dict(color="red", size=11),
        xanchor="right"  # Anchor to right to avoid overlap
    )
    
    # WHO elimination target
    who_target_year = 2030
    
    # Add the vertical line for WHO target
    fig.add_vline(
        x=who_target_year, 
        line_dash="dot", 
        line_color="green",
        line_width=2
    )
    
    # Add a clean annotation for WHO target
    fig.add_annotation(
        x=who_target_year,
        y=0.98,  # Near top of plot area
        yref="paper",
        text="<b>WHO 2030<br>Target</b>",
        showarrow=False,
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="green",
        borderwidth=1,
        borderpad=4,
        font=dict(color="green", size=11),
        xanchor="left"  # Anchor text to left
    )
    
    # Calculate appropriate x-axis range
    if len(fig.data) > 0:
        # Get all x values from traces
        all_x = []
        for trace in fig.data:
            if hasattr(trace, 'x') and trace.x is not None:
                all_x.extend([x for x in trace.x if x is not None])
        
        if all_x:
            x_min = min(all_x)
            x_max = max(max(all_x), who_target_year)  # Ensure WHO target is visible
    else:
        x_min = latest_data_year - 10
        x_max = who_target_year + 2
    
    fig.update_layout(
        xaxis_title="Year",
        yaxis_title="Number of Cases",
        yaxis_type="log",
        height=400,
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        xaxis=dict(
            range=[x_min - 1, x_max + 1] if 'x_min' in locals() and 'x_max' in locals() else [2010, 2032]
        )
    )
    
    return fig

#Mutation Timeline
def create_mutation_timeline(mutation_df, sequence_df, selected_virus, top_mutations=10):
    """Show emergence and spread of key mutations over time"""
    
    # Enrich mutation data with temporal information
    enriched_mutations = _enrich_mutation_df(mutation_df, sequence_df)
    
    if enriched_mutations.empty:
        return _empty_plot("No mutation data available")
    
    # Get top mutations by frequency
    top_muts = (enriched_mutations.groupby("mutation")["ID"]
                .nunique()
                .nlargest(top_mutations)
                .index.tolist())
    
    # Aggregate by year and mutation
    timeline_data = (enriched_mutations[enriched_mutations["mutation"].isin(top_muts)]
                    .groupby(["Year", "mutation"])
                    .size()
                    .reset_index(name="count"))
    
    # Calculate cumulative prevalence
    yearly_totals = sequence_df.groupby("Year").size().reset_index(name="total_sequences")
    timeline_data = timeline_data.merge(yearly_totals, on="Year", how="left")
    timeline_data["prevalence_pct"] = (timeline_data["count"] / timeline_data["total_sequences"]) * 100
    
    fig = px.scatter(timeline_data, 
                     x="Year", 
                     y="prevalence_pct",
                     color="mutation",
                     size="count",
                     hover_data={"count": True, "prevalence_pct": ":.2f"})
    
    # Add lines connecting the points
    for mutation in top_muts:
        mutation_data = timeline_data[timeline_data["mutation"] == mutation]
        fig.add_trace(go.Scatter(
            x=mutation_data["Year"],
            y=mutation_data["prevalence_pct"],
            mode='lines',
            line=dict(width=1, color='lightgray'),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    fig.update_layout(
        yaxis_title="Prevalence (%)",
        xaxis_title="Year",
        height=400,
        hovermode="closest"
    )
    
    return fig

#Transmission Cluster Map
def create_transmission_clusters(sequence_df, selected_virus, genetic_distance_threshold=0.05):
    """Identify and visualize transmission clusters based on temporal and geographic patterns"""
    
    if sequence_df.empty:
        return _empty_plot("No sequence data available for cluster analysis")
    
    # Use temporal and geographic patterns since we don't have genetic distance data
    cluster_data = sequence_df.copy()
    
    # Create clusters based on country and year proximity (realistic proxy)
    # Countries with sequences in consecutive years might indicate transmission clusters
    
    # Group by country and count sequences per year
    country_year_counts = (cluster_data.groupby(["Country_standard", "Year"])
                          .size()
                          .reset_index(name="sequence_count"))
    
    # Identify countries with increasing sequence counts (potential outbreaks)
    clusters = []
    for country in country_year_counts["Country_standard"].unique():
        country_data = country_year_counts[country_year_counts["Country_standard"] == country].sort_values("Year")
        
        if len(country_data) >= 2:
            # Calculate year-over-year growth
            country_data["growth"] = country_data["sequence_count"].pct_change()
            
            # Flag as cluster if significant growth detected
            significant_growth = country_data[country_data["growth"] > 0.5]  # 50% growth threshold
            
            for _, row in significant_growth.iterrows():
                clusters.append({
                    "Country_standard": country,
                    "Year": row["Year"],
                    "sequences": row["sequence_count"],
                    "growth_pct": row["growth"] * 100,
                    "cluster_type": "Emerging" if row["growth"] > 1.0 else "Growing"
                })
    
    if not clusters:
        # Fallback: show countries with highest recent sequencing activity
        recent_year = cluster_data["Year"].max()
        recent_data = cluster_data[cluster_data["Year"] == recent_year]
        if not recent_data.empty:
            country_counts = recent_data["Country_standard"].value_counts().head(10)
            for country, count in country_counts.items():
                clusters.append({
                    "Country_standard": country,
                    "Year": recent_year,
                    "sequences": count,
                    "growth_pct": 0,
                    "cluster_type": "Active Sequencing"
                })
    
    clusters_df = pd.DataFrame(clusters)
    
    if clusters_df.empty:
        return _empty_plot("No transmission patterns detected with current data")
    
    # Create cluster map with different colors for cluster types
    fig = px.scatter_geo(clusters_df,
                        locations="Country_standard",
                        locationmode="country names",
                        size="sequences",
                        color="cluster_type",
                        hover_name="Country_standard",
                        hover_data={
                            "Year": True, 
                            "sequences": True,
                            "growth_pct": ":.1f",
                            "cluster_type": True
                        },
                        color_discrete_map={
                            "Emerging": "#ff4444",
                            "Growing": "#ffaa00", 
                            "Active Sequencing": "#44ff44"
                        },
                        title=f"{selected_virus} Transmission Patterns and Sequencing Activity")
    
    fig.update_geos(
        projection_type="natural earth",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        fitbounds="locations"
    )
    
    fig.update_layout(
        height=500,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        )
    )
    
    return fig

#Replace Pie Chart with Treemap:
def create_country_treemap(df, selected_virus, selected_regions=None):
    """Treemap showing sequences by region > country"""
    
    if df.empty:
        return _empty_plot("No data available")
    
    # Aggregate data
    if selected_regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(selected_regions)]
    
    region_country_counts = (df.groupby(["WHO_Regions", "Country_standard"])
                            .size()
                            .reset_index(name="count"))
    
    fig = px.treemap(region_country_counts,
                     path=["WHO_Regions", "Country_standard"],
                     values="count",
                     color="count",
                     color_continuous_scale='Blues',
                     title=f"{selected_virus} Sequences by Region and Country")
    
    fig.update_layout(height=400)
    fig.update_traces(
        hovertemplate='<b>%{label}</b><br>Sequences: %{value:,}<br>Parent: %{parent}'
    )
    
    return fig

#Country barchart
def create_country_stacked_bar(df, selected_virus, selected_regions=None, selected_countries=None, top_n=10):
    """Create a stacked bar chart showing sequences by country with genotype breakdown"""
    
    if df.empty:
        return _empty_plot("No data available")
    
    # Apply region/country filters
    if selected_regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(selected_regions)]
    if selected_countries:
        df = df[df["Country_standard"].isin(selected_countries)]
    
    # Get top countries by total sequences
    country_totals = df.groupby('Country_standard').size().sort_values(ascending=False)
    top_countries = country_totals.head(top_n).index
    
    # Filter to top countries
    top_df = df[df['Country_standard'].isin(top_countries)].copy()
    
    # Aggregate by country and genotype
    stacked_data = (top_df.groupby(['Country_standard', 'Genotype'])
                    .size()
                    .reset_index(name='count'))
    
    if stacked_data.empty:
        return _empty_plot("No data available for stacked bar chart")
    
    # Get genotype colors based on virus type
    genotype_colors = HBV_GENOTYPE_COLORS if selected_virus.upper() == "HBV" else HCV_GENOTYPE_COLORS
    
    # Create stacked bar chart
    fig = px.bar(
        stacked_data,
        x='Country_standard',
        y='count',
        color='Genotype',
        labels={'count': 'Number of Sequences', 'Country_standard': 'Country'},
        color_discrete_map=genotype_colors
    )
    
    # Update layout for better readability
    fig.update_layout(
        height=500,
        xaxis_title="Country",
        yaxis_title="Number of Sequences",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            title="Genotype"
        ),
        margin=dict(r=150),  # Add margin for legend
        xaxis={'categoryorder': 'total descending'}  # Sort by total sequences
    )
    
    # Update hover template
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>Genotype: %{fullData.name}<br>Sequences: %{y:,}<extra></extra>"
    )
    
    return fig

#Priority Setting Tool
def create_priority_calculator(gap_df, ihme_df, selected_virus, weights=None):
    """Rank countries for sequencing investment based on multiple criteria"""
    
    if gap_df.empty:
        return _empty_plot("No data available for priority calculation"), pd.DataFrame()
    
    default_weights = {
        "burden": 0.4,
        "coverage_gap": 0.3, 
        "population": 0.2,
        "neighbor_sequencing": 0.1
    }
    weights = weights or default_weights
    
    # Calculate priority scores
    priority_data = gap_df.copy()
    
    # Normalize metrics (handle missing columns safely)
    for metric in ["burden", "coverage_gap", "observed_sequences"]:
        if metric in priority_data.columns:
            col_min = priority_data[metric].min()
            col_max = priority_data[metric].max()
            if col_max > col_min:  # Avoid division by zero
                priority_data[f"{metric}_norm"] = (priority_data[metric] - col_min) / (col_max - col_min)
            else:
                priority_data[f"{metric}_norm"] = 0.5  # Default value if all values are same
    
    # Proper population normalization
    population_col = priority_data.get("Population", pd.Series([1] * len(priority_data)))
    if hasattr(population_col, 'max'):    # Check if it's a Series with max method
        pop_max = population_col.max()
        population_norm = population_col / max(pop_max, 1)    # Avoid division by zero
    else:
        population_norm = 0     # Fallback if Population is not available
    
    # Calculate composite score
    priority_data["priority_score"] = (
        weights["burden"] * priority_data.get("burden_norm", 0) +
        weights["coverage_gap"] * priority_data.get("coverage_gap_norm", 0) +
        weights["population"] * population_norm
    )
    
    # Rank countries
    priority_data = priority_data.sort_values("priority_score", ascending=False)
    priority_data["rank"] = range(1, len(priority_data) + 1)
    
    # Create interactive table
    fig = go.Figure(data=[go.Table(
        header=dict(values=["Rank", "Country", "Priority Score", "Burden", "Coverage Gap", "Sequences"],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[priority_data["rank"], 
                          priority_data["Country_standard"],
                          priority_data["priority_score"].round(3),
                          priority_data.get("burden", 0).round(0),
                          priority_data.get("coverage_gap", 0).round(0),
                          priority_data.get("observed_sequences", 0)],
                   align='left'))
    ])
    
    fig.update_layout(
        height=400
    )
    
    return fig, priority_data

#---burden-lines-helper-----

MASTER_TICKS = np.array([3e5, 6e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7,
                         1e8, 2e8, 5e8], dtype=float)
MASTER_TEXT     = ["300k","600k","1M","2M","5M","10M","20M","50M",
                "100M","200M","500M"]

def virus_log_axis(selected_virus: str, values: np.ndarray) -> dict:
    """Return {range, tickvals, ticktext} for a single log axis tuned per virus."""
    v = (selected_virus or "HBV").upper()
    # sensible presets reflecting your data ranges
    presets = {
        "HBV": dict(floor=3e5, cap=6e8),   # 70M–400M + deaths ~500–700k
        "HCV": dict(floor=2e5, cap=3e8),   # 9M–200M + deaths ~200–600k
    }
    p = presets.get(v, dict(floor=2e5, cap=9e8))

    vals = np.asarray(values, float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        lo, hi = p["floor"], p["cap"]
    else:
        lo = max(p["floor"], np.nanmin(vals) * 0.9)        # small padding
        hi = min(p["cap"],     np.nanmax(vals) * 1.1)

    # snap to decades so the axis looks clean
    lo = 10 ** np.floor(np.log10(lo))
    hi = 10 ** np.ceil( np.log10(hi))

    mask = (MASTER_TICKS >= lo) & (MASTER_TICKS <= hi)
    tickvals = MASTER_TICKS[mask].tolist()
    ticktext = [t for t,m in zip(MASTER_TEXT, mask) if m]

    return dict(
        range=[float(np.log10(lo)), float(np.log10(hi))],
        tickvals=tickvals,
        ticktext=ticktext
    )

def make_burden_lines(
    ihme_df: pd.DataFrame,
    selected_virus: str,
    selected_continents: list = None,
    selected_countries: list = None,
    selected_years: list = None
) -> tuple[go.Figure, str]:
    
    """
    Create a line chart showing IHME burden data (prevalence, incidence, deaths) on single y-axis.
    Returns figure and title.
    """
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
    }
    cause_filter = cause_lookup.get((selected_virus or "HBV").upper())
    
    base = ihme_df[
        (ihme_df["sex"] == "Both") &
        (ihme_df["age"] == "All ages") &
        (ihme_df["cause"] == cause_filter)
    ].copy()
    
    # Apply year filtering if provided
    if selected_years and len(selected_years) == 2:
        y0, y1 = map(int, selected_years)
        base = base[(base["year"] >= y0) & (base["year"] <= y1)]
    
    if selected_continents:
        base = base[base["WHO_Regions"].isin(selected_continents)]
    if selected_countries:
        base = base[base["Country_standard"].isin(selected_countries)]

    def series(measure, metric="Number"):
        s = base[(base["measure"] == measure) & (base["metric"] == metric)].copy()
        if s.empty:
            return pd.DataFrame(columns=["year", f"{measure}_{metric}"])
        s["val"] = pd.to_numeric(s["val"], errors="coerce")
        s = s.dropna(subset=["val"])
        return (s.groupby("year", as_index=False)["val"]
                  .sum().rename(columns={"val": f"{measure}_{metric}"}))

    d_prev = series("Prevalence", "Number")
    d_inc  = series("Incidence",  "Number")
    d_dea  = series("Deaths",      "Number")

    from functools import reduce
    dfs = [d for d in [d_prev, d_inc, d_dea] if not d.empty]
    if dfs:
        burden_df = reduce(lambda L,R: pd.merge(L,R,on="year",how="outer"), dfs).sort_values("year")
    else:
        burden_df = pd.DataFrame(columns=["year","Prevalence_Number","Incidence_Number","Deaths_Number"])

    # Ensure numeric and mask non-positives for log plot
    for col in ["Prevalence_Number","Incidence_Number","Deaths_Number"]:
        if col in burden_df:
            burden_df[col] = pd.to_numeric(burden_df[col], errors="coerce")

    def ymask(col):
        return burden_df[col].mask(~(burden_df[col] > 0), None)

    # Build the figure
    burden_fig = go.Figure()

    # Prevalence
    if "Prevalence_Number" in burden_df:
        burden_fig.add_trace(go.Scatter(
            x=burden_df["year"], y=ymask("Prevalence_Number"),
            mode="lines+markers", name="Prevalence",
            hovertemplate="Year: %{x}<br>Prevalence: %{y:,.0f}<extra></extra>",
            line=dict(width=3)
        ))

    # Incidence
    if "Incidence_Number" in burden_df:
        burden_fig.add_trace(go.Scatter(
            x=burden_df["year"], y=ymask("Incidence_Number"),
            mode="lines+markers", name="Incidence",
            hovertemplate="Year: %{x}<br>Incidence: %{y:,.0f}<extra></extra>",
            line=dict(width=3)
        ))

    # Deaths
    if "Deaths_Number" in burden_df:
        burden_fig.add_trace(go.Scatter(
            x=burden_df["year"], y=ymask("Deaths_Number"),
            mode="lines+markers", name="Deaths",
            hovertemplate="Year: %{x}<br>Deaths: %{y:,.0f}<extra></extra>",
            line=dict(width=3)
        ))

    # Get values for axis scaling
    vals_for_axis = []
    for col in ("Prevalence_Number","Incidence_Number","Deaths_Number"):
        if col in burden_df.columns:
            vals_for_axis.append(burden_df[col].to_numpy())
    vals_for_axis = np.concatenate(vals_for_axis) if vals_for_axis else np.array([])
    
    # Use the helper function for consistent axis scaling
    axis = virus_log_axis(selected_virus, vals_for_axis)
    
    # Layout with single y-axis
    burden_fig.update_layout(
        height=400,
        margin=dict(t=50, b=20, l=60, r=20),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        xaxis=dict(
            title="Year",
            dtick=3, tickmode="linear",
            gridcolor="rgba(0,0,0,0.06)", showgrid=True
        ),
        yaxis=dict(
            title="Number of Cases (log scale)",
            type="log",
            autorange=False,
            range=axis["range"],
            tickmode="array",
            tickvals=axis["tickvals"],
            ticktext=axis["ticktext"],
            gridcolor="rgba(0,0,0,0.05)",
            zeroline=False
        ),
        plot_bgcolor="white",
        paper_bgcolor="white"
    )
    
    burden_title = f"Global burdens in {selected_virus}: Prevalence, Incidence & Deaths"
    
    return burden_fig, burden_title
    
def _empty_plot(message):
    """Create an empty plot with a message"""
    fig = go.Figure()
    fig.update_layout(
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[{
            "text": message,
            "xref": "paper", "yref": "paper",
            "x": 0.5, "y": 0.5, "showarrow": False,
            "font": {"size": 16}
        }],
        height=300
    )
    return fig

def _df_to_json(df: pd.DataFrame) -> str:
    if df is None or df.empty:
        return json.dumps({"columns": [], "data": []})
    return json.dumps({"columns": list(df.columns), "data": df.to_dict("records")})

def _df_from_json(payload: str) -> pd.DataFrame:
    try:
        obj = json.loads(payload or "{}")
        return pd.DataFrame(obj.get("data", []), columns=obj.get("columns", []))
    except Exception:
        return pd.DataFrame()

# === HELPERS FOR LAYOUT =====================================================
def build_indicators(virus):
    color_class = "bg-primary" if virus == "HBV" else "bg-warning"

    return dbc.Col([
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Total Whole Genomes", className="card-subtitle"),
                                html.H4(id="indicator-total", className="card-title")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-bar-chart-fill", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Countries"),
                                html.H4(id="indicator-countries")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-globe", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Genotypes"),
                                html.H4(id="indicator-genotypes")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-diagram-3-fill", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Years"),
                                html.H4(id="indicator-years")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-hourglass-split", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

        ], className="mb-4 shadow-sm"),
    ], className="mb-4", width=3)

# === APP SETUP ==============================================================
external_stylesheets = [
    "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css",
    dbc.themes.BOOTSTRAP
]

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

# === APP LAYOUT ==============================================================
def create_dashboard_layout():
    return dbc.Container([
        # --- page-level state stores ---
        dcc.Store(id="selected-virus", data="HBV"),
        dcc.Store(id="filtered-store"),
        dcc.Store(id="gap-store"),
        dcc.Store(id="computed-metrics-store"),
        dcc.Store(id="ihme-latest-store"),
        dcc.Store(id="priority-data-store"),
        dcc.Download(id="priority-download"),
        
        # Header Section
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Img(src="/assets/logo.png", height="40px", className="me-2"),
                    html.H4("Hepatitis Virus Sequence Dashboard", className="mb-0")
                ], className="d-flex align-items-center")
            ], width=4),
            
            dbc.Col([
                html.Div([
                    dbc.Button("HBV", id="btn-hbv", color="primary", size="lg", 
                              className="me-2", outline=False, n_clicks=1),
                    dbc.Button("HCV", id="btn-hcv", color="warning", size="lg", 
                              outline=True, n_clicks=0)
                ])
            ], width="auto"),
            
            dbc.Col([
                dbc.DropdownMenu(
                    label="Download",
                    children=[
                        dbc.DropdownMenuItem("Data", id="btn-download-data"),
                        dbc.DropdownMenuItem("Reports", href="/about#report-section"),
                    ],
                    color="success",
                    size="lg",
                    className="me-2"
                ),
                dbc.Toast("Your download is starting…", id="dl-toast", header="Download", is_open=False,
                  dismissable=True, icon="success", duration=3000, className="position-fixed top-0 end-0 m-3"),
                html.Div(id="download-trigger", style={"display": "none"}),
                dcc.Download(id="download-data"),
            ], width="auto"),
        ], className="mb-4 align-items-center", justify="between"),

        # Filters Section
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("FILTERS", className="text-muted mb-3"),
                        
                        dbc.Row([
                            dbc.Col([
                                html.Label("Year Range", className="fw-bold mb-2"),
                                html.Div(id="selected-years-display", className="text-primary fw-bold mb-2"),
                                dcc.RangeSlider(
                                    id="year-slider",
                                    step=1,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    className="mb-3"
                                ),
                                dbc.Button("Reset Time Range", id="btn-reset-time", 
                                         color="secondary", outline=True, size="sm")
                            ], width=12)
                        ], className="mb-3"),
                        
                        dbc.Row([
                            dbc.Col([
                                html.Label("WHO Region(s)", className="fw-bold mb-2"),
                                dcc.Dropdown(id="continent-dropdown", multi=True, 
                                           placeholder="All Regions", className="mb-3")
                            ], width=4),
                            
                            dbc.Col([
                                html.Label("Country(s)", className="fw-bold mb-2"),
                                dcc.Dropdown(id="country-dropdown", multi=True, 
                                           placeholder="All Countries", className="mb-3")
                            ], width=4),
                            
                            dbc.Col([
                                html.Label("Genotype(s)", className="fw-bold mb-2"),
                                dcc.Dropdown(id="genotype-dropdown", multi=True, 
                                           placeholder="Select virus first", className="mb-3")
                            ], width=4)
                        ])
                    ])
                ], className="mb-4 shadow-sm")
            ], width=12)
        ]),

        # Display Toggles
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                html.Label("Display Mode:", className="fw-bold me-2"),
                                dcc.RadioItems(
                                    id="display-mode",
                                    options=[
                                        {"label": " Raw Count", "value": "raw"},
                                        {"label": " Per Million", "value": "PerMillion"}
                                    ],
                                    value="raw",
                                    inline=True,
                                    className="g-3 align-items-end mb-2"
                                )
                            ], width=6),
                            
                            dbc.Col([
                                html.Label("Map Mode:", className="fw-bold me-2"),
                                dcc.RadioItems(
                                    id="map-mode",
                                    options=[
                                        {"label": " Sequences", "value": "sequences"},
                                        {"label": " Coverage", "value": "coverage"}, 
                                        {"label": "Epidemiology", "value": "epidemiology"},
                                    ],
                                    value="sequences",
                                    inline=True,
                                    inputStyle={"marginRight": "6px", "marginLeft": "12px"}
                                ),
                            ], width="auto"),
                        ], className="g-3 align-items-end mb-2"),
                    ])
                ], className="mb-3 shadow-sm")
            ], width=12)
        ]),

        # MAIN CONTENT
        dbc.Row([
            # Indicators
            build_indicators("HBV"),
            
            # Map Section
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.H4(id="map-title-main", className="mb-1"),
                            html.Small(id="map-title-sub", className="text-muted")
                        ], className="mb-3"),
                        
                        dbc.Row([
                            dbc.Col([
                                html.Label("IHME Metric:", className="fw-bold me-2"),
                                dcc.Dropdown(
                                    id="ihme-metric-type",
                                    options=[
                                        {"label": "Deaths", "value": "Deaths|Number"},
                                        {"label": "Incidence", "value": "Incidence|Number"},
                                        {"label": "Prevalence", "value": "Prevalence|Number"},
                                    ],
                                    value="Prevalence|Number",
                                    clearable=False,
                                    style={"width": "200px"}
                                )
                            ], width="auto", id="epidemiology-controls"),
                        ], className="mb-3 g-3"),
                        
                        dcc.Loading(
                            dcc.Graph(
                                id="genotype-map", 
                                config={'displayModeBar': True, 'displaylogo': False},
                                style={"height": "100%"}
                            ),
                            type="circle"
                        )
                    ])
                ], className="h-100")
            ], width=9)
        ], className="mb-4"),
        
        # FIRST ROW: Burden Forecast and Priority Ranking
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5("Burden Forecast with Projections", className="mb-3"),
                        dcc.Loading(
                            dcc.Graph(id="forecast-chart"),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5("Sequencing Priority Ranking", className="mb-3"),
                        dbc.Row([
                            dbc.Col([
                                dbc.Button(
                                    "Download Priority Table (CSV)",
                                    id="priority-download-btn",
                                    color="secondary",
                                    className="ms-2"
                                ),
                                dcc.Download(id="priority-download")
                            ], width="auto"),
                        ], className="g-3 mb-2"),
                        dcc.Loading(
                            dcc.Graph(id="priority-ranking"),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=6)
        ], className="mb-4"),
        
        # SECOND ROW: Time Series and Country Distribution
		dbc.Row([
			dbc.Col([
				dbc.Card([
					dbc.CardBody([
						html.H5(id="line-title-main", className="mb-3"),
						dcc.Loading(
							dcc.Graph(id="line-chart", style={"height": "100%"}),
							type="circle"
						)
					])
				], className="h-100 shadow-sm")
			], width=6),
			
			dbc.Col([
				dbc.Card([
					dbc.CardBody([
						html.H5("Sequences by Country and Genotype", className="mb-3"),
						dbc.Row([
							dbc.Col([
								html.Label("Top N Countries:", className="fw-bold me-2"),
								dcc.Dropdown(
									id="top-countries-count",
									options=[{"label": str(i), "value": i} for i in [10, 15, 20, 25]],
									value=10,
									clearable=False,
									style={"width": "150px"}
								)
							], width="auto")
						], className="mb-3"),
						dcc.Loading(
							dcc.Graph(id="country-barchart", style={"height": "100%"}),
							type="circle"
						)
					])
				], className="h-100 shadow-sm")
			], width=6)
		], className="mb-4"),        
        
        # THIRD ROW: Genotype and Mutations
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5(id="bar-title-main", className="mb-3"),
                        dcc.Loading(
                            dcc.Graph(id="genotype-bar-chart", style={"height": "100%"}),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5(id="mutation-title", className="mb-3"),
                        html.Label("Filter by:", className="fw-bold mb-2"),
                        dcc.Dropdown(
                            id="mutation-filter-dropdown",
                            multi=True,
                            placeholder="Select Drug (HBV) or Gene (HCV)",
                            className="mb-3"
                        ),
                        dcc.Loading(
                            dcc.Graph(id="mutation-barplot", style={"height": "100%"}),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=6)
        ], className="mb-4"),
        
        # FOURTH ROW: FULL WIDTH Mutation Timeline
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5("Mutation Timeline", className="mb-3"),
                        dbc.Row([
                            dbc.Col([
                                html.Label("Top N Mutations:", className="fw-bold me-2"),
                                dcc.Dropdown(
                                    id="top-mutations-count",
                                    options=[{"label": str(i), "value": i} for i in [5, 10, 15, 20]],
                                    value=10,
                                    clearable=False,
                                    style={"width": "150px"}
                                )
                            ], width="auto")
                        ], className="mb-2"),
                        dcc.Loading(
                            dcc.Graph(id="mutation-timeline"),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=12)
        ], className="mb-4"),
        
        # Footer
        html.Footer([
            dbc.Container([
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.Img(src="/assets/ceri_logo.png", 
                                    className="footer-logo mx-2",
                                    style={"height": "45px", "objectFit": "contain"}),
                            html.Img(src="/assets/CRICK.png", 
                                    className="footer-logo mx-2",
                                    style={"height": "45px", "objectFit": "contain"}),
                            html.Img(src="/assets/AHRI_logo.png", 
                                    className="footer-logo mx-2",
                                    style={"height": "45px", "objectFit": "contain"})
                        ], className="d-flex justify-content-center align-items-center flex-wrap")
                    ], width=12, className="mb-3")
                ]),
                
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.P([
                                f"© {datetime.now().year} Hepatitis Virus Sequence Dashboard. ",
                                html.Span("All rights reserved.", className="text-muted")
                            ], className="mb-1"),
                            html.P("Developed by Derek Tshiabuila, Vagner Fonseca, Eduan Wilkinson, Tulio de Oliveira", 
                                  className="mb-1 text-muted"),
                            html.A("GitHub Repository", 
                                  href="https://github.com/CERI-KRISP/Hepatitis-Dashboard.git", 
                                  target="_blank", 
                                  className="text-decoration-none text-primary",
                                  style={"fontSize": "0.9rem"})
                        ], className="text-center")
                    ], width=12)
                ])
            ], fluid=True)
        ], className="bg-light py-4 mt-5 border-top", 
           style={"marginTop": "2rem !important"})
        
    ], fluid=True, style={'backgroundColor': '#f8f9fa', 'minHeight': '100vh', 'padding': '20px'})
	
# === STATE STORES ==============================================================
# — Year bounds & dropdown option lists (fast) —
layout = create_dashboard_layout()	# contains the year-slider

@callback(Output("dl-toast", "is_open"), Input("btn-download-data", "n_clicks"), prevent_initial_call=True)
def _show_toast(n): return True


@callback(
	Output("selected-years-display", "children"),
	Input("year-slider", "value"),
	prevent_initial_call=True,		 # optional
)
def show_years(value):
	return f"{value[0]} – {value[1]}"

def init_controls(virus):
	data = get_data_store()
	base = data['hbv_data'] if (virus or "HBV")=="HBV" else data['hcv_data']
	if base.empty:
		return 2000, 2023, [2000, 2023], {}, [], [], []
	y0, y1 = int(base["Year"].min()), int(base["Year"].max())
	marks = {y: (str(y) if y%5==0 else "") for y in range(y0, y1+1)}
	cont_opts = [{"label": r, "value": r} for r in sorted(base["WHO_Regions"].dropna().unique()) if r!="Unknown"]
	country_opts = [{"label": c, "value": c} for c in sorted(base["Country_standard"].dropna().unique()) if c!="Unknown"]
	geno_opts = [{"label": g, "value": g} for g in sorted(base["Genotype"].dropna().unique())]
	return y0, y1, [y0, y1], marks, cont_opts, country_opts, geno_opts

# — Filtered sequence dataframe store —
@callback(
	Output("filtered-store", "data"),
	Input("selected-virus", "data"),
	Input("year-slider", "value"),
	Input("continent-dropdown", "value"),
	Input("country-dropdown", "value"),
	Input("genotype-dropdown", "value"),
)
def compute_filtered_store(virus, years, regions, countries, genotypes):
	data = get_data_store()	 # UPDATED
	if data['hbv_data'].empty and data['hcv_data'].empty:
		return _df_to_json(pd.DataFrame())

	base = data['hbv_data'] if (virus or "HBV") == "HBV" else data['hcv_data']
	if base.empty:
		return _df_to_json(pd.DataFrame())

	# year range
	if years and len(years) == 2:
		y0, y1 = years
	else:
		y0, y1 = int(base["Year"].min()), int(base["Year"].max())

	df = base[(base["Year"] >= y0) & (base["Year"] <= y1)].copy()
	if regions:
		df = df[df["WHO_Regions"].isin(regions)]
	if countries:
		df = df[df["Country_standard"].isin(countries)]
	if genotypes:
		df = df[df["Genotype"].isin(genotypes)]

	light = df[["Country_standard", "WHO_Regions", "Year", "Genotype"]].copy()
	return _df_to_json(light)


# — Burden-adjusted coverage store —
@callback(
	Output("gap-store", "data"),
	Input("filtered-store", "data"),
	Input("selected-virus", "data"),
	Input("ihme-metric-type", "value"),
)
def compute_gap_from_filtered(filtered_json, virus, ihme_metric_choice):
	seq_df = _df_from_json(filtered_json)
	data = get_data_store()	 # UPDATED
	if seq_df.empty:
		return _df_to_json(pd.DataFrame(columns=["Country_standard","coverage_gap"]))
	gap = compute_gap_df(
		virus=(virus or "HBV"),
		filtered_seq_df=seq_df,
		ihme_df=data["ihme_df"],  # UPDATED
		selected_years=[seq_df["Year"].min(), seq_df["Year"].max()],
		who_regions=None, countries=None,	  # already applied
		ihme_metric_choice=ihme_metric_choice or BURDEN_MEASURE_FALLBACK,
		target_per_10k=5.0,
	)
	# tidy & drop Unknown
	if "Country_standard" in gap.columns:
		gap = gap[gap["Country_standard"]!="Unknown"]
	return _df_to_json(gap)

# — Latest IHME per country store —
@callback(
	Output("ihme-latest-store", "data"),
	Input("selected-virus", "data"),
	Input("ihme-metric-type", "value"),
	Input("year-slider", "value"),
	Input("continent-dropdown", "value"),
	Input("country-dropdown", "value"),
)
def compute_ihme_latest_store(virus, metric, years, regions, countries):
	data = get_data_store()	 # UPDATED
	df = ihme_latest_by_country(
		ihme_df=data["ihme_df"],  # UPDATED
		virus=(virus or "HBV"),
		measure_metric=(metric or BURDEN_MEASURE_FALLBACK),
		regions=regions, countries=countries, years=years
	)
	return _df_to_json(df)
	
# — Indicator values —
@callback(
	Output("indicator-total", "children"),
	Output("indicator-countries", "children"),
	Output("indicator-genotypes", "children"),
	Output("indicator-years", "children"),
	Input("filtered-store", "data"),
	Input("year-slider", "value"),
)
def update_indicators(filtered_json, selected_years):
	df = _df_from_json(filtered_json)
	if df is None or df.empty:
		return "0", "0", "0", "N/A"

	# Totals
	total_genomes = len(df)
	unique_countries = df.get("Country_standard", pd.Series(dtype="object")).nunique()

	# Genotypes: exclude recombinants from the count, but flag presence
	g = df.get("Genotype", pd.Series(dtype="object")).astype("string").str.strip()
	recomb_mask = g.str.contains(r"recomb", case=False, na=False)  # catches 'Recombinant', 'Recombinants', etc.
	base_genotype_count = g[~recomb_mask].dropna().nunique()
	has_recomb = bool(recomb_mask.any())

	if has_recomb:
		genotypes_text = f"{base_genotype_count} + Recombinants"
	else:
		genotypes_text = f"{base_genotype_count}"

	# Years label
	years_text = f"{selected_years[0]} - {selected_years[1]}" \
		if selected_years and len(selected_years) == 2 else "All years"

	return f"{total_genomes:,}", str(unique_countries), genotypes_text, years_text


# === FIGURE CALLBACKS ==============================================================
# - Line trend -
@callback(
	Output("line-chart", "figure"),
	Output("line-title-main", "children"),
	Input("filtered-store", "data"),
	Input("selected-virus", "data"),
)
def render_line(filtered_json, virus):
	df = _df_from_json(filtered_json)
	if df.empty:
		empty = go.Figure()
		empty.update_layout(title="No data", xaxis={"visible": False}, yaxis={"visible": False})
		return empty, "No Data"
	selected_virus = virus or "HBV"
	return make_line_trend(df, selected_virus), f"{selected_virus.upper()} Whole Genomes Per Year"


# - Genotype bar -
@callback(
	Output("genotype-bar-chart", "figure"),
	Output("bar-title-main", "children"),
	Input("filtered-store", "data"),
	Input("selected-virus", "data"),
	Input("display-mode", "value"),
)
def render_genotype_bar(filtered_json, virus, display_mode):
	df = _df_from_json(filtered_json)
	data = get_data_store()	 # UPDATED
	fig = make_genotype_bar(
		filtered_df=df,
		population_df=data["population_df"],  # UPDATED
		selected_virus=(virus or "HBV").upper(),
		display_mode=(display_mode or "raw"),
	)
	return fig, f"Total {(virus or 'HBV').upper()} Sequences by Genotype"

# - Map -
# - Map callback with three distinct modes -
@callback(
    Output("genotype-map", "figure"),
    Output("map-title-main", "children"),
    Output("map-title-sub", "children"),
    Output("epidemiology-controls", "style"),  # Show/hide IHME metric selector
    Input("filtered-store", "data"),
    Input("gap-store", "data"), 
    Input("ihme-latest-store", "data"),
    Input("selected-virus", "data"),
    Input("display-mode", "value"),         # raw / PerMillion (for sequences mode only)
    Input("map-mode", "value"),             # sequences / coverage / epidemiology
    Input("ihme-metric-type", "value"),     # For epidemiology mode
    State("continent-dropdown", "value"),
    State("country-dropdown", "value"),
)
def render_map(filtered_json, gap_json, ihme_json, virus, display_mode, map_mode, 
               ihme_metric, regions, countries):
    selected_virus = (virus or "HBV")
    filtered = _df_from_json(filtered_json)
    gap = _df_from_json(gap_json)
    ihme = _df_from_json(ihme_json)
    data = get_data_store()
    
    # Show/hide epidemiology controls based on map mode - DEFINE THIS AT THE START
    epi_controls_style = {"display": "block"} if map_mode == "epidemiology" else {"display": "none"}

    # MODE 1: Coverage map
    if map_mode == "coverage":
        coverage_for_map = gap.copy()
        if not coverage_for_map.empty and "observed_sequences" in coverage_for_map.columns and "expected_sequences" in coverage_for_map.columns:
            coverage_for_map["coverage_ratio"] = (
                coverage_for_map["observed_sequences"] /
                coverage_for_map["expected_sequences"].replace({0: np.nan})
            )
        else:
            coverage_for_map = pd.DataFrame(columns=["coverage_ratio", "Country_standard"])
        
        fig = create_coverage_map(
            coverage_for_map, data["coord_lookup"], data.get("coords"),
            virus_type=selected_virus, who_regions=regions, countries=countries
        )
        title = f"{selected_virus} Burden-adjusted Sequencing Coverage"
        subtitle = "Coverage = sequences / estimated infections"
        return fig, title, subtitle, epi_controls_style

    # MODE 2: Epidemiology map (IHME data)
    elif map_mode == "epidemiology":
        if ihme.empty:
            # Empty state for epidemiology
            fig = _empty_world("No IHME epidemiology data available for current filters")
            title = f"{selected_virus} Epidemiology Map"
            subtitle = f"({ihme_metric}) - No data"
            return fig, title, subtitle, epi_controls_style
        
        # Use IHME data directly for epidemiology mode
        fig = create_world_map(
            country_data=ihme,
            country_genotype_counts=pd.DataFrame(),  # No genotype markers in epidemiology mode
            coord_lookup=data["coord_lookup"],
            virus_type=selected_virus,
            display_mode="ihme"  # Force IHME display mode
        )
        
        # Parse metric for title
        try:
            measure, metric = (ihme_metric or "Prevalence|Number").split("|")
            metric_display = f"{measure} ({metric})"
        except:
            metric_display = ihme_metric or "Prevalence (Number)"
            
        title = f"{selected_virus} Epidemiology Map"
        subtitle = f"{metric_display}"
        return fig, title, subtitle, epi_controls_style

    # MODE 3: Sequences map (default)
    else:
        # Prepare data for the sequence map
        if display_mode == "PerMillion":
            # Calculate per million values
            if not filtered.empty and "Population" in data and not data["population_df"].empty:
                pop_data = data["population_df"]
                country_data = filtered.groupby("Country_standard").size().reset_index(name="count")
                country_data = country_data.merge(
                    pop_data[["Country_standard", "Population"]].drop_duplicates(),
                    on="Country_standard",
                    how="left"
                )
                country_data["Metric_raw"] = (country_data["count"] / country_data["Population"]) * 1_000_000
            else:
                country_data = filtered.groupby("Country_standard").size().reset_index(name="Metric_raw")
        else:
            # Raw counts
            country_data = filtered.groupby("Country_standard").size().reset_index(name="Metric_raw")
        
        # Prepare genotype counts for sequence map
        country_genotype_counts = filtered.groupby(["Country_standard", "Genotype"]).size().reset_index(name="Count")
        
        # Create the sequence map
        fig = create_world_map(
            country_data,
            country_genotype_counts,
            data["coord_lookup"],
            virus_type=selected_virus,
            display_mode=display_mode or "raw"
        )
        
        # Set titles based on display mode
        if display_mode == "PerMillion":
            title = f"{selected_virus} Whole Genome Sequence Map"
            subtitle = "Sequences per million population"
        else:
            title = f"{selected_virus} Whole Genome Sequence Map" 
            subtitle = "Sequences (count)"
        
        return fig, title, subtitle, epi_controls_style
	
# - Mutation barplot -
@callback(
	Output("mutation-barplot", "figure"),
	Output("mutation-title", "children"),
	Input("filtered-store", "data"),
	Input("mutation-filter-dropdown", "value"),
	Input("selected-virus", "data"),
)
def render_mutation_bar(filtered_json, selected_filter, virus):
	seq_df = _df_from_json(filtered_json)
	selected_virus = (virus or "HBV").upper()
	data = get_data_store()	 # UPDATED

	# Total sequences (denominator)
	total_sequences = len(seq_df)

	# Choose mutation table + full source sequences
	if selected_virus == "HBV":
		mut = data["hbv_mut"].copy()
		seq_source = data["hbv_data"]
		facet_col = "drug"
	else:
		mut = data["hcv_mut"].copy()
		seq_source = data["hcv_data"]
		facet_col = "gene"

	# Enrich mutation rows with Country/Region/Year/Genotype
	enriched = _enrich_mutation_df(mut, seq_source)

	# Base span (ignore year filter): match Region/Country/Genotype only
	if not seq_df.empty:
		keys_base = ["Country_standard", "WHO_Regions", "Genotype"]
		base_for_span = enriched.merge(seq_df[keys_base].drop_duplicates(), on=keys_base, how="inner")
	else:
		base_for_span = enriched

	# Actual overlap for current selection (includes Year)
	if not seq_df.empty:
		keys_full = ["Country_standard", "WHO_Regions", "Year", "Genotype"]
		overlap = enriched.merge(seq_df[keys_full].drop_duplicates(), on=keys_full, how="inner")
	else:
		overlap = enriched

	# Apply filter by drug/gene
	if selected_filter:
		vals = [str(s).strip().lower() for s in selected_filter]
		overlap = overlap[overlap[facet_col].astype(str).str.strip().str.lower().isin(vals)]

	# Compute years_range for the *current* selection
	if not seq_df.empty and "Year" in seq_df.columns and seq_df["Year"].notna().any():
		y0, y1 = int(np.nanmin(seq_df["Year"])), int(np.nanmax(seq_df["Year"]))
		years_range = [y0, y1]
	else:
		years_range = None

	# Compute data_span from base_for_span (for "does years actually narrow the span?")
	if not base_for_span.empty and base_for_span["Year"].notna().any():
		ymin, ymax = int(np.nanmin(base_for_span["Year"])), int(np.nanmax(base_for_span["Year"]))
		data_span = [ymin, ymax]
	else:
		data_span = None

	# Other filters active? (is current selection a strict subset of all seqs)
	try:
		all_ids = set(seq_source["ID"]) if "ID" in seq_source.columns else set()
		sel_ids = set(seq_df["ID"]) if "ID" in seq_df.columns else set()
		other_active = bool(all_ids and len(sel_ids) < len(all_ids))
	except Exception:
		other_active = False

	# Build the figure + title in one place
	fig, title = make_mutation_bar(
		mutation_df=overlap,
		total_sequences=total_sequences,
		selected_virus=selected_virus,
		selected_filter=selected_filter,
		years_range=years_range,
		data_span=data_span,
		other_filters_active=other_active,
	)

	return fig, title

# Forecast Chart Callback
@callback(
    Output("forecast-chart", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_forecast_chart(virus, regions, countries):
    data = get_data_store()
    return create_forecast_chart(
        data["ihme_df"],
        virus or "HBV", 
        regions, 
        countries
    )
    
@callback(
    Output("mutation-timeline", "figure"),
    Input("selected-virus", "data"),
    Input("filtered-store", "data"),
    Input("top-mutations-count", "value"),
)
def update_mutation_timeline(virus, filtered_json, top_n):
    df = _df_from_json(filtered_json)
    selected_virus = virus or "HBV"
    data = get_data_store()
    
    # Get mutation data based on virus type
    if selected_virus.upper() == "HBV":
        mutation_data = data["hbv_mut"]
    else:
        mutation_data = data["hcv_mut"]
    
    return create_mutation_timeline(mutation_data, df, selected_virus, top_n)

# Country barchart
@callback(
    Output("country-barchart", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("top-countries-count", "value"),  # ADD THIS
)
def update_country_stacked_bar_callback(filtered_json, virus, regions, countries, top_n):
    df = _df_from_json(filtered_json)
    return create_country_stacked_bar(df, virus or "HBV", regions, countries, top_n=top_n or 15)

# Update priority ranking to be more responsive
@callback(
    Output("priority-ranking", "figure"),
    Output("priority-data-store", "data"),
    Input("gap-store", "data"),
    Input("selected-virus", "data"),
)
def update_priority_ranking_responsive(gap_json, virus):
    gap_df = _df_from_json(gap_json)
    data = get_data_store()
    
    # Use default weights instead of user inputs
    weights = {
        "burden": 0.4,
        "coverage_gap": 0.3,
        "population": 0.2,
        "neighbor_sequencing": 0.1
    }
    
    fig, priority_df = create_priority_calculator(
        gap_df, data["ihme_df"], virus or "HBV", weights
    )
    
    return fig, _df_to_json(priority_df)

# Priority Download Callback
@callback(
    Output("priority-download", "data"),
    Input("priority-download-btn", "n_clicks"),
    State("priority-data-store", "data"),
    State("selected-virus", "data"),
    prevent_initial_call=True,
)
def download_priority_table(n_clicks, priority_json, virus):
    if not n_clicks:
        return dash.no_update
    
    priority_df = _df_from_json(priority_json)
    if priority_df.empty:
        return dcc.send_string("No priority data available", "priority_data_empty.txt")
    
    # Create comprehensive CSV with priority data
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{virus}_sequencing_priority_ranking_{timestamp}.csv"
    
    return dcc.send_data_frame(priority_df.to_csv, filename, index=False)
	
@callback(
	Output("year-slider", "min"),
	Output("year-slider", "max"),
	Output("year-slider", "value"),
	Output("year-slider", "marks"),
	Output('continent-dropdown', 'options'),
	Output('country-dropdown', 'options'),
	Output('genotype-dropdown', 'options'),
	Input("selected-virus", "data"),
)
def init_controls(virus):
	data = get_data_store()	 # UPDATED
	base = data['hbv_data'] if (virus or "HBV")=="HBV" else data['hcv_data']
	if base.empty:
		return 2000, 2023, [2000, 2023], {}, [], [], []
	y0, y1 = int(base["Year"].min()), int(base["Year"].max())
	marks = {y: (str(y) if y % 5 == 0 else "") for y in range(y0, y1 + 1)}
	cont_opts = [{"label": r, "value": r} for r in sorted(base["WHO_Regions"].dropna().unique()) if r!="Unknown"]
	country_opts = [{"label": c, "value": c} for c in sorted(base["Country_standard"].dropna().unique()) if c!="Unknown"]
	geno_opts = [{"label": g, "value": g} for g in sorted(base["Genotype"].dropna().unique())]
	return y0, y1, [y0, y1], marks, cont_opts, country_opts, geno_opts

# === ACTIONS AND DOWNLOADS ==============================================================
# - Main data download with Taxa -
def _pick_mode(series: pd.Series):
	"""Pick a single representative value (mode → first non-null fallback)."""
	m = series.mode(dropna=True)
	if not m.empty:
		return m.iloc[0]
	s = series.dropna()
	return s.iloc[0] if not s.empty else None

def _build_keys(df: pd.DataFrame, is_main: bool) -> pd.DataFrame:
	"""
	Create harmonized join keys without renaming over existing columns.
	Avoids duplicate-column-name collisions and keeps groupby 1-D.
	"""
	# Drop duplicated column NAMES if any (keeps first occurrence)
	df = df.loc[:, ~df.columns.duplicated()].copy()

	# Country key (prefer Country_standard, fallback to Country)
	if "Country_standard" in df.columns and "Country" in df.columns:
		df["Country_key"] = df["Country_standard"].fillna(df["Country"])
	elif "Country_standard" in df.columns:
		df["Country_key"] = df["Country_standard"]
	elif "Country" in df.columns:
		df["Country_key"] = df["Country"]
	else:
		df["Country_key"] = pd.NA

	# Year key (accept Year or Date; coerce to numeric)
	year_num = None
	if "Year" in df.columns:
		year_num = pd.to_numeric(df["Year"], errors="coerce")
	if "Date" in df.columns:
		date_num = pd.to_numeric(df["Date"], errors="coerce")
		year_num = year_num.fillna(date_num) if year_num is not None else date_num
	df["Year_key"] = year_num if year_num is not None else pd.NA

	# Genotype key
	df["Genotype_key"] = df["Genotype"] if "Genotype" in df.columns else pd.NA

	# Normalize types for robust grouping/merging
	for c in ["Country_key", "Genotype_key"]:
		df[c] = df[c].astype("string")

	# Ensure expected columns exist on main DF (for downstream CSV)
	if is_main:
		for c in ["Country_standard", "WHO_Regions", "Year", "Genotype"]:
			if c not in df.columns:
				df[c] = pd.NA

	return df

# ----------------------------------------------------------------------
# Callback
# ----------------------------------------------------------------------
@callback(
	Output("download-data", "data"),
	Input("btn-download-data", "n_clicks"),
	State("filtered-store", "data"),
	State("selected-virus", "data"),
	State("year-slider", "value"),
	State("continent-dropdown", "value"),
	State("country-dropdown", "value"),
	State("genotype-dropdown", "value"),
	prevent_initial_call=True,
)
def download_main_data_with_taxa(n_clicks, filtered_json, virus, years, regions, countries, genotypes):
	# Only act on actual clicks
	if not n_clicks:
		raise PreventUpdate

	# Parse filtered data safely
	try:
		filtered_df = _df_from_json(filtered_json)
	except Exception:
		filtered_df = pd.DataFrame()

	selected_virus = (virus or "HBV").upper()
	data = get_data_store()	 # UPDATED

	# If empty → still return a small CSV so the click always downloads something
	if filtered_df.empty:
		empty_df = pd.DataFrame(columns=["Country_standard", "WHO_Regions", "Year", "Genotype", "Taxa"])
		return dcc.send_data_frame(
			empty_df.to_csv,
			f"{selected_virus}_no_data_available.csv",
			index=False
		)

	# Prepare main DF with join keys
	download_df = _build_keys(filtered_df.copy(), is_main=True)

	# Get summary data by virus
	summary_key = "hbv_summary_raw" if selected_virus == "HBV" else "hcv_summary_raw"
	summary_data = data.get(summary_key, pd.DataFrame())  # UPDATED

	if not summary_data.empty:
		summary_clean = _build_keys(summary_data.copy(), is_main=False)

		# Try progressively less-specific keys
		merge_strategies = [
			(["Genotype_key", "Country_key", "Year_key"], "exact_match"),
			(["Genotype_key", "Country_key"], "genotype_country"),
			(["Genotype_key"], "genotype_only"),
			(["Country_key"], "country_only"),
		]

		merged_successfully = False
		base_len = len(download_df)

		for merge_cols, _strategy in merge_strategies:
			# Collapse to ONE Taxa per key to avoid many-to-many merges
			taxa_map = (
				summary_clean
				.dropna(subset=["Taxa"])
				.groupby(merge_cols, dropna=False)["Taxa"]
				.agg(_pick_mode)
				.reset_index()
			)

			if taxa_map.empty or taxa_map["Taxa"].isna().all():
				continue

			try:
				merged = download_df.merge(taxa_map, on=merge_cols, how="left")
			except Exception:
				continue

			# Guard against accidental row multiplication
			if len(merged) > base_len * 2:
				# Skip this strategy if it explodes rows
				continue

			download_df = merged
			merged_successfully = True
			break

		if not merged_successfully:
			download_df["Taxa"] = "No Taxa match found"
	else:
		download_df["Taxa"] = "Summary data not available"

	# Put Taxa first if present
	if "Taxa" in download_df.columns:
		ordered = ["Taxa"] + [c for c in download_df.columns if c != "Taxa"]
		download_df = download_df.reindex(columns=ordered)

	# Construct filename
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	filename = f"{selected_virus}_sequence_data_with_taxa_{timestamp}.csv"

	# Return CSV
	return dcc.send_data_frame(download_df.to_csv, filename, index=False)
	
def merge_taxa_information(main_df, summary_df, virus_type):
	"""
	Merge Taxa information from summary data into main dataframe
	"""
	if summary_df.empty or main_df.empty:
		main_df["Taxa"] = "No summary data available" if summary_df.empty else "No matches found"
		return main_df
	
	summary_clean = summary_df.rename(columns={
		"Country": "Country_standard",
		"Date": "Year"
	}).copy()
	
	summary_clean["Year"] = pd.to_numeric(summary_clean["Year"], errors="coerce")
	
	# Try multiple merge strategies
	strategies = [
		(["Genotype", "Country_standard", "Year"], "exact_match"),
		(["Genotype", "Country_standard"], "genotype_country"), 
		(["Genotype"], "genotype_only"),
		(["Country_standard"], "country_only"),
	]
	
	for merge_cols, strategy in strategies:
		if all(col in main_df.columns and col in summary_clean.columns for col in merge_cols):
			try:
				merged = main_df.merge(
					summary_clean[merge_cols + ["Taxa"]].drop_duplicates(),
					on=merge_cols,
					how="left"
				)
				if merged["Taxa"].notna().any():
					print(f"Taxa merge successful with {strategy}")
					return merged
			except Exception as e:
				print(f"Taxa merge failed with {strategy}: {e}")
				continue
	
	main_df["Taxa"] = "No match found"
	return main_df
	
@callback(
	Output("download-data", "data", allow_duplicate=True),
	Input("btn-download", "n_clicks"),
	State("filtered-store", "data"),
	State("selected-virus", "data"),
	prevent_initial_call=True,
)
def download_main_data_simple(n_clicks, filtered_json, virus):
	if not n_clicks or n_clicks == 0:
		return dash.no_update
	
	filtered_df = _df_from_json(filtered_json)
	selected_virus = virus or "HBV"
	summary_key = "hbv_summary_raw" if selected_virus == "HBV" else "hcv_summary_raw"
	summary_data = data_store.get(summary_key, pd.DataFrame())
	
	# Merge Taxa information
	download_df = merge_taxa_information(filtered_df, summary_data, selected_virus)
	
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	filename = f"{selected_virus}_sequence_data_with_taxa_{timestamp}.csv"
	
	return dcc.send_data_frame(download_df.to_csv, filename, index=False)
	


# === SMALL UI TOGGLES ==============================================================
@callback(
	Output("selected-virus", "data"),
	Output("btn-hbv", "outline"),
	Output("btn-hcv", "outline"),
	Input("btn-hbv", "n_clicks"),
	Input("btn-hcv", "n_clicks"),
	prevent_initial_call=True
)

def update_virus(btn_hbv_clicks, btn_hcv_clicks):
	ctx = callback_context
	if not ctx.triggered:
		return dash.no_update, dash.no_update, dash.no_update
	
	triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
	
	if triggered_id == "btn-hbv":
		return "HBV", False, True
	else:
		return "HCV", True, False


# - Mutation section toggle -
@callback(
	Output("hbv-mutation-section", "style"),
	Output("hcv-mutation-section", "style"),
	Input("selected-virus", "data")
)
def toggle_mutation_sections(selected_virus):
	selected_virus = selected_virus or "HBV"
	if selected_virus == "HBV":
		return {"display": "block"}, {"display": "none"}
	elif selected_virus == "HCV":
		return {"display": "none"}, {"display": "block"}
	else:
		return {"display": "none"}, {"display": "none"}
		
# - Mutation filter options -
# DROP-IN REPLACEMENT
@callback(
	Output("mutation-filter-dropdown", "options"),
	Input("selected-virus", "data"),
	prevent_initial_call=False,
)
def update_mutation_filter_options(virus):
	data = get_data_store() or {}  # UPDATED
	selected = (virus or "HBV").upper()

	# Pick the right mutations table
	mut = data.get("hbv_mut") if selected == "HBV" else data.get("hcv_mut")
	if mut is None or len(mut) == 0:
		return []  # nothing to show yet

	# Column names differ: HBV uses 'drug', HCV uses 'gene'
	col = "drug" if selected == "HBV" else "gene"
	if col not in mut.columns:
		# Be defensive: try to discover a plausible column
		for candidate in ["drug", "gene", "Drug", "Gene"]:
			if candidate in mut.columns:
				col = candidate
				break
		else:
			return []

	opts = sorted(pd.Series(mut[col]).dropna().astype(str).unique())
	return [{"label": v, "value": v} for v in opts]

		
# - Year reset -button
@callback(
	Output("year-slider", "value", allow_duplicate=True),
	Input("btn-reset-time", "n_clicks"),
	State("selected-virus", "data"),
	prevent_initial_call=True
)
def reset_time_range(n_clicks, selected_virus):
	if n_clicks is None or n_clicks == 0:
		raise PreventUpdate
		
	data = get_data_store()
	df = data['hbv_data'] if selected_virus == "HBV" else data['hcv_data']
	min_year = int(df["Year"].min())
	max_year = int(df["Year"].max())
	
	return [min_year, max_year]
