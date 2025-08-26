import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from datetime import datetime
import pycountry
from functools import lru_cache
import numpy as np
from dash import callback_context
from concurrent.futures import ThreadPoolExecutor
import re
import json


# Pre-compile regex patterns
ACCESSION_PATTERN = re.compile(r"([A-Z0-9]+\.\d+)")
WHITESPACE_PATTERN = re.compile(r"\s+")

# Pre-defined constants and mappings
SPECIAL_CASES = {
	'USA': 'United States',
	'US': 'United States',
	'CapeVerde': 'Cape Verde',
	'Cabo Verde': 'Cape Verde',
	'Democratic Republic of the Congo': 'Congo [DRC]',
	'Congo, Dem. Rep.': 'Congo [DRC]',
	'Congo, Rep.': 'Congo [Republic]',
	'Congo': 'Congo [Republic]',
	'Federal Republic of Nigeria': 'Nigeria',
	'Republic of Niger': 'Niger',
	'Republic of Serbia': 'Serbia',
	'Myanmar': 'Myanmar [Burma]',
	'VietNam': 'Viet Nam',
	'Vietnam': 'Viet Nam',
	'Iran, Islamic Republic of': 'Iran',
	'Iran (Islamic Republic of)': 'Iran',
	'Iran, Islamic Rep.': 'Iran',
	'Iran': 'Iran',
	'Korea, Republic of': 'South Korea',
	'South Korea': 'South Korea',
	'North Korea': 'North Korea',
	"Lao People's Democratic Republic": 'Laos',
	'Lao PDR': 'Laos',
	'Laos': 'Laos',
	'Côte d\'Ivoire': 'Ivory Coast',
	'Cote d\'Ivoire': 'Ivory Coast',
	'Bolivia (Plurinational State of)': 'Bolivia, Plurinational State of',
	'Bolivia': 'Bolivia, Plurinational State of',
	'Venezuela (Bolivarian Republic of)': 'Venezuela, Bolivarian Republic of',
	'Venezuela': 'Venezuela, Bolivarian Republic of',
	'Tanzania, United Republic of': 'Tanzania',
	'United Republic of Tanzania': 'Tanzania',
	'Tanzania': 'Tanzania',
	'Swaziland': 'Eswatini',
	'Turkey': 'Turkey',
	'St. Lucia': 'Saint Lucia',
	'St. Kitts and Nevis': 'Saint Kitts and Nevis',
	'St. Vincent and the Grenadines': 'Saint Vincent and the Grenadines',
	'Hong Kong': 'China',
    'Hong-Kong': 'China',
    'HongKong': 'China',
    'Hong Kong SAR': 'China',
    'Hong Kong SAR, China': 'China',
    'Hong Kong, SAR China': 'China',
    'Hong Kong Special Administrative Region': 'China',
    'Hong Kong S.A.R.': 'China',
    'HKSAR': 'China',
    'HK': 'China',     # sometimes appears as a short code in metadata
    'HKG': 'China',    # ISO3 code occasionally found in CSVs
	'Venezuela (Bolivarian Republic of)': 'Venezuela',
	'Democratic People\'s Republic of Korea': 'North Korea',
	'Republic of Korea': 'South Korea',
	'Taiwain, Province of China': 'Taiwan, Province of China',
	'Taiwan (Province of China)': 'Taiwan, Province of China',
	'Taiwan': 'Taiwan, Province of China',
    'Micronesia (Federated States of)': 'Micronesia, Federated States of',
    'Micronesia, Federated States of': 'Micronesia, Federated States of',
    'Micronesia': 'Micronesia, Federated States of',
    'Namibia': 'Namibia',  # This should already map to itself
    'Venezuela, Bolivarian Republic of': 'Venezuela',
    'Venezuela (Bolivarian Republic of)': 'Venezuela',
	'United States Virgin Islands': 'U.S. Virgin Islands',
	'Sao Tome and Principe': 'Sao Tome and Principe',
	'United States Virgin Islands': 'U.S. Virgin Islands',
	'Venezuela': 'Venezuela',
	'North Macedonia': 'Macedonia [FYROM]',
    'Macedonia': 'Macedonia [FYROM]',
    'Republic of North Macedonia': 'Macedonia [FYROM]',
	'Uknown': None,
	'unknown': None,
    'Not provided': None,
    'not provided': None,
}

ENCODINGS = ['utf-8-sig', 'utf-8', 'latin1', 'iso-8859-1', 'cp1252']

@lru_cache(maxsize=5000)
def standardize_country_cached(name):
	if not isinstance(name, str) or not name.strip():
		return None
	
	name = name.strip()
	
	# First use special case mapping
	if name in SPECIAL_CASES:
		return SPECIAL_CASES[name]
	
	# Then try fuzzy match
	try:
		match = pycountry.countries.search_fuzzy(name)
		if match:
			return match[0].name
	except:
		pass
	
	# If still no match, return original name
	return name

def load_with_encoding(file_path, sep='\t'):
	for encoding in ENCODINGS:
		try:
			return pd.read_csv(file_path, sep=sep, encoding=encoding, 
							 on_bad_lines='skip', low_memory=False)
		except Exception:
			continue
	raise ValueError(f"Could not read {file_path} with provided encodings.")

def optimize_dataframe_dtypes(df):
	"""Optimize dataframe memory usage by converting to appropriate dtypes"""
	for col in df.columns:
		if df[col].dtype == 'object':
			# Convert to category if low cardinality
			if df[col].nunique() / len(df[col]) < 0.5:
				df[col] = df[col].astype('category')
		elif df[col].dtype in ['int64', 'float64']:
			# Downcast numeric types
			if pd.api.types.is_integer_dtype(df[col]):
				df[col] = pd.to_numeric(df[col], downcast='integer')
			else:
				df[col] = pd.to_numeric(df[col], downcast='float')
	return df

def merge_population_nearest(counts_df: pd.DataFrame, pop_df: pd.DataFrame, tol_years: int = 3) -> pd.DataFrame:
	counts = counts_df.copy()
	pop = pop_df.copy()

	# Precompute standardized country names
	for df in (counts, pop):
		df["Country_standard"] = df["Country_standard"].astype(str)
		df["Country_standard"] = df["Country_standard"].str.strip()
		df["Country_standard"] = df["Country_standard"].str.replace(WHITESPACE_PATTERN, " ", regex=True)
		df["Year"] = pd.to_numeric(df["Year"], errors="coerce")

	pop["Population"] = pd.to_numeric(pop["Population"], errors="coerce")

	counts = counts.dropna(subset=["Country_standard", "Year"])
	pop = pop.dropna(subset=["Country_standard", "Year", "Population"])

	counts["Year"] = counts["Year"].astype("int32")
	pop["Year"] = pop["Year"].astype("int32")

	out = []
	# Use vectorized operations where possible
	for country, g in counts.groupby("Country_standard", sort=False, observed=True):
		gp = pop[pop["Country_standard"] == country]
		if gp.empty:
			g["Population"] = np.nan
			out.append(g)
			continue
		
		# Sort once and use merge_asof
		g_sorted = g.sort_values("Year", kind='stable')
		gp_sorted = gp[["Year", "Population"]].sort_values("Year", kind='stable')
		merged = pd.merge_asof(g_sorted, gp_sorted, on="Year", direction="nearest", tolerance=tol_years)
		out.append(merged)

	return pd.concat(out, ignore_index=True)

def merge_population_nearest_two_pass(counts_df, pop_df, tol_years_first=3, tol_years_wide=50):
	first = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_first)
	if first["Population"].isna().any():
		widened = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_wide)
		first["Population"] = first["Population"].fillna(widened["Population"])
	return first

def load_recombinant_ids(file_path):
	"""Load recombinant IDs efficiently"""
	with open(file_path) as f:
		return set(line.strip() for line in f)

def extract_accession_id(taxa_series):
	"""Vectorized accession ID extraction"""
	return taxa_series.str.extract(ACCESSION_PATTERN, expand=False)

def process_genotype_column(df, recombs_set):
	"""Efficiently process genotype column"""
	mask = df['Taxa'].isin(recombs_set)
	df.loc[mask, 'Genotype'] = 'Recombinant'
	return df

def merge_who_region(df, coords_df):
	"""Optimized WHO region merging"""
	# Create lookup dictionary for faster merging
	region_lookup = coords_df.set_index('Country_standard')['WHO'].rename('WHO_Regions')
	df = df.join(region_lookup, on='Country_standard', how='left')
	df["WHO_Regions"] = df["WHO_Regions"].fillna("Unknown").astype('category')
	return df

def standardize_country_column(df, column_name):
	"""Standardize country column in a dataframe"""
	df['Country_standard'] = df[column_name].apply(standardize_country_cached)
	return df

def get_countries_missing_who_region(data_store):
    """
    Returns a sorted list of all unique countries with a missing WHO region
    across all loaded datasets (HBV, HCV, IHME).
    
    Parameters:
    data_store (dict): The dictionary containing all preprocessed dataframes
    
    Returns:
    list: Sorted list of country names with missing WHO regions
    """
    hbv_data = data_store['hbv_data']
    hcv_data = data_store['hcv_data']
    ihme_df = data_store['ihme_df']
    
    # Get unique unknown countries from each dataset
    hbv_unknown = set(hbv_data[hbv_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    hcv_unknown = set(hcv_data[hcv_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    ihme_unknown = set(ihme_df[ihme_df["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    
    # Combine all
    all_unknown = hbv_unknown | hcv_unknown | ihme_unknown
    
    return sorted(all_unknown)

# Usage example after your data is loaded:
# missing_countries = get_countries_missing_who_region(data_store)
# print("Countries needing WHO region mapping:", missing_countries)

def load_and_preprocess_data():
	print("📥 Loading data...")
	
	# Load data sequentially for better error handling
	hbv_data = load_with_encoding("HBV_summary.tsv").dropna()
	hcv_data = load_with_encoding("HCV_summary.tsv").dropna()
	coords = load_with_encoding("WHO_regions_countries_coordinates.txt").dropna()

	# Load recombinant IDs
	hbv_recombs = load_recombinant_ids("hbv_recombs.txt")
	hcv_recombs = load_recombinant_ids("hcv_recombs.txt")

	# Process genotype columns efficiently
	hbv_data = process_genotype_column(hbv_data, hbv_recombs)
	hcv_data = process_genotype_column(hcv_data, hcv_recombs)

	# Load additional data
	population_df = pd.read_csv("population_by_country_year.csv")
	ihme_df = pd.read_csv("IHME-GBD_2021_DATA-9e7ec2c0-1.csv")
	hbv_cov = pd.read_csv("hbv_burden_adjusted_coverage.csv")
	hcv_cov = pd.read_csv("hcv_burden_adjusted_coverage.csv")

	print("🧹 Standardizing country names...")
	
	# Standardize country names for each dataframe
	hbv_data = standardize_country_column(hbv_data, 'Country')
	hcv_data = standardize_country_column(hcv_data, 'Country')
	coords = standardize_country_column(coords, 'name')
	population_df = standardize_country_column(population_df, 'Country_standard')
	ihme_df = standardize_country_column(ihme_df, 'location')
	hbv_cov = standardize_country_column(hbv_cov, 'Country_standard')
	hcv_cov = standardize_country_column(hcv_cov, 'Country_standard')

	# Sanitize strings efficiently
	for df in [hbv_data, hcv_data, coords]:
		df["Country_standard"] = df["Country_standard"].astype(str).str.strip()

	# Drop any coords with missing standardization
	coords = coords.dropna(subset=["Country_standard"])

	print("📌 Merging WHO Regions...")
	### Diagnose why countries missing
	print("📌 Merging WHO Regions...")

	# DEBUG: Check the state of the data BEFORE merging
	print("\n🔍 DEBUG: Checking a few key countries in each dataframe:")
	
	# Check what 'Namibia' looks like in each DataFrame after standardization
	namibia_in_coords = coords[coords['Country_standard'].str.contains('Namibia', na=False, case=False)]
	namibia_in_hbv = hbv_data[hbv_data['Country_standard'].str.contains('Namibia', na=False, case=False)]
	namibia_in_hcv = hcv_data[hcv_data['Country_standard'].str.contains('Namibia', na=False, case=False)]
	
	print("'Namibia' in coords:", namibia_in_coords['Country_standard'].unique() if not namibia_in_coords.empty else "NOT FOUND")
	print("'Namibia' in hbv_data:", namibia_in_hbv['Country_standard'].unique() if not namibia_in_hbv.empty else "NOT FOUND")
	print("'Namibia' in hcv_data:", namibia_in_hcv['Country_standard'].unique() if not namibia_in_hcv.empty else "NOT FOUND")
		
	# Check Micronesia
	micronesia_in_coords = coords[coords['Country_standard'].str.contains('Micronesia', na=False, case=False)]
	micronesia_in_hbv = hbv_data[hbv_data['Country_standard'].str.contains('Micronesia', na=False, case=False)]
	micronesia_in_hcv = hcv_data[hcv_data['Country_standard'].str.contains('Micronesia', na=False, case=False)]
	
	print("'Micronesia' in coords:", micronesia_in_coords['Country_standard'].unique() if not micronesia_in_coords.empty else "NOT FOUND")
	print("'Micronesia' in hbv_data:", micronesia_in_hbv['Country_standard'].unique() if not micronesia_in_hbv.empty else "NOT FOUND")
	print("'Micronesia' in hcv_data:", micronesia_in_hcv['Country_standard'].unique() if not micronesia_in_hcv.empty else "NOT FOUND")
	
	# Now proceed with the merge
	hbv_data = merge_who_region(hbv_data, coords)
	hcv_data = merge_who_region(hcv_data, coords)
	ihme_df = merge_who_region(ihme_df, coords)

	# --- DIAGNOSTICS: Find countries with missing WHO regions ---
	print("\n🔍 Checking for countries with missing WHO regions...")
	
	# Get unique unknown countries from each dataset
	hbv_unknown_countries = hbv_data[hbv_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique()
	hcv_unknown_countries = hcv_data[hcv_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique()
	ihme_unknown_countries = ihme_df[ihme_df["WHO_Regions"] == "Unknown"]["Country_standard"].unique()
	
	# Combine and get the unique set
	all_unknown_countries = set(hbv_unknown_countries) | set(hcv_unknown_countries) | set(ihme_unknown_countries)
	
	# Print the results
	print(f"⚠️  Found {len(all_unknown_countries)} unique countries with 'Unknown' WHO region:")
	for country in sorted(all_unknown_countries):
		print(f"   - {country}")
	
	# Also check if these countries exist in the coords dataframe but lack a WHO region
	print("\n🔍 Checking if these countries exist in coordinates file but lack a WHO region:")
	for country in sorted(all_unknown_countries):
		coords_match = coords[coords['Country_standard'] == country]
		if not coords_match.empty:
			who_region = coords_match['WHO'].iloc[0] if 'WHO' in coords_match.columns else 'N/A'
			print(f"   - {country}: Present in coords file, WHO value = '{who_region}'")
	
	print("📆 Processing dates...")
	for df in [hbv_data, hcv_data]:
		df["Year"] = pd.to_numeric(df["Date"], errors='coerce').fillna(0).astype('int32')
		df.dropna(subset=["Year"], inplace=True)
		# Keep Date column if needed elsewhere, otherwise remove it
		if 'Date' in df.columns:
			df.drop(columns=["Date"], inplace=True)

	# Fix population years and drop invalid rows
	population_df["Year"] = pd.to_numeric(population_df["Year"], errors="coerce").astype("Int32")
	population_df["Population"] = pd.to_numeric(population_df["Population"], errors="coerce")
	population_df = (
		population_df
		.dropna(subset=["Year", "Population"])
		.sort_values(["Country_standard", "Year"])
		.drop_duplicates(["Country_standard", "Year"], keep="last")
	)
	population_df.loc[population_df["Population"] <= 0, "Population"] = np.nan

	print("📊 Precomputing groupings...")
	# Use more efficient grouping
	grouping_cols = ["Year", "Genotype", "WHO_Regions", "Country_standard"]
	
	hbv_grouped = (hbv_data.groupby(grouping_cols, observed=True)
				   .size()
				   .reset_index(name="Count"))
	
	hcv_grouped = (hcv_data.groupby(grouping_cols, observed=True)
				   .size()
				   .reset_index(name="Count"))

	# Merge with population data
	hbv_grouped = pd.merge(hbv_grouped, population_df, on=["Country_standard", "Year"], how="left")
	hcv_grouped = pd.merge(hcv_grouped, population_df, on=["Country_standard", "Year"], how="left")

	# Calculate per million
	hbv_grouped["PerMillion"] = (hbv_grouped["Count"] / hbv_grouped["Population"]) * 1_000_000
	hcv_grouped["PerMillion"] = (hcv_grouped["Count"] / hcv_grouped["Population"]) * 1_000_000

	print("📍 Creating coordinate lookup...")
	coord_lookup = (
		coords.drop_duplicates(subset='Country_standard')
			  .set_index('Country_standard')[['latitude', 'longitude']]
			  .to_dict('index')
	)
	# Add manual overrides
	coord_lookup.update({
		"Nigeria": {"latitude": 9.0820, "longitude": 7.49508},
		"New Caledonia": {"latitude": -21.450553, "longitude": 165.505710}
	})

	print("🧬 Loading mutation data...")
	# Load mutation data
	hbv_mut = pd.read_csv("hbv_resistant_mutations.csv")
	hcv_mut = pd.read_csv("hcv_resistance_mutations.csv")

	# Extract accession IDs efficiently
	hbv_data["ID"] = extract_accession_id(hbv_data["Taxa"])
	hcv_data["ID"] = extract_accession_id(hcv_data["Taxa"])

	# Handle mutation file column names
	if "ID" not in hbv_mut.columns and "Taxa" in hbv_mut.columns:
		hbv_mut["ID"] = extract_accession_id(hbv_mut["Taxa"])
	
	if "ID" not in hcv_mut.columns and "Taxa" in hcv_mut.columns:
		hcv_mut["ID"] = extract_accession_id(hcv_mut["Taxa"])

	# Merge mutation data with main data
	merge_cols = ["ID", "Country_standard", "Year", "Genotype", "WHO_Regions"]
	
	hbv_mut = hbv_mut.merge(
		hbv_data[merge_cols],
		on="ID", how="left"
	)
	
	hcv_mut = hcv_mut.merge(
		hcv_data[merge_cols],
		on="ID", how="left"
	)

	# Drop rows with missing country/year
	hbv_mut = hbv_mut.dropna(subset=["Country_standard", "Year"])
	hcv_mut = hcv_mut.dropna(subset=["Country_standard", "Year"])

	print("✅ Finalizing mutation stats...")
	hbv_mut_stats = hbv_mut.groupby(["Country_standard", "drug", "mutation"], observed=True).size().reset_index(name="Count")
	hcv_mut_stats = hcv_mut.groupby(["Country_standard", "gene", "mutation"], observed=True).size().reset_index(name="Count")

	# Compute log coverage for plotting
	for df in (hbv_cov, hcv_cov):
		df["Coverage_log10"] = np.where(
			df["Coverage_ratio"] > 0,
			np.log10(df["Coverage_ratio"]),
			np.nan
		)

	# Optimize dataframe memory usage
	dataframes_to_optimize = [
		hbv_data, hcv_data, coords, hbv_grouped, hcv_grouped,
		population_df, ihme_df, hbv_mut, hcv_mut, hbv_mut_stats,
		hcv_mut_stats, hbv_cov, hcv_cov
	]
	
	for df in dataframes_to_optimize:
		optimize_dataframe_dtypes(df)

	print("✅ Done.")
	return {
		'hbv_data': hbv_data,
		'hcv_data': hcv_data,
		'coords': coords,
		'hbv_grouped': hbv_grouped,
		'hcv_grouped': hcv_grouped,
		'coord_lookup': coord_lookup,
		'population_df': population_df,
		'ihme_df': ihme_df,
		'hbv_mut': hbv_mut,
		'hcv_mut': hcv_mut,
		'hbv_mut_stats': hbv_mut_stats,
		'hcv_mut_stats': hcv_mut_stats,
		'cov_hbv': hbv_cov,
		'cov_hcv': hcv_cov
	}

# Load all data once at startup
data_store = load_and_preprocess_data()

# Get the list of countries with missing WHO regions
missing_countries = get_countries_missing_who_region(data_store)

print(f"\nFound {len(missing_countries)} countries with missing WHO regions:")
for country in missing_countries:
    print(f"  - {country}")

# You can also export this list to a file if needed
with open("missing_who_regions.txt", "w") as f:
    for country in missing_countries:
        f.write(f"{country}\n")

# Genotype colors - defined once
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
	"Recombinant": "#e00603"
}

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

external_stylesheets = [
	"https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css",
	dbc.themes.BOOTSTRAP
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

app.layout = html.Div([
	dcc.Store(id="selected-virus", data="HBV"),
	dcc.Store(id='filtered-data-store'),
	dcc.Store(id='computed-metrics-store'),
	
	dbc.Container([
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
				dbc.Button("Download Data", id="btn-download", color="success", size="lg"),
				dcc.Download(id="download-data")
			], width="auto"),
		], className="mb-4 align-items-center", justify="between"),

		# Filters Section - Improved layout
		dbc.Row([
			dbc.Col([
				dbc.Card([
					dbc.CardBody([
						html.H6("FILTERS", className="text-muted mb-3"),
						
						# Year Slider
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
						
						# Dropdown Filters - Improved layout
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

		# Display Toggles - Better organization
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
									className="d-inline"
								)
							], width=6),
							
							dbc.Col([
								html.Label("Map Mode:", className="fw-bold me-2"),
								dcc.RadioItems(
									id="map-mode",
									options=[
										{"label": " Sequences", "value": "sequences"},
										{"label": " Coverage", "value": "coverage"}
									],
									value="sequences",
									inline=True,
									className="d-inline"
								)
							], width=6)
						])
					])
				], className="mb-3 shadow-sm")
			], width=12)
		]),

		# Main Content
		dbc.Row([
			# Indicators
			build_indicators(virus="HBV"),
			
			# Map Section
			dbc.Col([
				dbc.Card([
					dbc.CardBody([
						html.Div([
							html.H4(id="map-title-main", className="mb-1"),
							html.Small(id="map-title-sub", className="text-muted")
						], className="mb-3"),
						
						# IHME Metric Selector + Toggle
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
							], width="auto"),
							dbc.Col([
								dbc.Checklist(
									id="show-ihme-toggle",
									options=[{"label": "Show IHME on Sequences map", "value": "ihme"}],
									value=[],	# default OFF → sequences map shows sequences
									switch=True,
								)
							], width="auto", className="pt-4")
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
		
		# Middle Row - Burden and Gap Analysis
		dcc.Store(id="gap-store"),
		dbc.Row([
			dbc.Col([
				dbc.Card([
					dbc.CardBody([
						html.H5(id="burden-title-main", className="mb-3"),
						dcc.Loading(
							dcc.Graph(id="global-burden-lineplot"),
							type="circle"
						)
					])
				], className="h-100 shadow-sm")
			], width=6),
			
			dbc.Col([
				dbc.Card([
					dbc.CardBody([
						html.H5(id="gap-title-main", className="mb-2"),
					
						# Small control row for the bar
						dbc.Row([
							dbc.Col([
								html.Label("Order:", className="fw-bold me-2"),
								dcc.Dropdown(
									id="coverage-order",
									options=[
										{"label": "Lowest coverage (largest gaps)", "value": "lowest"},
										{"label": "Highest coverage (best covered)", "value": "highest"},
									],
									value="lowest", clearable=False, style={"width": "290px"}
								)
							], width="auto"),
					
							dbc.Col([
								html.Label("Top N:", className="fw-bold me-2"),
								dcc.Input(
									id="coverage-top-n",
									type="number", min=5, max=50, step=5, value=20,
									style={"width": "100px"}
								)
							], width="auto"),
					
							# Download button + Download component
							dbc.Col([
								dbc.Button(
									"Download table (CSV)",
									id="coverage-download-btn",
									color="secondary",
									className="ms-2"
								),
								dcc.Download(id="coverage-download")
							], width="auto"),
						], className="g-3 mb-2"),
					
						dcc.Loading(
							dcc.Graph(
								id="undersequenced-bar",
								config={"displayModeBar": True, "displaylogo": False},
								style={"height": "450px"}
							),
							type="circle"
						),
					])
				], className="h-100 shadow-sm")
			], width=6)
		], className="mb-4"),
		
		# Time Series and Pie Chart
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
						html.H5(id="pie-title-main", className="mb-3"),
						dcc.Loading(
							dcc.Graph(id="country-pie-chart", style={"height": "100%"}),
							type="circle"
						)
					])
				], className="h-100 shadow-sm")
			], width=6)
		], className="mb-4"),
		
		# Bottom Row - Genotype and Mutations
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
						html.H5(id="mutation-title-main", className="mb-3"),
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
		], className="mb-4")
		
	], fluid=True),
	
	# Footer Section - Optimized
	html.Footer([
		dbc.Container([
			dbc.Row([
				# Logo Section - Centered and responsive
				dbc.Col([
					html.Div([
						html.Img(src=app.get_asset_url('sun_logo.png'), 
								className="footer-logo mx-2",
								style={"height": "45px", "objectFit": "contain"}),
						html.Img(src=app.get_asset_url('ceri_logo.png'), 
								className="footer-logo mx-2",
								style={"height": "45px", "objectFit": "contain"}),
						html.Img(src=app.get_asset_url('CRICK.png'), 
								className="footer-logo mx-2",
								style={"height": "45px", "objectFit": "contain"}),
						html.Img(src=app.get_asset_url('AHRI_logo.png'), 
								className="footer-logo mx-2",
								style={"height": "45px", "objectFit": "contain"})
					], className="d-flex justify-content-center align-items-center flex-wrap")
				], width=12, className="mb-3")
			]),
			
			dbc.Row([
				# Text Content - Properly centered
				dbc.Col([
					html.Div([
						html.P([
							f"© {datetime.now().year} Hepatitis Virus Sequence Dashboard. ",
							html.Span("All rights reserved.", className="text-muted")
						], className="mb-1"),
						html.P("Developed by Derek Tshiabuila, Vagner Fonseca, Eduan Wilkinson, Tulio de Oliveira", 
							  className="mb-1 text-muted"),
						html.A("GitHub Repository", 
							  href="https://github.com/your-repo", 
							  target="_blank", 
							  className="text-decoration-none text-primary",
							  style={"fontSize": "0.9rem"})
					], className="text-center")
				], width=12)
			])
		], fluid=True)
	], className="bg-light py-4 mt-5 border-top", 
	   style={"marginTop": "2rem !important"})
	
	], style={'backgroundColor': '#f8f9fa', 'minHeight': '100vh', 'padding': '20px'})


@app.callback(
	Output("hbv-mutation-section", "style"),
	Output("hcv-mutation-section", "style"),
	Input("selected-virus", "data")
)
def toggle_mutation_sections(selected_virus):
	if selected_virus == "HBV":
		return {"display": "block"}, {"display": "none"}
	elif selected_virus == "HCV":
		return {"display": "none"}, {"display": "block"}
	else:
		return {"display": "none"}, {"display": "none"}

@app.callback(
	Output("mutation-filter-dropdown", "options"),
	Output("mutation-filter-dropdown", "placeholder"),
	Input("selected-virus", "data")
)
def update_mutation_filter_options(virus):
	if virus == "HBV":
		options = [{"label": drug, "value": drug} for drug in sorted(data_store["hbv_mut"]["drug"].dropna().unique())]
		return options, "Select Drug"
	elif virus == "HCV":
		options = [{"label": gene, "value": gene} for gene in sorted(data_store["hcv_mut"]["gene"].dropna().unique())]
		return options, "Select Gene"
	else:
		return [], "Select Drug or Gene"

@app.callback(
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

# --- Burden-adjusted coverage helpers ---

BURDEN_MEASURE_FALLBACK = "Prevalence|Number"  # used if dropdown missing

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
		return pd.DataFrame(columns=["Country_standard", "observed_sequences", "burden", "expected_sequences", "coverage_gap"])

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
	df["expected_sequences"] = (df["burden"] / 10000.0) * target_per_10k
	df["coverage_gap"] = df["expected_sequences"] - df["observed_sequences"]
	df.loc[df["coverage_gap"] < 0, "coverage_gap"] = 0	# clip negative
	
	who_map = (
		ihme_df[["Country_standard", "WHO_Regions"]]
		.dropna(subset=["Country_standard"])
		.drop_duplicates("Country_standard")
	)
	
	df = df.merge(who_map, on="Country_standard", how="left")
	return df

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

@app.callback(
	Output("gap-store", "data"),
	Input("selected-virus", "data"),
	Input("year-slider", "value"),
	Input("continent-dropdown", "value"),
	Input("country-dropdown", "value"),
	Input("genotype-dropdown", "value"),
	Input("ihme-metric-type", "value"),
	prevent_initial_call=False,
)
def compute_and_store_gap(virus, selected_years, continents, countries, genotypes, ihme_metric_choice):
	# ---- fallbacks ----
	virus = (virus or "HBV").upper()
	years = selected_years or [
		int(data_store["hbv_data"]["Year"].min() if virus == "HBV" else data_store["hcv_data"]["Year"].min()),
		int(data_store["hbv_data"]["Year"].max() if virus == "HBV" else data_store["hcv_data"]["Year"].max()),
	]
	ihme_metric_choice = ihme_metric_choice or BURDEN_MEASURE_FALLBACK

	# ---- filter sequence DF once ----
	base = data_store["hbv_data"] if virus == "HBV" else data_store["hcv_data"]
	df = base[(base["Year"] >= years[0]) & (base["Year"] <= years[1])].copy()
	if continents:
		df = df[df["WHO_Regions"].isin(continents)]
	if countries:
		df = df[df["Country_standard"].isin(countries)]
	if genotypes:
		df = df[df["Genotype"].isin(genotypes)]

	# ---- build gap_df ----
	gap_df = compute_gap_df(
		virus=virus,
		filtered_seq_df=df,
		ihme_df=data_store["ihme_df"],
		selected_years=years,
		who_regions=continents,
		countries=countries,
		ihme_metric_choice=ihme_metric_choice,
		target_per_10k=5.0,
	)

	# ---- normalize schema / numerics / WHO region ----
	expected_cols = [
		"Country_standard", "WHO_Regions",
		"observed_sequences", "burden", "expected_sequences", "coverage_gap"
	]

	if not isinstance(gap_df, pd.DataFrame) or gap_df.empty:
		gap_df = pd.DataFrame(columns=expected_cols)
	else:
		gap_df = gap_df.copy()

		# add WHO_Regions if missing
		if "WHO_Regions" not in gap_df.columns:
			who_map = (
				data_store["ihme_df"][["Country_standard", "WHO_Regions"]]
				.dropna(subset=["Country_standard"])
				.drop_duplicates("Country_standard")
			)
			gap_df = gap_df.merge(who_map, on="Country_standard", how="left")

		# coerce numerics
		for col in ["observed_sequences", "burden", "expected_sequences", "coverage_gap"]:
			if col in gap_df.columns:
				gap_df[col] = pd.to_numeric(gap_df[col], errors="coerce")

		# drop invalids and placeholders
		if "coverage_gap" in gap_df.columns:
			gap_df = gap_df.dropna(subset=["coverage_gap"])
		if "Country_standard" in gap_df.columns:
			gap_df = gap_df[(gap_df["Country_standard"].notna()) & (gap_df["Country_standard"] != "Unknown")]

	# (optional) keep only expected cols if you want a tidy payload
	# gap_df = gap_df[[c for c in expected_cols if c in gap_df.columns]]

	return _df_to_json(gap_df)



@app.callback(
    Output("undersequenced-bar", "figure"),
    Input("gap-store", "data"),
    Input("selected-virus", "data"),
    Input("coverage-order", "value"),
    Input("coverage-top-n", "value"),
)
def render_gap_bar(gap_json, virus, order, top_n):
    virus = (virus or "HBV").upper()
    order = order or "lowest"
    top_n = int(top_n or 20)

    df = _df_from_json(gap_json)
    if df.empty:
        return go.Figure().update_layout(
            height=450,
            annotations=[dict(text="No data available for current filters",
                              x=0.5, y=0.5, showarrow=False)]
        )

    # --- normalize expected columns & dtypes ---
    # allow for different casings or missing columns
    colmap = {c.lower(): c for c in df.columns}
    def col(name):  # case-insensitive getter
        return colmap.get(name.lower(), name)

    for c in ["observed_sequences", "expected_sequences", "coverage_gap"]:
        if col(c) in df.columns:
            df[col(c)] = pd.to_numeric(df[col(c)], errors="coerce")

    # compute Coverage_Ratio if missing (case-insensitively)
    if "coverage_ratio" in colmap:
        df["Coverage_Ratio"] = pd.to_numeric(df[col("coverage_ratio")], errors="coerce")
    else:
        # create it from observed/expected
        obs = df[col("observed_sequences")] if col("observed_sequences") in df.columns else np.nan
        exp = df[col("expected_sequences")] if col("expected_sequences") in df.columns else np.nan
        df["Coverage_Ratio"] = np.where(pd.to_numeric(exp, errors="coerce") > 0,
                                        pd.to_numeric(obs, errors="coerce") / pd.to_numeric(exp, errors="coerce"),
                                        np.nan)

    # basic clean
    if col("coverage_gap") in df.columns:
        df = df.dropna(subset=[col("coverage_gap")])
    df["Coverage_Ratio"] = pd.to_numeric(df["Coverage_Ratio"], errors="coerce")

    # --- select / order top N ---
    if order == "lowest":
        # “under-sequenced” = biggest additional genomes needed
        df = df.sort_values(col("coverage_gap"), ascending=False).head(top_n)
        df = df.sort_values(col("coverage_gap"), ascending=True)  # for horizontal bar (small→large)
        x_col = col("coverage_gap")
        x_label = "Estimated additional genomes needed"
        title = f"Top {len(df)} under-sequenced countries ({virus})"
    else:
        # “best covered” = highest coverage ratio
        df = df.sort_values("Coverage_Ratio", ascending=False).head(top_n)
        df = df.sort_values("Coverage_Ratio", ascending=True)      # small→large left→right
        x_col = "Coverage_Ratio"
        x_label = "Coverage ratio (Observed / Expected)"
        title = f"Top {len(df)} best-covered countries ({virus})"

    color = "#3182bd" if virus == "HBV" else "#d94801"
    fig = px.bar(
        df,
        x=x_col,
        y="Country_standard",
        orientation="h",
        color_discrete_sequence=[color],
        labels={x_col: x_label, "Country_standard": ""},
        title=title,
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



@app.callback(
	Output("coverage-download", "data"),
	Input("coverage-download-btn", "n_clicks"),
	State("gap-store", "data"),
	State("coverage-order", "value"),
	State("coverage-top-n", "value"),
	State("selected-virus", "data"),
	prevent_initial_call=True,
)

def download_coverage_table(n_clicks, gap_json, order, top_n, virus):
	if not n_clicks:
		return dash.no_update

	df = _df_from_json(gap_json)
	if df.empty:
		empty = pd.DataFrame(columns=[
			"Country_standard","WHO_Regions","observed_sequences","Estimated_Infections",
			"expected_sequences","coverage_gap","Coverage_Ratio","Coverage_per10k"
		])
		return dcc.send_data_frame(empty.to_csv, "coverage_table_empty.csv", index=False)

	order = (order or "lowest")
	top_n = int(top_n or 20)
	virus = (virus or "HBV").upper()

	# --- case-insensitive column lookup ---
	colmap = {c.lower(): c for c in df.columns}
	def col(name): return colmap.get(name.lower(), name)

	# numeric coercions
	for c in ["observed_sequences","burden","expected_sequences","coverage_gap"]:
		if col(c) in df.columns:
			df[col(c)] = pd.to_numeric(df[col(c)], errors="coerce")

	# remove placeholders & invalids
	if col("Country_standard") in df.columns:
		df = df[df[col("Country_standard")] != "Unknown"]
	if col("coverage_gap") in df.columns:
		df = df.dropna(subset=[col("coverage_gap")])

	# derived metrics (robust naming)
	df["Estimated_Infections"] = df.get(col("burden"))
	df["Coverage_Ratio"] = np.where(
		pd.to_numeric(df.get(col("expected_sequences")), errors="coerce") > 0,
		pd.to_numeric(df.get(col("observed_sequences")), errors="coerce") /
		pd.to_numeric(df.get(col("expected_sequences")), errors="coerce"),
		np.nan
	)
	df["Coverage_per10k"] = np.where(
		pd.to_numeric(df["Estimated_Infections"], errors="coerce") > 0,
		pd.to_numeric(df.get(col("observed_sequences")), errors="coerce") /
		pd.to_numeric(df["Estimated_Infections"], errors="coerce") * 10_000,
		np.nan
	)

	# order + take top N
	if order == "lowest":  # largest gaps
		df_sorted = df.sort_values(col("coverage_gap"), ascending=False).head(top_n)
		tag = "lowest"
	else:				   # best covered = highest ratio
		df_sorted = df.sort_values("Coverage_Ratio", ascending=False).head(top_n)
		tag = "highest"

	# rounding for CSV readability
	for c in ["observed_sequences","expected_sequences","coverage_gap","Estimated_Infections"]:
		if c in df_sorted.columns:
			df_sorted[c] = pd.to_numeric(df_sorted[c], errors="coerce").round(0).astype("Int64")
	if "Coverage_Ratio" in df_sorted.columns:
		df_sorted["Coverage_Ratio"] = pd.to_numeric(df_sorted["Coverage_Ratio"], errors="coerce").round(3)
	if "Coverage_per10k" in df_sorted.columns:
		df_sorted["Coverage_per10k"] = pd.to_numeric(df_sorted["Coverage_per10k"], errors="coerce").round(2)

	cols = [c for c in [
		"Country_standard",
		"WHO_Regions",
		"observed_sequences",
		"Estimated_Infections",
		"expected_sequences",
		"coverage_gap",
		"Coverage_Ratio",
		"Coverage_per10k",
	] if c in df_sorted.columns]

	filename = f"{virus}_coverage_table_{tag}_top{min(top_n, len(df_sorted))}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
	return dcc.send_data_frame(df_sorted[cols].to_csv, filename, index=False)




@app.callback(
	Output("gap-title-main", "children"),
	Input("selected-virus", "data"),
)
def set_gap_card_title(virus):
	virus = (virus or "HBV").upper()
	return f"Burden‑adjusted coverage — {virus}"

	
@app.callback(
	Output("genotype-bar-chart", "figure", allow_duplicate=True),
	[Input("genotype-bar-chart", "restyleData"),
	 Input("genotype-bar-chart", "relayoutData")],
	State("genotype-bar-chart", "figure"),
	prevent_initial_call=True
)
def keep_log_and_xaxis(_restyle, _relayout, fig):
	if not fig:
		return fig

	# --- persistent log y ---
	rng = fig.get("layout", {}).get("meta", {}).get("y_full_log_range")
	if rng:
		fig["layout"]["yaxis"].update(
			type="log", autorange=False, range=rng, tickformat="~s",
			gridcolor="rgba(0,0,0,0.08)", linecolor="rgba(0,0,0,0.25)", zeroline=False
		)

	# --- remove hidden categories from xaxis ---
	visible = []
	for tr in fig.get("data", []):
		if tr.get("visible", True) != "legendonly":
			nm = tr.get("name") or (tr.get("x")[0] if tr.get("x") else None)
			if nm: visible.append(nm)

	all_cats = fig.get("layout", {}).get("meta", {}).get("all_categories", [])
	if visible:
		new_order = [c for c in all_cats if c in set(visible)]
	else:
		new_order = all_cats
	fig["layout"]["xaxis"].update(categoryorder="array", categoryarray=new_order)

	return fig
	
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


@app.callback(
	Output("indicator-total", "children"),
	Output("indicator-countries", "children"),
	Output("indicator-genotypes", "children"),
	Output("indicator-years", "children"),
	Output("year-slider", "min"),
	Output("year-slider", "max"),
	Output("year-slider", "value"),
	Output("year-slider", "marks"),
	Output("line-chart", "figure"),
	Output('continent-dropdown', 'options'),
	Output("country-dropdown", "options"),
	Output("genotype-dropdown", "options"),
	Output("country-pie-chart", "figure"),
	Output("genotype-bar-chart", "figure"),
	Output("genotype-map", "figure"),
	Output("selected-years-display", "children"),
	Output("map-title-main", "children"),
	Output("map-title-sub", "children"),
	Output("line-title-main", "children"),
	Output("pie-title-main", "children"),
	Output("bar-title-main", "children"),
	Output("mutation-title-main", "children"),
	Output("global-burden-lineplot", "figure"),
	Output("burden-title-main", "children"),
	Input("selected-virus", "data"),
	Input("year-slider", "value"),
	Input("btn-reset-time", "n_clicks"),
	Input("continent-dropdown", "value"),
	Input("country-dropdown", "value"),
	Input("genotype-dropdown", "value"),
	Input("display-mode", "value"),
	Input("map-mode","value"),
	Input("ihme-metric-type", "value"),
	Input("show-ihme-toggle", "value"),
	)


def update_indicators(selected_virus, selected_years, btn_reset_clicks,
					  selected_continents, selected_countries, selected_genotypes,
					  display_mode, map_mode, ihme_metric_type, show_ihme_toggle):

	df = data_store['hbv_data'] if selected_virus == "HBV" else data_store['hcv_data']

	# --- resolve selected_years (incl. reset button) ---
	ctx = callback_context
	triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None

	if triggered_id == "btn-reset-time":
		min_year_selected = int(df["Year"].min())
		max_year_selected = int(df["Year"].max())
		selected_years = [min_year_selected, max_year_selected]
	else:
		if not selected_years or len(selected_years) != 2:
			min_year_selected = int(df["Year"].min())
			max_year_selected = int(df["Year"].max())
			selected_years = [min_year_selected, max_year_selected]
		else:
			min_year_selected, max_year_selected = selected_years

	# --- filtered_df (apply all UI filters once) ---
	filtered_df = df[
		(df["Year"] >= min_year_selected) & (df["Year"] <= max_year_selected)
	].copy()

	if selected_continents:
		filtered_df = filtered_df[filtered_df["WHO_Regions"].isin(selected_continents)]
	if selected_countries:
		filtered_df = filtered_df[filtered_df["Country_standard"].isin(selected_countries)]
	if selected_genotypes:
		filtered_df = filtered_df[filtered_df["Genotype"].isin(selected_genotypes)]

	# --- top indicators, slider marks, line chart, pie, genotype bar (unchanged from your working versions) ---
	total_genomes = len(filtered_df)
	unique_countries = filtered_df['Country_standard'].nunique()
	unique_genotypes = filtered_df['Genotype'].nunique()

	year_min_full_range = int(df["Year"].min())
	year_max_full_range = int(df["Year"].max())

	all_years = list(range(year_min_full_range, year_max_full_range + 1))
	marks = {year: {'label': str(year), 'style': {'color': '#77b0b1'}} if year % 5 == 0 else '' for year in all_years}
	marks[min_year_selected] = {'label': str(min_year_selected), 'style': {'color': '#f50', 'font-weight': 'bold'}}
	marks[max_year_selected] = {'label': str(max_year_selected), 'style': {'color': '#f50', 'font-weight': 'bold'}}

	line_data = (
		filtered_df.groupby(["Year", "Genotype"])
		.size().reset_index(name="Genome Sequences")
	)
	
	genotype_colors = HBV_GENOTYPE_COLORS if selected_virus == "HBV" else HCV_GENOTYPE_COLORS
	
	fig = px.line(
		line_data,
		x="Year",
		y="Genome Sequences",
		color="Genotype",
		color_discrete_map=genotype_colors,
		markers=True
	)
	
	if not line_data.empty:
		year_min = int(line_data["Year"].min())
		year_max = int(line_data["Year"].max())
		fig.update_layout(
			xaxis=dict(
				tickmode="array",
				tickvals=list(range(year_min, year_max+1, 4))
			)
		)
	
	# ---- compute a fixed log range from ALL data (not affected by legend toggles) ----
	y_pos = line_data["Genome Sequences"].dropna()
	if len(y_pos):
		y_min = float(np.nanmax([1.0, np.nanmin(y_pos)]))	 # at least 1
		y_max = float(np.nanmax(y_pos))
		lo = np.log10(max(0.8, y_min/1.5))					 # small padding
		hi = np.log10(y_max*1.5)							 # headroom
	else:
		lo, hi = 0, 2  # fallback
	
	full_log_range = [lo, hi]
	
	fig.update_layout(
		xaxis_title="Year",
		yaxis=dict(
			title="Number of Sequences (log scale)",
			type="log",
			autorange=False,				# don't re-autoscale on legend toggle
			range=full_log_range,			# <- persistent
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
		# Store the canonical range in layout.meta so a callback can re-apply it:
		meta={"y_full_log_range": full_log_range}
	)
	all_continents_options = [
	{'label': r, 'value': r}
	for r in sorted(df["WHO_Regions"].dropna().unique())
	if r != "Unknown"
	]
	all_countries_options = [
	{'label': c, 'value': c}
	for c in sorted(df["Country_standard"].dropna().unique())
	if c != "Unknown"
	]
	all_genotypes_options  = [
	{'label': c, 'value': c}
	for c in sorted(df["Genotype"].dropna().unique())
	if c != "Unknown"
	]
	
	# --- genotype bar ---
	cyg = (
		filtered_df.groupby(["Country_standard", "Year", "Genotype"], as_index=False)
				   .size().rename(columns={"size": "Count"})
	)
	
	# 2) bring in population by nearest year (re-use your helper)
	# NOTE: this expects data_store["population_df"] to have columns: Country_standard, Year, Population
	cyg = merge_population_nearest_two_pass(
		cyg,
		data_store["population_df"],
		tol_years_first=3,	 # your earlier values
		tol_years_wide=50
	)
	
	# 3) aggregate to genotype level
	agg = (
		cyg.groupby("Genotype", as_index=False)
		   .agg(Total=("Count", "sum"),
				Pop=("Population", "sum"))
	)
	
	# 4) compute per-million (guard against zero/NaN population)
	agg["PerMillion"] = np.where(
		(agg["Pop"].notna()) & (agg["Pop"] > 0),
		(agg["Total"] / agg["Pop"]) * 1_000_000.0,
		np.nan
	)
	
	# 5) choose series depending on display_mode
	y_col = "PerMillion" if display_mode == "PerMillion" else "Total"
	y_title = "Sequences per Million" if display_mode == "PerMillion" else "Number of Sequences"
	
	# For log axis: replace 0 with None so Plotly doesn’t error on log(0)
	agg[y_col] = agg[y_col].where(agg[y_col] > 0, None)
	
	# Keep genotypes ordered alphabetically (or change to Total-desc if preferred)
	agg = agg.sort_values("Genotype")
	
	# Compute log-safe min/max from this agg
	vals = agg[y_col].dropna().astype(float)
	if len(vals):
		ymin = float(vals.min())
		ymax = float(vals.max())
		# safe padding for log axis
		low_pad	 = max(0.1, ymin / 1.2)	  # ~20% below min, never < 0.1
		high_pad = ymax * 1.3			  # ~30% above max
		full_log_range = [np.log10(low_pad), np.log10(high_pad)]
	else:
		full_log_range = [0, 2]	 # fallback
		
	bar_fig = px.bar(
		agg,
		x="Genotype",
		y=y_col,
		color="Genotype",
		color_discrete_map=genotype_colors,
	)
	
	bar_fig.update_traces(
		hovertemplate="%{x}<br>%{y:,} sequences<extra></extra>"
	)
	
	# Use explicit category order (all genotypes)
	cat_order = agg["Genotype"].tolist()
	bar_fig.update_xaxes(categoryorder="array", categoryarray=cat_order)
	
	bar_fig.update_layout(
		xaxis_title="Genotype",
		yaxis=dict(
			title=y_title + " (log scale)",
			type="log",
			autorange=False,
			range=full_log_range,		  # persistent with padding below min
			tickformat="~s",
			gridcolor="rgba(0,0,0,0.08)",
			zeroline=False,
			linecolor="rgba(0,0,0,0.25)"
		),
		bargap=0.25,					  # optional aesthetics
		plot_bgcolor="white",
		paper_bgcolor="white",
		meta={
			"y_full_log_range": full_log_range,
			"all_categories": agg["Genotype"].tolist()
		}
	)
	
	# 6) layout: clear background + log axis
	# hide zeros on log axis
	bar_fig.update_traces(hovertemplate="%{x}<br>%{y:,} sequences<extra></extra>")
	bar_fig.for_each_trace(lambda t: t.update(y=None) if t.y is not None and np.any(np.array(t.y)<=0) else ())
	
	bar_fig.update_layout(
		xaxis_title="Genotype",
		yaxis=dict(
			title="Number of Sequences (log scale)",
			type="log",
			tickformat="~s",			  # 10, 100, 1k, 10k…
			gridcolor="rgba(0,0,0,0.08)"
		),
		plot_bgcolor="white",
		paper_bgcolor="white",
		height=400,
		margin=dict(t=100, b=0, l=0, r=0),
		legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center", yanchor="bottom")
	)
	
	# ==============================
	# NEW: compute burden-adjusted gap_df (once)
	# ==============================
	ihme_df = data_store["ihme_df"]
	gap_df = compute_gap_df(
		virus=selected_virus,
		filtered_seq_df=filtered_df,
		ihme_df=ihme_df,
		selected_years=[min_year_selected, max_year_selected],
		who_regions=selected_continents,
		countries=selected_countries,
		ihme_metric_choice=ihme_metric_type or BURDEN_MEASURE_FALLBACK,
		target_per_10k=5.0
	)
	gap_title_main = f"Lowest burden-adjusted coverage — {selected_virus.upper()}"
	#underseq_fig = undersequenced_bar(gap_df, virus=selected_virus, n=20)
	
	# --- PIE CHART (auto-switch by map-mode) ---
	def _nice_pct_labels(pcts):
		return [f"{p:.1f}%" for p in pcts] if len(pcts) else []
	
	if map_mode == "coverage":
		# Use IHME expected infections (Number). Prefer the dropdown *if* it's Number; else fallback.
		choice = ihme_metric_type or "Prevalence|Number"
		try:
			ihme_measure, ihme_metric = choice.split("|")
		except Exception:
			ihme_measure, ihme_metric = "Prevalence", "Number"
		if ihme_metric != "Number":
			ihme_measure, ihme_metric = "Prevalence", "Number"
	
		cause_lookup = {
			"HBV": "Total burden related to hepatitis B",
			"HCV": "Total burden related to hepatitis C",
		}
		cause_filter = cause_lookup.get(selected_virus.upper(), "")
	
		ihme_df = data_store["ihme_df"].copy()
		base = ihme_df[
			(ihme_df["sex"] == "Both") &
			(ihme_df["age"] == "All ages") &
			(ihme_df["cause"] == cause_filter) &
			(ihme_df["measure"] == ihme_measure) &
			(ihme_df["metric"] == "Number")
		].copy()
	
		# Align with UI filters (optional but keeps views consistent)
		if selected_continents:
			base = base[base["WHO_Regions"].isin(selected_continents)]
		if selected_countries:
			base = base[base["Country_standard"].isin(selected_countries)]
	
		# Use latest year in the selected range (or latest overall if empty)
		if selected_years and len(selected_years) == 2:
			y0, y1 = int(selected_years[0]), int(selected_years[1])
			base = base[(base["year"] >= y0) & (base["year"] <= y1)]
		if not base.empty:
			base = base[base["year"] == base["year"].max()]
	
		if base.empty:
			# Empty, make a blank pie
			pie_fig = px.pie(pd.DataFrame({"Country": [], "Value": []}), names="Country", values="Value", hole=0.5)
			pie_title_main = f"Top 10 countries by estimated infections (IHME) — {selected_virus.upper()}"
		else:
			per_country = (
				base.groupby("Country_standard", as_index=False)["val"]
					.sum().rename(columns={"Country_standard": "Country", "val": "Estimated"})
			)
			per_country = per_country.sort_values("Estimated", ascending=False)
	
			total_est = float(per_country["Estimated"].sum()) if not per_country.empty else 0.0
			top10 = per_country.head(10).copy()
			top10["SharePct"] = (top10["Estimated"] / total_est * 100.0) if total_est > 0 else 0.0
	
			pie_fig = px.pie(
				top10,
				names="Country",
				values="Estimated",
				hole=0.5,
				color="Country"
			)
			pie_fig.update_traces(
				text=_nice_pct_labels(top10["SharePct"].tolist()),
				texttemplate="%{text}",
				textinfo="text",
				sort=False,
				marker=dict(line=dict(color="#fff", width=1)),
				hovertemplate="<b>%{label}</b><br>"
							  "Share of expected infections: %{customdata[0]:.1f}%<br>"
							  "Estimated infections: %{value:,.0f}<extra></extra>",
				customdata=top10[["SharePct"]]
			)
	
			if len(per_country) > 10:
				excluded_n = len(per_country) - 10
				excluded_sum = int(per_country["Estimated"].iloc[10:].sum())
				pie_fig.add_annotation(
					text=f"<i>+ {excluded_n} more countries ({excluded_sum:,} infections not shown)</i>",
					x=0.5, y=-0.12, showarrow=False, font=dict(size=10, color="gray")
				)
	
			pie_title_main = (
				f"Top 10 countries by estimated infections (IHME burden-adjusted) — {selected_virus.upper()}"
				f"<br><span style='font-size:0.85em;color:gray'>Measure: {ihme_measure} (Number); "
				f"Year: {int(base['year'].max()) if not base.empty else '—'}</span>"
			)
	
	else:
		# Sequences (current mode) pie
		vc = filtered_df["Country_standard"].value_counts(dropna=False)
		country_counts = vc.rename_axis("Country").reset_index(name="Count")  # → columns: Country, Count
		country_counts["Count"] = pd.to_numeric(country_counts["Count"], errors="coerce").fillna(0).astype(int)
		
		if country_counts.empty:
			pie_fig = px.pie(pd.DataFrame({"Country": [], "Count": []}), names="Country", values="Count", hole=0.5)
			pie_title_main = f"Top 10 countries by {selected_virus.upper()} sequences (current selection)"
		else:
			total_count = int(country_counts["Count"].sum())
			top10 = country_counts.head(10).copy()
			top10["SharePct"] = (top10["Count"] / total_count * 100.0) if total_count > 0 else 0.0
		
			pie_fig = px.pie(top10, names="Country", values="Count", hole=0.5, color="Country")
			pie_fig.update_traces(
				text=[f"{p:.1f}%" for p in top10["SharePct"]],
				texttemplate="%{text}",
				textinfo="text",
				sort=False,
				marker=dict(line=dict(color="#fff", width=1)),
				hovertemplate="<b>%{label}</b><br>"
							  "Share of total sequences: %{customdata[0]:.1f}%<br>"
							  "Sequences: %{value:,}<extra></extra>",
				customdata=top10[["SharePct"]],
			)
			
			pie_fig.update_layout(
				legend=dict(
					orientation="v",	  # vertical legend
					yanchor="top",		  # anchor legend's top
					y=1,				  # place at top of chart
					xanchor="left",		  # anchor to the left side
					x=-0.3				  # move it outside the plot to the left
				)
			)
		
			if len(country_counts) > 10:
				excluded_n = len(country_counts) - 10
				excluded_sum = int(country_counts["Count"].iloc[10:].sum())
				pie_fig.add_annotation(
					text=f"<i>+ {excluded_n} more countries ({excluded_sum:,} sequences not shown)</i>",
					x=0.5, y=-0.12, showarrow=False, font=dict(size=10, color="gray")
				)
		
			pie_title_main = f"Top 10 countries by {selected_virus.upper()} sequences (current selection)"


	# ==============================
	# Map branch: coverage vs sequences
	# ==============================
	if map_mode == "coverage":
		# (unchanged) coverage ratio map
		coverage_for_map = gap_df.copy()
		coverage_for_map["coverage_ratio"] = (
			coverage_for_map["observed_sequences"] /
			coverage_for_map["expected_sequences"].replace({0: np.nan})
		)
		map_fig = create_coverage_map(coverage_for_map, data_store["coord_lookup"], virus_type=selected_virus)
		map_title_main = f"{selected_virus.upper()} Burden-adjusted Sequencing Coverage"
		map_title_sub  = "Coverage ratio = sequences / estimated infections"
	
	else:
		# NEW: use IHME metric in sequences mode too
		ihme_latest = ihme_latest_by_country(
			ihme_df=data_store["ihme_df"],
			virus=selected_virus,
			measure_metric=(ihme_metric_type or "Prevalence|Number"),
			regions=selected_continents,
			countries=selected_countries,
			years=selected_years
		)
	
		# USE THE CORRECT VARIABLE NAME FROM CALLBACK PARAMETER
		if not ihme_latest.empty and show_ihme_toggle:
			# DEBUG: Print column names to see what's available
			print("IHME DataFrame columns:", ihme_latest.columns.tolist())
			print("IHME DataFrame sample:")
			print(ihme_latest.head())
			
			# Show IHME burden map (prevalence/incidence/deaths) in sequences mode
			# Format data for create_world_map
			ihme_mapped = ihme_latest.copy()
			
			# Check which column contains the actual values - adjust based on what you see in the print output
			# Common column names might be: 'value', 'val', 'metric_value', 'prevalence', etc.
			value_column = None
			possible_columns = ['val', 'value', 'metric_value', 'prevalence', 'incidence', 'number', 'measure_value']
			
			for col in possible_columns:
				if col in ihme_mapped.columns:
					value_column = col
					break
			
			if value_column is None:
				# If none of the expected columns exist, use the first numeric column
				numeric_columns = ihme_mapped.select_dtypes(include=[np.number]).columns
				if len(numeric_columns) > 0:
					value_column = numeric_columns[0]
				else:
					# Fallback to sequence data if no numeric column found
					print("No numeric column found in IHME data, falling back to sequences")
					show_ihme_toggle = False
			
			if value_column:
				ihme_mapped["Metric_raw"] = ihme_mapped[value_column]  # Use the correct column name
				ihme_mapped["Metric"] = ihme_mapped[value_column]  # Keep original for display
				
				# For IHME data, we need to use a different display approach
				map_fig = create_world_map(
					ihme_mapped,
					pd.DataFrame(columns=["Country_standard", "Genotype", "Count"]),
					data_store["coord_lookup"],
					virus_type=selected_virus,
					display_mode="ihme"	 # SPECIAL MODE FOR IHME
				)
	
				# Title/subtitle using dropdown selection
				try:
					mm = (ihme_metric_type or "Prevalence|Number").split("|")
					ihme_measure, ihme_metric = mm[0], mm[1]
				except Exception:
					ihme_measure, ihme_metric = "Prevalence", "Number"
	
				latest_year_txt = int(ihme_latest["year"].max())
				map_title_main = f"{selected_virus.upper()} IHME Global Burden Map"
				map_title_sub  = f"{ihme_measure} ({ihme_metric}) — latest year: {latest_year_txt}"
	
		if not show_ihme_toggle or ihme_latest.empty:  # Moved this to handle the fallback case
			# Fallback to your original sequences map (raw / per million)
			if display_mode == "PerMillion":
				tmp = filtered_df[["Country_standard", "Year"]].copy()
				tmp["Country_standard"] = tmp["Country_standard"].astype(str)
				tmp["Year"] = pd.to_numeric(tmp["Year"], errors="coerce").astype("Int64")
				tmp = tmp.dropna(subset=["Year"]).astype({"Year": "int64"})
	
				country_year_counts = (
					tmp.groupby(["Country_standard", "Year"], as_index=False)
					   .size().rename(columns={"size": "Total"})
					   .sort_values(["Country_standard", "Year"])
				)
	
				country_year_counts = merge_population_nearest_two_pass(
					country_year_counts,
					data_store["population_df"],
					tol_years_first=3,
					tol_years_wide=50
				)
	
				country_year_counts["PerMillion"] = np.where(
					(country_year_counts["Population"].notna()) & (country_year_counts["Population"] > 0),
					(country_year_counts["Total"] / country_year_counts["Population"]) * 1_000_000,
					np.nan
				)
	
				latest_year_data = (
					country_year_counts.sort_values(["Country_standard", "Year"])
					.groupby("Country_standard", as_index=False)
					.tail(1)
					.reset_index(drop=True)
				)
	
				latest_year_data["Metric_raw"] = latest_year_data["PerMillion"]
				latest_year_data["Metric"] = latest_year_data["Metric_raw"].apply(
					lambda x: np.log10(x) if (np.isfinite(x) and x > 0) else np.nan
				)
			else:
				latest_year_data = (
					filtered_df.groupby("Country_standard", as_index=False)
							   .size().rename(columns={"size": "Total"})
				)
				latest_year_data["Metric_raw"] = latest_year_data["Total"].fillna(0)
				latest_year_data["Metric"] = latest_year_data["Metric_raw"]
	
			country_genotype_counts = (
				filtered_df.groupby(["Country_standard", "Genotype"], as_index=False)
						   .size().rename(columns={"size": "Count"})
			)
	
			map_fig = create_world_map(
				latest_year_data,
				country_genotype_counts,
				data_store["coord_lookup"],
				virus_type=selected_virus,
				display_mode=display_mode
			)
			map_title_main = f"{selected_virus.upper()} Whole Genome Sequence Map"
			map_title_sub  = "(per million)" if display_mode == "PerMillion" else "Sequences (count)"
	# === IHME burden (numbers for all three series) ===
	cause_lookup = {
		"HBV": "Total burden related to hepatitis B",
		"HCV": "Total burden related to hepatitis C",
	}
	cause_filter = cause_lookup.get(selected_virus.upper())
	base = ihme_df[
		(ihme_df["sex"] == "Both") &
		(ihme_df["age"] == "All ages") &
		(ihme_df["cause"] == cause_filter)
	].copy()
	if selected_continents:
		base = base[base["WHO_Regions"].isin(selected_continents)]
	if selected_countries:
		base = base[base["Country_standard"].isin(selected_countries)]

	def series(measure, metric="Number"):
		s = base[(base["measure"] == measure) & (base["metric"] == metric)].copy()
		if s.empty:
			return pd.DataFrame(columns=["year", f"{measure}_{metric}"])
		s["val"] = pd.to_numeric(s["val"], errors="coerce").fillna(0)
		return s.groupby("year", as_index=False)["val"].sum().rename(columns={"val": f"{measure}_{metric}"})

	d_prev = series("Prevalence", "Number")
	d_inc  = series("Incidence",  "Number")
	d_dea  = series("Deaths",	  "Number")

	from functools import reduce
	dfs = [d for d in [d_prev, d_inc, d_dea] if not d.empty]
	burden_df = reduce(lambda L, R: pd.merge(L, R, on="year", how="outer"), dfs).sort_values("year") if dfs else \
				pd.DataFrame(columns=["year","Prevalence_Number","Incidence_Number","Deaths_Number"])

	# Ensure numeric and clean zeros for log axis
	for col in ["Prevalence_Number", "Incidence_Number", "Deaths_Number"]:
		if col in burden_df.columns:
			burden_df[col] = pd.to_numeric(burden_df[col], errors="coerce").fillna(0)
	
	# Left-axis series (must be >0 for log)
	left_cols = [c for c in ["Prevalence_Number", "Incidence_Number"] if c in burden_df.columns]
	left_vals = []
	for c in left_cols:
		v = burden_df[c].values.astype(float)
		left_vals.extend(v[v > 0])	# only positive values count for log range
	
	# Set a persistent log range across both left series
	if len(left_vals):
		vmin = np.nanmin(left_vals)
		vmax = np.nanmax(left_vals)
		# pad 10x downwards if small, upwards by ~1.5× to create headroom
		lo = max(1.0, vmin / 10.0)
		hi = vmax * 1.5
		log_range = [np.log10(lo), np.log10(hi)]
	else:
		# sensible fallback if no data
		log_range = [0, 2]	# 10^0..10^2 (won’t really be used)
		
	# Build the figure
	burden_fig = go.Figure()
	
	# Prevalence (left, log)
	if "Prevalence_Number" in burden_df.columns:
		y_prev = burden_df["Prevalence_Number"].mask(burden_df["Prevalence_Number"] <= 0, None)
		burden_fig.add_trace(go.Scatter(
			x=burden_df["year"], y=y_prev,
			mode="lines+markers", name="Prevalence (number)",
			hovertemplate="Year: %{x}<br>Prevalence: %{y:,.0f}<extra></extra>"
		))
	
	# Incidence (left, log)
	if "Incidence_Number" in burden_df.columns:
		y_inc = burden_df["Incidence_Number"].mask(burden_df["Incidence_Number"] <= 0, None)
		burden_fig.add_trace(go.Scatter(
			x=burden_df["year"], y=y_inc,
			mode="lines+markers", name="Incidence (number)",
			hovertemplate="Year: %{x}<br>Incidence: %{y:,.0f}<extra></extra>"
		))
	
	# Deaths (right, linear)
	if "Deaths_Number" in burden_df.columns:
		burden_fig.add_trace(go.Scatter(
			x=burden_df["year"], y=burden_df["Deaths_Number"],
			mode="lines+markers", name="Deaths (number)",
			yaxis="y2",
			hovertemplate="Year: %{x}<br>Deaths: %{y:,.0f}<extra></extra>"
		))
	
	# Layout: log left axis (fixed range), linear right axis
	burden_fig.update_layout(
		height=400,
		margin=dict(t=50, b=20, l=60, r=60),
		legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
		xaxis=dict(
			title="Year",
			dtick=3, tickmode="linear",
			gridcolor="rgba(0,0,0,0.06)", showgrid=True
		),
		yaxis=dict(
			title="Prevalence / Incidence (number, log scale)",
			type="log",
			autorange=False,			# Disable auto-rescaling
			range=[7, 9],				# Fixed log range (10^7 to 10^9)
			tickmode="array",			# Force these ticks
			tickvals=[1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9],
			ticktext=["10M","20M","50M","100M","200M","500M","1G"],
			gridcolor="rgba(0,0,0,0.05)",
			#linecolor="rgba(0,0,0,0.25)",
			zeroline=False,
			fixedrange=True			   # LOCK THE SCALE (no zoom/double-click reset)
		),
		yaxis2=dict(
			title="Deaths (number)",
			overlaying="y", side="right",
			tickformat="~s",		  # e.g., 650k
			gridcolor="rgba(0,0,0,0.06)", showgrid=False
		),
		plot_bgcolor="white",
		paper_bgcolor="white"
	)
	
	burden_title_main = f"Global burdens in {selected_virus}: Prevalence, Incidence & Deaths"

	# titles
	selected_years_text = f"({min_year_selected} - {max_year_selected})"
	line_title_main = f"{selected_virus.upper()} Whole Genomes Per Year"
	pie_title_main	= f"Top 10 Countries by {selected_virus.upper()} Whole Genome Sequence Count"
	bar_title_main	= f"Total {selected_virus.upper()} Whole Genome Sequences by Genotype"
	mutation_title_main = f"{selected_virus} Drug Resistance Mutation Counts"

	# --- return in EXACT Output order ---
	return (
		f"{total_genomes:,}",			 # indicator-total.children
		unique_countries,				 # indicator-countries.children
		unique_genotypes,				 # indicator-genotypes.children
		f"{min_year_selected} - {max_year_selected}",  # indicator-years.children
		year_min_full_range,			 # year-slider.min
		year_max_full_range,			 # year-slider.max
		selected_years,					 # year-slider.value
		marks,							 # year-slider.marks
		fig,							 # line-chart.figure
		all_continents_options,			 # continent-dropdown.options
		all_countries_options,			 # country-dropdown.options
		all_genotypes_options,			 # genotype-dropdown.options
		pie_fig,						 # country-pie-chart.figure
		bar_fig,						 # genotype-bar-chart.figure
		map_fig,						 # genotype-map.figure
		f"({min_year_selected} - {max_year_selected})",	 # selected-years-display.children
		map_title_main,					 # map-title-main.children
		map_title_sub,					 # map-title-sub.children
		f"{selected_virus.upper()} Whole Genomes Per Year",	 # line-title-main.children
		pie_title_main,					 # pie-title-main.children
		f"Total {selected_virus.upper()} Whole Genome Sequences by Genotype",  # bar-title-main.children
		f"{selected_virus} Drug Resistance Mutation Counts",  # mutation-title-main.children
		burden_fig,						 # global-burden-lineplot.figure
		f"Global burdens in {selected_virus}: Prevalence, Incidence & Deaths"  # burden-title-main.children
	)

def _empty_world(message: str) -> go.Figure:
	fig = go.Figure()
	fig.update_layout(
		height=700,
		margin=dict(t=30, b=30, l=10, r=10),
		annotations=[dict(text=message, xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)],
	)
	fig.update_geos(
		projection_type="natural earth",
		showcountries=True,
		countrycolor="rgba(0,0,0,0.2)",
		showsubunits=True,
		domain=dict(x=[0, 1], y=[0.2, 1]),
	)
	return fig


def create_world_map(
	country_data: pd.DataFrame, 
	country_genotype_counts: pd.DataFrame, 
	coord_lookup: dict[str, dict[str, float]], 
	virus_type: str = "HBV", 
	display_mode: str = "raw"
) -> go.Figure:
	
	df = country_data.copy()
	if "Metric_raw" not in df.columns:
		df["Metric_raw"] = np.nan
	valid = df.dropna(subset=["Metric_raw"]).copy()

	fig = go.Figure()
	if valid.empty:
		return _empty_world("No country data available for current filters")

	if display_mode == "PerMillion":
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
	
	elif display_mode == "ihme":
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
			coords: dict[str, Any] | None = coord_lookup.get(country)
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


def create_coverage_map(cov_df: pd.DataFrame, coord_lookup: dict[str, any], virus_type: str = "HBV") -> go.Figure:
	valid = cov_df.copy()
	valid = valid[valid["Country_standard"].notna()]
	valid = valid[pd.to_numeric(valid["coverage_ratio"], errors="coerce").notna()]

	valid["coverage_capped"] = valid["coverage_ratio"].clip(upper=2.0)

	bins = [0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.01]
	labels = ["<0.10", "0.10–0.25", "0.25–0.50", "0.50–0.75", "0.75–1.0", "1.0–1.5", "≥1.5"]
	valid["bin"] = pd.cut(valid["coverage_capped"], bins=bins, labels=labels, include_lowest=True, right=False)
	bin_to_idx = {lab: i for i, lab in enumerate(labels)}
	valid["z"] = valid["bin"].map(bin_to_idx).astype(float)

	colors = ["#3b0f70", "#8c2981", "#de4968", "#f66e5b", "#fca636", "#f7d13d", "#e4f17a"]
	colorscale = [[i / (len(colors) - 1), c] for i, c in enumerate(colors)]

	fig = go.Figure()

	if not valid.empty:
		fig.add_trace(
			go.Choropleth(
				locations=valid["Country_standard"],
				locationmode="country names",
				z=valid["z"].astype(float).to_numpy(),
				zmin=0,
				zmax=len(labels) - 1,
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
						f"Observed: {r.get('observed_sequences', np.nan):,.0f}<br>"
						f"Expected: {r.get('expected_sequences', np.nan):,.0f}<br>"
						f"Coverage: {r['coverage_ratio']:.2f}×"
					),
					axis=1,
				),
				hoverinfo="text",
				marker_line_color="rgba(0,0,0,0.25)",
				marker_line_width=0.5,
				name="",
			)
		)
	else:
		fig.add_annotation(text="No coverage data for current filters", x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False)

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
	
	
@app.callback(
	Output("coverage-undersequenced-bar", "figure", allow_duplicate=True),
	[Input("coverage-undersequenced-bar", "restyleData"),
	 Input("coverage-undersequenced-bar", "relayoutData")],
	State("coverage-undersequenced-bar", "figure"),
	prevent_initial_call=True
)
def persist_undersequenced_yaxis(_restyle, _relayout, fig):
	if not fig:
		return fig

	rng = fig.get("layout", {}).get("meta", {}).get("y_full_log_range")
	if rng:
		fig["layout"]["yaxis"].update(
			type="log",
			autorange=False,
			range=rng,
			tickformat="~s",
			gridcolor="rgba(0,0,0,0.08)",
			zeroline=False,
			linecolor="rgba(0,0,0,0.25)"
		)
	return fig
	
@app.callback(
	Output("line-chart", "figure", allow_duplicate=True),
	[Input("line-chart", "restyleData"),
	 Input("line-chart", "relayoutData")],
	State("line-chart", "figure"),
	prevent_initial_call=True
)
def persist_line_log_axis(_restyle, _relayout, fig):
	if not fig:
		return fig

	# If double‑click set autorange=True or changed type, force back to log + fixed range
	saved = fig.get("layout", {}).get("meta", {}).get("y_full_log_range")
	if saved:
		fig["layout"].setdefault("yaxis", {})
		fig["layout"]["yaxis"].update(
			type="log",
			autorange=False,
			range=saved,
			tickformat="~s",
			gridcolor="rgba(0,0,0,0.08)",
			linecolor="rgba(0,0,0,0.25)",
			zeroline=False
		)
	return fig

@app.callback(
	Output("mutation-barplot", "figure"),
	[
		Input("year-slider", "value"),
		Input("continent-dropdown", "value"),
		Input("country-dropdown", "value"),
		Input("genotype-dropdown", "value"),
		Input("mutation-filter-dropdown", "value"),
		Input("selected-virus", "data")
	]
)
def update_mutation_barplot(selected_years, selected_continents, selected_countries,
						  selected_genotypes, selected_filter, virus, color=None):
						  
	if color is None:
		color = "#3182bd" if virus == "HBV" else "#d94801" if virus == "HCV" else "#3182bd"

	fig = go.Figure()
	
	if virus == "HBV":
		df = data_store["hbv_mut"].copy()
		col_name = "drug"
	else:
		df = data_store["hcv_mut"].copy()
		col_name = "gene"
	
	# Apply filters
	if selected_years:
		df = df[(df["Year"] >= selected_years[0]) & (df["Year"] <= selected_years[1])]
	if selected_continents:
		df = df[df["WHO_Regions"].isin(selected_continents)]
	if selected_countries:
		df = df[df["Country_standard"].isin(selected_countries)]
	if selected_genotypes:
		df = df[df["Genotype"].isin(selected_genotypes)]
	if selected_filter:
		df = df[df[col_name].notna()]
		selected_filter_lower = [str(s).strip().lower() for s in selected_filter]
		df = df[df[col_name].astype(str).str.strip().str.lower().isin(selected_filter_lower)]
	
	print(f"Final DataFrame shape after filtering: {df.shape}")
	
	if df.empty:
		fig.update_layout(
			title="No data available for selected filters",
			xaxis={"visible": False},
			yaxis={"visible": False},
			annotations=[{
				"text": "No mutations found for selected filters",
				"xref": "paper",
				"yref": "paper",
				"showarrow": False,
				"font": {"size": 16}
			}],
		)
		return fig
	
	mutation_counts = (
		df["mutation"]
		.value_counts()
		.reset_index()
	)
	mutation_counts.columns = ["Mutation", "Count"]
	mutation_counts = mutation_counts.sort_values("Count", ascending=False).head(20)  # Limit to top 20
	
	fig = px.bar(
		mutation_counts,
		x="Mutation",
		y="Count",
		labels={"Count": "Total Sequences", "Mutation": "Mutation"},
		color_discrete_sequence=[color]
	)
	fig.update_layout(
		xaxis_tickangle=-45,
		height=450,
		margin=dict(t=100, b=0, l=0, r=0),
		plot_bgcolor="white",
		paper_bgcolor="white",
		xaxis=dict(
			showgrid=False,
			linecolor="rgba(0,0,0,0.25)",	# subtle axis line
			#ticks="outside"
		),
		yaxis=dict(
			title="Total Sequences",
			gridcolor="rgba(0,0,0,0.08)",	# soft gray gridlines
			zeroline=False
		)
	)
	
	return fig

Output("global-burden-lineplot", "figure"),





def download_virus_data(n_clicks, selected_virus):
	df = data_store['hbv_data'] if selected_virus == "HBV" else data_store['hcv_data']
	return dcc.send_data_frame(df.to_csv, filename=f"{selected_virus}_data.csv", index=False)

if __name__ == '__main__':
	app.run(debug=True)