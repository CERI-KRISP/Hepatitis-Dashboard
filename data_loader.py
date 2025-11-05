import pandas as pd
import numpy as np
import re
from functools import lru_cache
import pycountry
import os

# Configuration
ENCODINGS = ['utf-8-sig', 'utf-8', 'latin1', 'iso-8859-1', 'cp1252', 'windows-1252', 'mac_roman']

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
	'HK': 'China',	   # sometimes appears as a short code in metadata
	'HKG': 'China',	   # ISO3 code occasionally found in CSVs
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

# === DATA LOADING ============================================================
# Pre-compile regex patterns
ACCESSION_PATTERN = re.compile(r"([A-Z0-9]+\.\d+)")
WHITESPACE_PATTERN = re.compile(r"\s+")

# Get the data directory path
def get_data_path(file_name):
    """Get absolute path to data file"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, 'data')
    return os.path.join(data_dir, file_name)

def load_with_encoding(file_name, sep='\t'):
    """Load file with encoding fallback and better error handling"""
    file_path = get_data_path(file_name)
    
    # Check if file exists first
    if not os.path.exists(file_path):
        print(f"⚠️  File not found: {file_path}")
        return pd.DataFrame()  # Return empty DataFrame instead of crashing
    
    for encoding in ENCODINGS:
        try:
            print(f"Trying encoding {encoding} for {file_name}...")
            df = pd.read_csv(file_path, sep=sep, encoding=encoding, 
                           on_bad_lines='skip', low_memory=False)
            print(f"✅ Successfully loaded {file_name} with {encoding} encoding")
            return df
        except Exception as e:
            print(f"Failed with encoding {encoding} for {file_name}: {e}")
            continue
    
    print(f"❌ Could not read {file_path} with any encoding")
    return pd.DataFrame()  # Return empty DataFrame instead of crashing

def load_recombinant_ids(file_name):
    """Load recombinant IDs efficiently"""
    file_path = get_data_path(file_name)
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"⚠️  Recombinant IDs file not found: {file_path}")
        return set()  # Return empty set
    
    try:
        with open(file_path) as f:
            return set(line.strip() for line in f)
    except Exception as e:
        print(f"❌ Error loading recombinant IDs: {e}")
        return set()

def load_csv_file(file_name, **kwargs):
    """Load CSV file with error handling"""
    file_path = get_data_path(file_name)
    
    if not os.path.exists(file_path):
        print(f"⚠️  File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path, **kwargs)
    except Exception as e:
        print(f"❌ Error loading {file_name}: {e}")
        return pd.DataFrame()

def process_genotype_column(df, recombs_set):
	"""Efficiently process genotype column"""
	mask = df['Taxa'].isin(recombs_set)
	df.loc[mask, 'Genotype'] = 'Recombinant'
	return df

def standardize_country_column(df, column_name):
	"""Standardize country column in a dataframe"""
	df['Country_standard'] = df[column_name].apply(standardize_country_cached)
	return df
	
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

def load_and_preprocess_data():
    print("📥 Loading data...")
    
    try:
        # Load main data files
        hbv_data = load_with_encoding("HBV_summary.tsv").dropna()
        hcv_data = load_with_encoding("HCV_summary.tsv").dropna()
        coords = load_with_encoding("WHO_regions_countries_coordinates.txt").dropna()

        # Check if files loaded successfully
        if hbv_data.empty:
            print("❌ HBV data is empty - check if file exists")
        if hcv_data.empty:
            print("❌ HCV data is empty - check if file exists")
        if coords.empty:
            print("❌ Coordinates data is empty - check if file exists")

        # Load recombinant IDs
        hbv_recombs = load_recombinant_ids("hbv_recombs.txt")
        hcv_recombs = load_recombinant_ids("hcv_recombs.txt")

        # Process genotype columns
        if not hbv_data.empty:
            hbv_data = process_genotype_column(hbv_data, hbv_recombs)
        if not hcv_data.empty:
            hcv_data = process_genotype_column(hcv_data, hcv_recombs)

        # Load additional data with better error handling
        population_df = load_csv_file("population_by_country_year.csv")
        ihme_df = load_csv_file("IHME-GBD_2021_DATA-9e7ec2c0-1.csv")
        hbv_cov = load_csv_file("hbv_burden_adjusted_coverage.csv")
        hcv_cov = load_csv_file("hcv_burden_adjusted_coverage.csv")

        # Load mutation data (handle missing files gracefully)
        hbv_mut = load_csv_file("hbv_resistant_mutations.csv")
        hcv_mut = load_csv_file("hcv_resistance_mutations.csv")

        print("🧹 Standardizing country names...")
        
        # Only process non-empty dataframes
        if not hbv_data.empty:
            hbv_data = standardize_country_column(hbv_data, 'Country')
        if not hcv_data.empty:
            hcv_data = standardize_country_column(hcv_data, 'Country')
        if not coords.empty:
            coords = standardize_country_column(coords, 'name')
        if not population_df.empty:
            population_df = standardize_country_column(population_df, 'Country_standard')
        if not ihme_df.empty:
            ihme_df = standardize_country_column(ihme_df, 'location')
        if not ihme_df.empty:
            for c in ["measure", "sex", "age", "cause", "metric"]:
                if c in ihme_df.columns:
                    ihme_df[c] = ihme_df[c].astype(str).str.strip()
            if "year" in ihme_df.columns:
                ihme_df["year"] = pd.to_numeric(ihme_df["year"], errors="coerce")
            for c in ["val", "upper", "lower"]:
                if c in ihme_df.columns:
                    ihme_df[c] = pd.to_numeric(ihme_df[c], errors="coerce")
        if not hbv_cov.empty:
            hbv_cov = standardize_country_column(hbv_cov, 'Country_standard')
        if not hcv_cov.empty:
            hcv_cov = standardize_country_column(hcv_cov, 'Country_standard')
    
        # Sanitize strings efficiently for non-empty dataframes
        for df in [hbv_data, hcv_data, coords]:
            if not df.empty:
                df["Country_standard"] = df["Country_standard"].astype(str).str.strip()
    
        # Drop any coords with missing standardization
        if not coords.empty:
            coords = coords.dropna(subset=["Country_standard"])
    
        print("📌 Merging WHO Regions...")
        
        # Only merge if we have data
        if not hbv_data.empty and not coords.empty:
            hbv_data = merge_who_region(hbv_data, coords)
        if not hcv_data.empty and not coords.empty:
            hcv_data = merge_who_region(hcv_data, coords)
        if not ihme_df.empty and not coords.empty:
            ihme_df = merge_who_region(ihme_df, coords)
    
        print("📆 Processing dates...")
        for df in [hbv_data, hcv_data]:
            if not df.empty:
                df["Year"] = pd.to_numeric(df["Date"], errors='coerce').fillna(0).astype('int32')
                df.dropna(subset=["Year"], inplace=True)
                if 'Date' in df.columns:
                    df.drop(columns=["Date"], inplace=True)
    
        # Fix population data if available
        if not population_df.empty:
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
        # Only group if we have data
        hbv_grouped = pd.DataFrame()
        hcv_grouped = pd.DataFrame()
        
        if not hbv_data.empty:
            grouping_cols = ["Year", "Genotype", "WHO_Regions", "Country_standard"]
            hbv_grouped = (hbv_data.groupby(grouping_cols, observed=True)
                           .size()
                           .reset_index(name="Count"))
            
            # Merge with population if available
            if not population_df.empty:
                hbv_grouped = pd.merge(hbv_grouped, population_df, on=["Country_standard", "Year"], how="left")
                hbv_grouped["PerMillion"] = (hbv_grouped["Count"] / hbv_grouped["Population"]) * 1_000_000
        
        if not hcv_data.empty:
            grouping_cols = ["Year", "Genotype", "WHO_Regions", "Country_standard"]
            hcv_grouped = (hcv_data.groupby(grouping_cols, observed=True)
                           .size()
                           .reset_index(name="Count"))
            
            if not population_df.empty:
                hcv_grouped = pd.merge(hcv_grouped, population_df, on=["Country_standard", "Year"], how="left")
                hcv_grouped["PerMillion"] = (hcv_grouped["Count"] / hcv_grouped["Population"]) * 1_000_000
    
        print("📍 Creating coordinate lookup...")
        coord_lookup = {}
        if not coords.empty:
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
    
        print("🧬 Processing mutation data...")
        # Only process mutation data if files were found
        hbv_mut_stats = pd.DataFrame()
        hcv_mut_stats = pd.DataFrame()
        
        if not hbv_mut.empty and not hbv_data.empty:
            # Extract accession IDs
            hbv_data["ID"] = extract_accession_id(hbv_data["Taxa"])
            
            # Handle mutation file column names
            if "ID" not in hbv_mut.columns and "Taxa" in hbv_mut.columns:
                hbv_mut["ID"] = extract_accession_id(hbv_mut["Taxa"])
            
            # Merge mutation data with main data
            merge_cols = ["ID", "Country_standard", "Year", "Genotype", "WHO_Regions"]
            hbv_mut = hbv_mut.merge(hbv_data[merge_cols], on="ID", how="left")
            hbv_mut = hbv_mut.dropna(subset=["Country_standard", "Year"])
            
            # Create mutation stats
            hbv_mut_stats = hbv_mut.groupby(["Country_standard", "drug", "mutation"], observed=True).size().reset_index(name="Count")
        
        if not hcv_mut.empty and not hcv_data.empty:
            hcv_data["ID"] = extract_accession_id(hcv_data["Taxa"])
            
            if "ID" not in hcv_mut.columns and "Taxa" in hcv_mut.columns:
                hcv_mut["ID"] = extract_accession_id(hcv_mut["Taxa"])
            
            merge_cols = ["ID", "Country_standard", "Year", "Genotype", "WHO_Regions"]
            hcv_mut = hcv_mut.merge(hcv_data[merge_cols], on="ID", how="left")
            hcv_mut = hcv_mut.dropna(subset=["Country_standard", "Year"])
            
            hcv_mut_stats = hcv_mut.groupby(["Country_standard", "gene", "mutation"], observed=True).size().reset_index(name="Count")
    
        print("✅ Done loading data.")
        
        # Return all data, even if some dataframes are empty
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
            'cov_hcv': hcv_cov,
            'hbv_summary_raw': hbv_data.copy(),
            'hcv_summary_raw': hcv_data.copy()
        }
        
    except Exception as e:
        print(f"❌ Critical error loading data: {e}")
        # Return empty data structure instead of crashing
        return {
            'hbv_data': pd.DataFrame(),
            'hcv_data': pd.DataFrame(),
            'coords': pd.DataFrame(),
            'hbv_grouped': pd.DataFrame(),
            'hcv_grouped': pd.DataFrame(),
            'coord_lookup': {},
            'population_df': pd.DataFrame(),
            'ihme_df': pd.DataFrame(),
            'hbv_mut': pd.DataFrame(),
            'hcv_mut': pd.DataFrame(),
            'hbv_mut_stats': pd.DataFrame(),
            'hcv_mut_stats': pd.DataFrame(),
            'cov_hbv': pd.DataFrame(),
            'cov_hcv': pd.DataFrame(),
            'hbv_summary_raw': pd.DataFrame(),
            'hcv_summary_raw': pd.DataFrame()
        }