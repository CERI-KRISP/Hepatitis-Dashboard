import requests
import pandas as pd

def fetch_latest_population_data():
    url = "http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv"
    print("Downloading data from World Bank...")
    
    response = requests.get(url)
    with open("population_data.zip", "wb") as f:
        f.write(response.content)

    # Extract ZIP and load relevant file
    import zipfile
    with zipfile.ZipFile("population_data.zip", "r") as zip_ref:
        zip_ref.extractall("population_data")

    # Find the correct data file (usually starts with 'API_SP.POP.TOTL')
    import glob
    data_file = glob.glob("population_data/API_SP.POP.TOTL*.csv")[0]

    df = pd.read_csv(data_file, skiprows=4)

    # Melt and filter for recent year
    df_melted = df.melt(id_vars=["Country Name", "Country Code"], var_name="Year", value_name="Population")
    df_melted["Year"] = pd.to_numeric(df_melted["Year"], errors="coerce")
    latest_year = df_melted["Year"].max()

    df_latest = df_melted[df_melted["Year"] == latest_year].dropna()
    df_latest.rename(columns={"Country Name": "Country_standard"}, inplace=True)

    return df_latest[["Country_standard", "Year", "Population"]]

# Example usage
population_df = fetch_latest_population_data()
population_df.to_csv("population_by_country_year.csv", index=False)
