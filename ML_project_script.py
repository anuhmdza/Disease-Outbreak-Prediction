# Auto-generated script from ML_project_enhanced.ipynb

from census import Census
from us import states
import pandas as pd

# Replace with your Census API key
CENSUS_API_KEY = "fe6843e76d70b7f6bcc3ad4885d926a06f69c5c4"

c = Census(CENSUS_API_KEY)

# Collect ACS data for all years
acs_all_years = []

for year in range(2016, 2023):
    data = c.acs5.get(
        ("NAME", "B19013_001E", "B25064_001E", "B27010_001E"),
        geo={"for": "state:*"},
        year=year
    )
    for row in data:
        row['year'] = year
    acs_all_years.extend(data)

# Convert to DataFrame
acs_df = pd.DataFrame(acs_all_years)

# Rename columns
acs_df.rename(columns={
    "NAME": "state_full",
    "B19013_001E": "median_household_income",
    "B25064_001E": "median_gross_rent",
    "B27010_001E": "health_insurance_coverage"
}, inplace=True)

# Add state abbreviations
fips_to_abbr = {str(s.fips).zfill(2): s.abbr for s in states.STATES}
acs_df['state'] = acs_df['state'].map(fips_to_abbr)

acs_df = acs_df[['state', 'year', 'median_household_income', 'median_gross_rent', 'health_insurance_coverage']]
# Save to CSV
acs_df.to_csv("census_data_2016_2022.csv", index=False)
print(" Census data saved as census_data_2016_2022.csv")





import pandas as pd
import requests
import os

# ----------------------------------------
# STEP 1: Load and Filter Inventory
# ----------------------------------------

inv_cols = ["station_id", "latitude", "longitude", "element", "first_year", "last_year"]
inventory = pd.read_csv("ghcnd-inventory.txt", delim_whitespace=True, names=inv_cols)
inv_filtered = inventory[inventory["element"].isin(["TMAX", "TMIN", "PRCP"])]

# Only stations with all 3 core elements
counts = inv_filtered.groupby("station_id")["element"].nunique().reset_index()
counts = counts[counts["element"] == 3]
station_elements = inv_filtered.merge(counts[["station_id"]], on="station_id", how="inner")
recent_stations = station_elements[station_elements["last_year"] >= 2022]

# ----------------------------------------
# STEP 2: Station Metadata and Best Stations
# ----------------------------------------

colspecs = [(0, 11), (12, 20), (21, 30), (31, 37), (38, 40), (41, 71)]
colnames = ["station_id", "lat", "lon", "elev", "state", "name"]
stations = pd.read_fwf("ghcnd-stations.txt", colspecs=colspecs, names=colnames)

merged = recent_stations.merge(stations, on="station_id")
us_states = [
    'AL','AK','AZ','AR','CA','CO','CT','DE','FL','GA','HI','ID','IL','IN','IA','KS','KY','LA','ME','MD',
    'MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC',
    'SD','TN','TX','UT','VT','VA','WA','WV','WI','WY','DC'
]
merged = merged[merged["state"].isin(us_states)]
best_stations = merged.sort_values(by=["state", "elev"]).groupby("state").first().reset_index()
best_stations.to_csv("best_stations_per_state.csv", index=False)

# ----------------------------------------
# STEP 3: Download .dly Files
# ----------------------------------------

os.makedirs("dly_files", exist_ok=True)
base_url = "https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all"

for station_id in best_stations["station_id"]:
    file_path = f"dly_files/{station_id}.dly"
    if not os.path.exists(file_path):
        url = f"{base_url}/{station_id}.dly"
        r = requests.get(url)
        if r.status_code == 200:
            with open(file_path, "wb") as f:
                f.write(r.content)
            print(f" Downloaded {station_id}.dly")
        else:
            print(f"Failed to download {station_id}")

# ----------------------------------------
# STEP 4: Parse .dly Files for 2016–2022
# ----------------------------------------

def parse_dly_file(file_path, station_id, state, years=range(2016, 2023)):
    records = []
    with open(file_path, "r") as file:
        for line in file:
            year = int(line[11:15])
            month = int(line[15:17])
            element = line[17:21]

            if year not in years or element not in ["TMAX", "TMIN", "PRCP"]:
                continue

            for day in range(1, 32):
                value_str = line[21 + (day - 1) * 8 : 26 + (day - 1) * 8]
                try:
                    value = int(value_str[:5])
                    if value == -9999:
                        continue
                    date = pd.to_datetime(f"{year}-{month:02d}-{day:02d}", errors="coerce")
                    if pd.notna(date):
                        records.append({
                            "date": date,
                            "state": state,
                            "element": element,
                            "value": value / 10
                        })
                except:
                    continue
    return pd.DataFrame(records)

# Parse all files
all_dfs = []
for _, row in best_stations.iterrows():
    path = f"dly_files/{row['station_id']}.dly"
    if os.path.exists(path):
        df = parse_dly_file(path, row["station_id"], row["state"])
        all_dfs.append(df)

df_all = pd.concat(all_dfs)

# ----------------------------------------
# STEP 5: Pivot + Monthly Aggregation
# ----------------------------------------

df_all["year"] = df_all["date"].dt.year
df_all["month"] = df_all["date"].dt.month

monthly = df_all.groupby(["state", "year", "month", "element"])["value"].mean().unstack().reset_index()
monthly.columns.name = None
monthly.rename(columns={"PRCP": "Precip_mm", "TMAX": "Max_Temp_C", "TMIN": "Min_Temp_C"}, inplace=True)

monthly.to_csv("climate_monthly_2016_2022.csv", index=False)
print(" All done! Saved: climate_monthly_2016_2022.csv")


import pandas as pd

# Try reading it as a tab-separated file instead of Excel
nndss_raw = pd.read_csv("NNDSS Annual Summary Data 2016-2022 (2).xls", sep="\t", engine="python", encoding="utf-8")

# Display column names
print(nndss_raw.columns.tolist())

# Quick look
nndss_raw


import pandas as pd

# Load the file (tab-delimited)
nndss_df = pd.read_csv("NNDSS Annual Summary Data 2016-2022 (2).xls", sep="\t", engine="python")

# Keep only the columns we need
nndss_df = nndss_df[["Year", "Regions/States", "Disease", "Case Count"]]

# Rename for clarity
nndss_df.columns = ["year", "state", "disease", "cases"]

# Filter to diseases of interest
diseases_of_interest = ["Diphtheria", "Tetanus", "Pertussis", "Measles", "Mumps", "Rubella"]
nndss_df = nndss_df[nndss_df["disease"].isin(diseases_of_interest)]

# Remove rows where state is US Territories or Totals
us_states = [
    'Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware',
    'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky',
    'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi',
    'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico',
    'New York', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania',
    'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont',
    'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming'
]
nndss_df = nndss_df[nndss_df["state"].isin(us_states)]

# Map state names to abbreviations
import us
state_abbrev = {state.name: state.abbr for state in us.states.STATES}
nndss_df["state"] = nndss_df["state"].map(state_abbrev)

# Convert case count to numeric
nndss_df["cases"] = pd.to_numeric(nndss_df["cases"], errors="coerce").fillna(0).astype(int)

# Final clean dataframe
nndss_cleaned = nndss_df[["state", "year", "disease", "cases"]]

# Save it
nndss_cleaned.to_csv("nndss_dtp_mmr_cleaned.csv", index=False)
print(" Cleaned and saved as 'nndss_dtp_mmr_cleaned.csv'")


import pandas as pd

# Load dataset
vacc = pd.read_csv("Vaccination_Coverage_and_Exemptions_among_Kindergartners_20250326.csv")

# Filter to only 'States' geography
vacc = vacc[vacc["Geography Type"] == "States"]

# Extract year from 'School Year' like "2021-22"
vacc["year"] = vacc["School Year"].str.extract(r"(\d{4})").astype(float)

# Filter to 2016–2022
vacc = vacc[vacc["year"].between(2016, 2022)]

# Clean and keep relevant columns
vacc = vacc[["Geography", "year", "Vaccine/Exemption", "Estimate (%)"]]
vacc.rename(columns={"Geography": "state"}, inplace=True)

# Convert Estimate to numeric (force errors to NaN)
vacc["Estimate (%)"] = pd.to_numeric(vacc["Estimate (%)"], errors="coerce")

# Pivot the table: One column per vaccine or exemption type
vacc_pivot = vacc.pivot_table(
    index=["state", "year"],
    columns="Vaccine/Exemption",
    values="Estimate (%)",
    aggfunc="mean"  # Handle duplicates
).reset_index()

# Save result
vacc_pivot.to_csv("cleaned_vaccination_data_2016_2022.csv", index=False)
print(" Vaccination data cleaned and saved.")


import pandas as pd

# -------------------------------
# STEP 1: Load datasets
# -------------------------------

vacc = pd.read_csv("cleaned_vaccination_data_2016_2022.csv")
nndss = pd.read_csv("nndss_dtp_mmr_cleaned.csv")
climate = pd.read_csv("climate_monthly_2016_2022.csv")
census = pd.read_csv("census_data_2016_2022.csv")

# -------------------------------
# STEP 2: Preprocess vaccination data
# -------------------------------

# Ensure correct types
vacc['year'] = vacc['year'].astype(int)

# Map state names to abbreviations
state_mapping = {
    'Alabama': 'AL', 'Alaska': 'AK', 'Arizona': 'AZ', 'Arkansas': 'AR', 'California': 'CA',
    'Colorado': 'CO', 'Connecticut': 'CT', 'Delaware': 'DE', 'Florida': 'FL', 'Georgia': 'GA',
    'Hawaii': 'HI', 'Idaho': 'ID', 'Illinois': 'IL', 'Indiana': 'IN', 'Iowa': 'IA',
    'Kansas': 'KS', 'Kentucky': 'KY', 'Louisiana': 'LA', 'Maine': 'ME', 'Maryland': 'MD',
    'Massachusetts': 'MA', 'Michigan': 'MI', 'Minnesota': 'MN', 'Mississippi': 'MS', 'Missouri': 'MO',
    'Montana': 'MT', 'Nebraska': 'NE', 'Nevada': 'NV', 'New Hampshire': 'NH', 'New Jersey': 'NJ',
    'New Mexico': 'NM', 'New York': 'NY', 'North Carolina': 'NC', 'North Dakota': 'ND', 'Ohio': 'OH',
    'Oklahoma': 'OK', 'Oregon': 'OR', 'Pennsylvania': 'PA', 'Rhode Island': 'RI', 'South Carolina': 'SC',
    'South Dakota': 'SD', 'Tennessee': 'TN', 'Texas': 'TX', 'Utah': 'UT', 'Vermont': 'VT',
    'Virginia': 'VA', 'Washington': 'WA', 'West Virginia': 'WV', 'Wisconsin': 'WI', 'Wyoming': 'WY'
}
vacc['state'] = vacc['state'].map(state_mapping)

# -------------------------------
# STEP 3: Preprocess NNDSS outbreak data
# -------------------------------

nndss['year'] = nndss['year'].astype(int)
nndss['state'] = nndss['state'].str.upper()

# Aggregate by state and year
nndss_grouped = (
    nndss.groupby(['state', 'year'])['cases']
    .sum()
    .reset_index()
)

# -------------------------------
# STEP 4: Preprocess climate data
# -------------------------------

climate['year'] = climate['year'].astype(int)
climate['state'] = climate['state'].str.upper()

# Aggregate monthly to yearly per state
climate_grouped = (
    climate.groupby(['state', 'year']).agg({
        'Precip_mm': 'mean',
        'Max_Temp_C': 'mean',
        'Min_Temp_C': 'mean'
    }).reset_index()
)

# -------------------------------
# STEP 5: Preprocess census data
# -------------------------------

census['year'] = census['year'].astype(int)
census['state'] = census['state'].str.upper()

# -------------------------------
# STEP 6: Merge all data
# -------------------------------

merged = climate_grouped.merge(census, on=['state', 'year'], how='left')
merged = merged.merge(vacc, on=['state', 'year'], how='left')
merged = merged.merge(nndss_grouped, on=['state', 'year'], how='left')

# Fill missing cases with 0
merged['cases'] = merged['cases'].fillna(0).astype(int)

# -------------------------------
# STEP 7: Export final dataset
# -------------------------------

merged.to_csv("final_model_dataset_2016_2022.csv", index=False)
print(" Final dataset saved as 'final_model_dataset_2016_2022.csv'")


merged
merged.isna().sum()

percent= merged.isnull().sum() * 100/len(merged)
percent

import pandas as pd

# Load your final dataset
df = pd.read_csv("final_model_dataset_2016_2022.csv")

# Columns to impute
vaccine_cols = [
    'DTP, DTaP, or DT', 'MMR', 'MMR (PAC)', 'Polio', 'Varicella',
    'Hepatitis B', 'Exemption'
]
numeric_cols = [
    'Precip_mm', 'Max_Temp_C', 'Min_Temp_C',
    'median_household_income', 'median_gross_rent', 'health_insurance_coverage'
]

cols_to_fill = vaccine_cols + numeric_cols

# 1. Create missing indicators for vaccine columns
for col in vaccine_cols:
    df[f"{col}_missing"] = df[col].isna().astype(int)

# 2. Fill missing values with median per state
df[cols_to_fill] = df.groupby('state')[cols_to_fill].transform(lambda x: x.fillna(x.median()))

# 3. Optional: Save the imputed dataset
df.to_csv("final_model_dataset_2016_2022_imputed_by_state.csv", index=False)
print(" Missing values filled per state and saved as 'final_model_dataset_2016_2022_imputed_by_state.csv'")


df





import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load all cleaned datasets
age_df = pd.read_csv("data/Cleaned_NNDSS_Age (1).csv")
ethnicity_df = pd.read_csv("data/Cleaned_NNDSS_Ethnicity (1).csv")
month_df = pd.read_csv("data/Cleaned_NNDSS_Month (1).csv")
race_df = pd.read_csv("data/Cleaned_NNDSS_Race (1).csv")
sex_df = pd.read_csv("data/Cleaned_NNDSS_Sex (1).csv")

# Display head of one sample dataframe
month_df.head()


# Convert to datetime if needed
month_df['Month'] = pd.to_datetime(month_df['Month'])

# Plot a few diseases over time
diseases = month_df['Disease'].unique()[:4]
plt.figure(figsize=(12, 6))
for disease in diseases:
    subset = month_df[month_df['Disease'] == disease]
    plt.plot(subset['Month'], subset['Cases'], label=disease)

plt.title("Monthly Case Trends for Selected Diseases")
plt.xlabel("Month")
plt.ylabel("Cases")
plt.legend()
plt.tight_layout()
plt.show()
