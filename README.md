# 🧬 Hepatitis B and C Virus Dashboard

This repository contains the code and files for the **interactive Hepatitis Genomics Dashboard**, which integrates HBV and HCV whole-genome sequences, genotype distributions, resistance mutations, recombination data, population denominators, and global burden estimates.

---

## 🔧 Installation

Clone the repository and set up the environment:

```bash
git clone <your-repo-url>
cd <your-repo-name>
conda env create -f requirements.yaml
conda activate hepatitis-dashboard
```

## 📊 Data Preparation

The dashboard integrates multiple datasets. Before running, prepare the required input files:

1. Retrieve HBV and HCV genomic data and metadata
   
  ```bash
  python hepatitis_data_collection.py
  ```

2. Multiple sequence alignment & phylogenetics

  - Perform MSA (e.g., using MAFFT)
  
  - Build maximum-likelihood trees (e.g., IQ-TREE)
  
  - Color tree tips by reference genotypes

4. Mutation analysis

  - Use the geno2pheno online tool

6. Recombination analysis
 
  - Use RDP5
  
  - Save recombinant accession IDs into text files

7. Population data
   
  ```bash
  python population_retriever.py
  ```

6. IHME Global Burden of Disease data

   - Download prevalence, incidence, and deaths for HBV and HCV from the IHME GBD dataset (https://vizhub.healthdata.org/gbd-results/)

## ▶️ Running the Dashboard

Once the data are prepared:

1. Generate dashboard-ready files:
   
   ```bash
   python generate_dash_data.py
   ```
3. Launch the app:
   
   ```bash
   python Hepatitis_dash.py
   ```

By default, the app runs locally at:

👉 http://127.0.0.1:8050

📂 Repository Structure

├── Hepatitis_Dash.py              # Main application file
├── requirements.txt               # Python dependencies
├── Procfile                      # For deployment (if needed)
├── runtime.txt                   # Python version specification
├── .gitignore
├── assets/
│   ├── style.css
│   ├── images/                   # Store images here
│   └── custom.js                # Any custom JavaScript
├── data/
│   ├── hepatitis_data_collection.py
│   ├── generate_dash_data.py
│   ├── population_retriever.py
│   ├── WHO_regions_countries_coordinates.txt
│   └── processed_data/           # For cached/processed data
├── pages/
│   ├── __init__.py              # Make pages a Python package
│   ├── dashboard.py             # Main dashboard page
│   ├── about.py                 # About page
│   ├── resources.py             # Resources page
│   └── contact.py               # Contact page
├── utils/
│   ├── __init__.py
│   ├── data_loader.py           # Move data_loader here
│   └── helpers.py               # Any helper functions
└── README.md
## 📖 Citation

If you use this dashboard or data in your work, please cite:

Tshiabuila D., Vagner D., Darren M., et al. Interactive Hepatitis B and C Virus Genomics Dashboard (2025).

