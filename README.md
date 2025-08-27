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

- Hepatitis_dash.py – main dashboard application

- generate_dash_data.py – precomputes and formats dashboard data

- hepatitis_data_collection.py – fetches HBV/HCV genomic data and metadata

- population_retriever.py – retrieves country-level population denominators

- hbv_resistant_mutations.csv / hcv_resistance_mutations.csv – drug resistance mutation data

- hbv_burden_adjusted_coverage.csv / hcv_burden_adjusted_coverage.csv – burden-adjusted sequencing coverage estimates

- requirements.yaml – reproducible conda environment

## 📖 Citation

If you use this dashboard or data in your work, please cite:

Tshiabuila D., Vagner D., Darren M., et al. Interactive Hepatitis B and C Virus Genomics Dashboard (2025).

