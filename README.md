# Hepatitis-spatial-transcriptomics

## Liver Spatial Transcriptomics Analysis Pipeline

This repository contains a comprehensive Python pipeline for analyzing spatial transcriptomics data from liver samples with three different conditions: Autoimmune Hepatitis (AIH), Donor (D), and Seronegative hepatitis (SN). The pipeline performs neighborhood enrichment analysis and Ripley's spatial statistics to understand cellular spatial organization patterns.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Data Structure](#data-structure)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)
- [Analysis Pipeline](#analysis-pipeline)
- [Output](#output)

## Overview

This pipeline analyzes spatial transcriptomics data from liver samples to:
- Perform neighborhood enrichment analysis at multiple spatial radii (30, 50, 80 microns)
- Calculate Ripley's G function for spatial clustering analysis
- Compare spatial patterns between different liver disease conditions
- Generate statistical comparisons between conditions using Kruskal-Wallis and Mann-Whitney U tests

The analysis focuses on 12 distinct cell types:
- B cell
- Epithelial
- Hepatic stellate
- Hepatocyte
- Macrophage
- Monocyte
- NK cell
- Neutrophil
- Plasma cell
- T cell
- Type 1 LSEC (Liver Sinusoidal Endothelial Cell)
- Type 2 LSEC

## Features

- **Neighborhood Enrichment Analysis**: Quantifies spatial co-localization patterns between cell types
- **Multi-radius Analysis**: Analyzes spatial patterns at 30μm, 50μm, and 80μm radii
- **Ripley's G Function**: Statistical analysis of spatial clustering patterns
- **Cross-condition Comparison**: Statistical comparison of spatial patterns between disease conditions
- **Field-of-view (FOV) Analysis**: Per-FOV analysis for robust statistical testing
- **Comprehensive Visualization**: Heatmaps, boxplots, and Ripley function plots
- **Parenchyma vs Non-Parenchyma Analysis**: Specialized analysis for parenchyma/non-parenchyma cells

## Requirements

### Python Version
- Python 3.10.18 (packaged by conda-forge)

### Core Dependencies
- scanpy==1.11.3
- squidpy==1.6.5
- pandas==2.3.1
- numpy==2.2.6
- matplotlib>=3.5.0
- seaborn>=0.11.0
- scipy>=1.9.0
- pointpats>=2.2.0
- pathlib (built-in)

## Installation

1. **Clone the repository:**
   ```bash
   git clone <https://github.com/kyliesavoye/hepatitis-spatial-transcriptomics.git>
   cd liver-spatial-analysis
   ```

2. **Create a conda environment:**
   ```bash
   conda create -n spatial_liver python=3.10
   conda activate spatial_liver
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

   Or install manually:
   ```bash
   pip install scanpy==1.11.3 squidpy==1.6.5 pandas==2.3.1 numpy==2.2.6
   pip install matplotlib seaborn scipy pointpats
   ```

## Data Structure

The pipeline expects the following directory structure:

```
"/path/to/your/data/"
├── counts/
│   ├── AIH_1607_expression_mat_fixed.csv
│   ├── AIH_6079_expression_mat_fixed.csv
│   ├── D_6414_expression_mat_fixed.csv
│   ├── D_6446_expression_mat_fixed.csv
│   ├── D_7678_expression_mat_fixed.csv
│   ├── SN_606_expression_mat_fixed.csv
│   ├── SN_1042_expression_mat_fixed.csv
│   └── SN_2739_expression_mat_fixed.csv
├── meta/
│   ├── AIH_1607_spatial_metadata_fixed.csv
│   ├── AIH_6079_spatial_metadata_fixed.csv
│   ├── D_6414_spatial_metadata_fixed.csv
│   ├── D_6446_spatial_metadata_fixed.csv
│   ├── D_7678_spatial_metadata_fixed.csv
│   ├── SN_606_spatial_metadata_fixed.csv
│   ├── SN_1042_spatial_metadata_fixed.csv
│   └── SN_2739_spatial_metadata_fixed.csv
├── AIH_output/
├── D_output/
├── SN_output/
└── comparison_output/
```

**Note**: You will need to update the `root` path in `funcs.py` or the `get_paths()` function to match your local data directory structure.

## Usage

### Individual Condition Analysis

Run analysis for each condition separately:

```bash
# Analyze AIH samples
python Liver_analysis_AIH.py

# Analyze D samples  
python Liver_analysis_D.py

# Analyze SN samples
python Liver_analysis_SN.py
```

### Cross-condition Comparison

After running individual analyses, perform comparative analysis:

```bash
# Compare all conditions
python Liver_analysis_comparison.py
```

### Customizing Analysis Parameters

To modify analysis parameters, edit the relevant functions in `funcs.py`:

- **Spatial radius**: Modify `radius` parameter in `neighbourhood_enrichment()` calls
- **Cell type colors**: Update the `palette` list in each main script
- **Output directories**: Modify paths in `get_paths()` function

## File Descriptions

### Core Analysis Files

- **`funcs.py`**: Central module containing all analysis functions and utilities
  - `get_paths()`: Manages dataset-specific file paths
  - `neighbourhood_enrichment()`: Performs spatial neighborhood enrichment analysis
  - `run_neighbourhood_enrichment_per_fov()`: FOV-level neighborhood analysis
  - `average_neighbourhood_plot()`: Creates averaged heatmaps across samples
  - `plot_ripley_all_clusters()`: Ripley's G function analysis and plotting
  - `canonical_pair()`: Standardizes cell type pair naming

### Condition-Specific Analysis Scripts

- **`amberLiverDA_AIH.py`**: Autoimmune Hepatitis sample analysis
  - Analyzes 2 AIH samples (1607, 6079)
  - Generates neighborhood enrichment heatmaps at multiple radii
  - Computes Ripley's G statistics with random simulations
  - Performs parenchyma vs non-parenchyma cells analysis

- **`amberLiverDA_D.py`**: Steatotic liver disease sample analysis  
  - Analyzes 3 D samples (6414, 6446, 7678)
  - Similar analysis pipeline as AIH
  - Includes immune cell subset analysis

- **`amberLiverDA_SN.py`**: Steatohepatitis sample analysis
  - Analyzes 3 SN samples (606, 1042, 2739)
  - Comprehensive spatial analysis pipeline
  - Parenchyma/non-parenchyma cells specific analysis

### Comparison Analysis

- **`amberLiverDA_comparison.py`**: Cross-condition comparative analysis
  - Loads results from individual condition analyses  
  - Generates statistical comparisons using Kruskal-Wallis and Mann-Whitney U tests
  - Creates comparative boxplots for key cell type interactions
  - Focuses on macrophage-related interactions

### Function-Specific Modules

- **`funcs_AIH.py`**: AIH-specific analysis functions
- **`funcs_D.py`**: D condition-specific functions  
- **`funcs_SN.py`**: SN condition-specific functions
- **`funcs_comparison.py`**: Comparison analysis utilities

## Analysis Pipeline

### 1. Data Loading and Preprocessing
- Load expression matrices and spatial metadata using Squidpy
- Convert cell type annotations to categorical format
- Quality control and data validation

### 2. Neighborhood Enrichment Analysis
- Calculate spatial neighbors at multiple radii (30μm, 50μm, 80μm)
- Compute enrichment Z-scores for all cell type pairs
- Generate condition-specific and averaged heatmaps

### 3. Spatial Statistics
- Perform Ripley's G function analysis for each cell type
- Compare observed patterns against random simulations (n=99)
- Analyze both all cell types and immune cell subsets

### 4. Field-of-View Analysis  
- Run per-FOV neighborhood enrichment for statistical robustness
- Enable sample-level and FOV-level comparisons
- Support for meta-analysis across conditions

### 5. Statistical Comparison
- Kruskal-Wallis test for overall condition differences
- Mann-Whitney U test for pairwise condition comparisons
- Focus on clinically relevant cell type interactions

## Output

The pipeline generates several types of output files:

### Neighborhood Enrichment
- `{condition}_neighborhood_enrichment_average_radius{X}.png`: Averaged heatmaps
- `{condition}_neighborhood_enrichment_average_radius{X}.csv`: Z-score matrices
- Individual sample heatmaps for each radius

### Ripley's Analysis  
- `{sample}_ripley_G_sq.png`: Individual Ripley's G plots with confidence intervals
- `{condition}_average_ripley_G_function_with_random_expectation.png`: Averaged plots
- `{condition}_average_ripley_G_immunecells.png`: Immune cell focused plots

### Comparative Analysis
- `macrophage_cell_pairs_boxplots_with_stats.png`: Statistical comparison boxplots
- Per-FOV Z-score data: `{condition}_combined_zscore_per_fov.csv`

### Parenchyma vs Non-parenchyma Analysis
- `{sample}_ripley_G_sq_PvsNP.png`: parenchyma/non-parenchyma cell specific spatial analysis

## Data Format Requirements

### Expression Matrix (CSV)
- Rows: Cells (with unique cell IDs)
- Columns: Genes/features  
- Index: Cell barcodes/IDs

### Spatial Metadata (CSV)
Required columns:
- `final_cell_types`: Cell type annotations
- `fov`: Field of view identifier  
- `P_VS_NP`: parenchyma vs non-parenchyma classification
- Spatial coordinates (x, y coordinates)

## Troubleshooting

### Common Issues

1. **Path Errors**: Update the `root` path in `funcs.py` to match your data location
2. **Memory Issues**: For large datasets, consider processing fewer samples simultaneously
3. **Missing Dependencies**: Ensure all required packages are installed with correct versions
4. **Data Format**: Verify that CSV files have the expected column names and format

### Performance Optimization

- Use SSD storage for faster I/O operations
- Consider increasing available RAM for large datasets
- Parallelize FOV processing if needed

---

**Version Information:**
- Python: 3.10.18
- Last updated: September 2025
- Tested on: Linux/Unix systems
