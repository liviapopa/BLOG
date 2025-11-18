# BLOG
Associated code for the BLOG manuscript. This reposidory contains R scripts for **analyzing real metabolite-like data** and performing **univariate and multivariate simulations**, including validation and visualization. The workflow is modular, allowing users to run specific analyses or simulations independently.

## Repository Structure
| File | Description |
|------|-------------|
| `clean_real_data_rankings.R` | Computes Bayes factor rankings for features in a dataset using Zellner’s g-prior. Supports datasets with features × samples × timepoints. |
| `clean_simulations.R` | Runs sequential univariate and multivariate simulations for testing analysis pipelines. |
| `cleaned_simulations_parallel.R` | Parallelized version of simulation scripts for faster computation. |
| `clean_univariate_validation_plots.R` | Generates validation plots for univariate simulations. |
| `clean_multivariate_validation_plots.R` | Generates validation plots for multivariate simulations. |
| `clean_heatmaps_univariate_multivariate_simulation.R` | Produces heatmaps to visualize simulation results. |

## Requirements
The following R packages are required to run the scripts in this repository:

### Data Manipulation
- `dplyr`  
- `tidyr`  
- `purrr`  
- `tibble`  
- `reshape2`  
- `abind`  

### Plotting
- `ggplot2`  
- `cowplot`  
- `grid`  
- `plot.matrix`  

### Statistics / Modeling
- `MASS`  
- `MBSGS`  

### Parallel Computation
- `doParallel`  
- `foreach`  
### Installation

You can install all packages using:

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "tibble", "reshape2", "abind",
  "ggplot2", "cowplot", "plot.matrix", "MASS", "doParallel", "foreach"
))
```

# MBSGS is available on CRAN (or GitHub if necessary)
install.packages("MBSGS")

## Usage

### Real Data Analysis
1. Place your dataset in `.RData` or `.csv` format.
2. Edit `clean_real_data_rankings.R` to specify:
   - `data_path` – path to your dataset
   - `dat_array` – variable name containing a 3D array `[features x samples x timepoints]`
   - `output_file` – filename for Bayes factor rankings
3. Run the script:
```r
source("clean_real_data_rankings.R")
```

### Simulations
- Univariate and multivariate simulations can be run with `clean_simulations.R' or the parallel version `cleaned_simulations_parallel.R`.
- Validation plots and heatmaps are generated using the corresponding plotting scripts

## Output
- CSV files containing Bayes factor rankings for features.
- Validation plots and heatmaps for simulation studies.

## License
This repository is released under the [MIT License](LICENSE).

