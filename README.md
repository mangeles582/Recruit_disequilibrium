

# Plant canopies promote climatic disequilibrium in Mediterranean recruit communities

This repository contains a series of R scripts developed for the study titled "Plant canopies promote climatic disequilibrium in Mediterranean recruit communities", authored by: Perez-Navarro MA, Lloret F, Molina-Venegas R, Alcántara JM and Verdú M. The author of these scripts is Perez-Navarro MA. These scripts are used to download species occurrences, filter occurrences, perform multivariate and univariate niche modelling, characterize climatic disequilibrium, and conduct statistical analyses.

## Getting Started

Before running the scripts, ensure you have R and RStudio installed on your computer. Additionally, some scripts require specific data files and dependencies. Follow the setup instructions below to prepare your environment.

### Prerequisites

- R and RStudio
- Required R packages: (dismo, stringr,ggplot, ggmap, maps, tidyverse, readr, magrittr, rgbif, taxize, forcats, sp, raster, ade4, RColorBreewer, ks, data.table, gridExtra, viridis, lsmeans, nlme, lme4, lmerTest, MASS, DHARMa, effectsize, effects, MuMIn, parameters)
- Download and place climatic layers from CHELSA as described in section 1.2. A download script is included in the R script files.

### Installation

1. Clone this repository to your local machine:
   ```
   git clone https://github.com/yourusername/niche-modelling-research.git
   ```
2. Install required R packages by running:
   ```
   install.packages(c("dplyr", "sf", "raster"))
   ```
   (Adjust the package list as needed.)

## Script Descriptions

### Script 01 - Download Species Occurrences

- **Purpose**: Downloads species occurrences data.
- **Required to run**: `niche_modelling/data/field_raw_data.csv`

### Script 02 - Species Occurrence Filtering

- **Purpose**: Filters the downloaded species occurrences.
- **Required to run**:
  - Climatic layers from CHELSA (see section 1.2 for download instructions).
  - Each species dataframe from `./niche_modelling/occurrences/unclean/gbif_` or directly from `0061950-200613084148143.csv` which can be downloaded from [GBIF](https://doi.org/10.15468/dl.9c6h5v).

### Script 03a - Multivariate Niche Modelling

- **Advice**: Run from section 4. The rest may take a few hours, especially on slow computers. Heavy files (`dudi_sp` and `all_sp_climate`) are not present.
- **Required to run**:
  1. `field_raw_data`
  2. `Dudi_sp`
  3. `perc95_niche_chelsa_noseasonality`
  4. `all_sp_climate.csv`

### Script 03b - Univariate Niche Modelling

- **Advice**: Run from section 4. The rest may take a few hours, especially on slow computers.
- **Required to run**:
  1. `field_raw_data`
  2. `populations_climate`
  3. `perc95_niche_chelsa_univariate`
  4. `all_sp_climate.csv`

### Script 04a - Multivariate CD Characterization

- **Purpose**: Characterizes centroid distances in a multivariate niche modelling context.
- **Required tables**: 
  - `Niche_modelling/Output/centroid_distances.csv`
  - `Niche_modelling/Output/perc95_distances_in0.csv`
  - `Niche_modelling/Output/perc90_distances_in0.csv`

### Script 04b - Univariate CD Characterization

- **Purpose**: Characterizes centroid distances in a univariate niche modelling context.
- **Required tables**: Various outputs from scripts 03b and 04a (see detailed list in the provided information).

### Script 05 and 06 - Statistical Analyses

- **Purpose**: Performs statistical analyses to obtain main text figures 3 to 6.
- **Required tables**: Various outputs from previous scripts (see detailed list in the provided information).

## Running the Scripts

Each script can be run independently, provided the required files are in place. It is advised to follow the script order as listed for a smooth workflow.

## Contributing

Please read [CONTRIBUTING.md](https://github.com/yourusername/niche-modelling-research/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/yourusername/niche-modelling-research/LICENSE.md) file for details.

## Acknowledgments

- Mention any collaborators, data providers, or any other acknowledgments.

---

Make sure to replace placeholders (like `https://github.com/yourusername/niche-modelling-research.git`) with actual links and details relevant to your project. This README template is designed to be comprehensive and user-friendly,
