

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
   install.packages(c("dismo", "stringr", "ggplot", "ggmap", "maps", "tidyverse", "readr", "magrittr", "rgbif", "taxize", "forcats", "sp", "raster", "ade4", "RColorBreewer", "ks", "data.table", "gridExtra", "viridis", "lsmeans", "nlme", "lme4", "lmerTest", "MASS", "DHARMa", "effectsize", "effects", "MuMIn", "parameters"))
   ```
 
## Script Descriptions

### Script 01 - Download Species Occurrences

- **Purpose**: Downloads species occurrences data.
- **Required files**: `niche_modelling/data/field_raw_data.csv`

### Script 02 - Species Occurrence Filtering

- **Purpose**: Filters the downloaded species occurrences.
- **Required files**:
  1. Climatic layers from CHELSA (see section 1.2 for download instructions).
  2. Dataframe containing all species occurrences (`./niche_modelling/occurrences/unclean/0061950-200613084148143.csv`) previously obtained in script 01 which must be directly downloaded from [GBIF](https://doi.org/10.15468/dl.9c6h5v).

### Script 03a - Multivariate Niche Modelling

- **Purpose**: Characterizes each species niches in the multivariate space.
- **Advice**: Run from section 4. The rest may take a few hours, especially on slow computers. Heavy files (`all_sp_climate.csv`) are not present and can be obtained running 02 script.
- **Required files**:
  1. `niche_modelling/data/field_raw_data.csv`
  2. `niche_mokdelling/occurrences/clean/all_sp_climate.csv`

### Script 03b - Univariate Niche Modelling

- **Purpose**: Characterizes each species niches in the univariate space.
- **Advice**: Run from section 3. The rest may take a few hours, especially on slow computers. Heavy files (`all_sp_climate.csv`) are not present and can be obtained running 02 script.
- **Required files**:
  1. `niche_modelling/data/field_raw_data.csv`
  2. `niche_modelling/data/populations_climate.csv`
  3. `niche_modelling/output/perc95_niche_chelsa_univariate.csv`
  4. `niche_mokdelling/occurrences/clean/all_sp_climate.csv`

### Script 04a - Multivariate CD Characterization

- **Purpose**: Characterizes centroid distances in a multivariate niche space.
- **Advice**: Only section 3 (for niche surface discarding lowest 95 percentile) is required to obtain tables for later analyses.
- **Required files**: 
  1. `niche_modelling/output/centroid_distances.csv`
  2. `niche_modelling/output/perc95_distances_in0.csv`
  3. `niche_modelling/output/perc90_distances_in0.csv`

### Script 04b - Univariate CD Characterization
- **Purpose**: Characterizes centroid distances in a univariate niche modelling context.
- **Advice**: Only section 3 (for niche surface discarding lowest 95 percentile) is required to obtain tables for later analyses. In addition, the script is prepared to run all biovariables, but only bio01, bio06 and bio12 are required for ulterior analyses.
- **Required files**: Various outputs from scripts 03b and 04a (see detailed list in the provided information).
- Based on centroid:
  1. `niche_modelling/Output/centroid_distances_univariate.csv` obtained from 03b
  2. `niche_modelling/Output/perc95_distances_in0.csv` obtained from 03b
  3. `statistical_analyses/Output/disequilibrium_nurse.csv` obtained from 04a
  4. `statistical_analyses/Output/disequilibrium_centroid_univariate.csv` obtained from 04b
- Based on mode:
  1. `niche_modelling/output/centroid_distances_univariate.csv` from 03b
  2. `niche_modelling/output/perc95_distances_in0.csv` from 03b
  3. `statistical_analyses/output/disequilibrium_mode_centroid.csv` from 04a
  4. `statistical_analyses/output/disequilibrium_mode_univariate.csv` from 04b
- Based on niche surface 95 perc:
  1. `niche_modelling/output/perc95_distances_in0_univariate.csv`from 03b
  2. `statistical_analyses/output/disequilibrium_perc95_in0.csv` from 04a
  3. `statistical_analyses/output/disequilibrium_perc_95_in0_univariate.csv` 04b
- Based on niche surface 90 perc:
  1. `niche_modelling/output/perc90_distances_in0_univariate.csv` from 03b
  2. `statistical_analyses/output/disequilibrium_perc90_in0.csv` from 04a
  3. `statistical_analyses/output/disequilibrium_perc_90_in0_univariate.csv` 04b

## Scripts 05 and 06 - Statistical Analyses

- **Purpose**: Performs statistical analyses to obtain main text figures 3 to 6.
- **Advice**: Both Scripts 05 and 06 can be completely run with these provided tables obtained from the rest of the R scripts. They allow obtaining main text figures 3 to 6.
- **Required files**:
  1. `./statistical_analyses/Output/disequilibrium_mode_centroid.csv`
  2. `./statistical_analyses/Output/disequilibrium_nurse.csv`
  3. `./statistical_analyses/Output/disequilibrium_perc90_in0.csv`
  4. `./statistical_analyses/Output/disequilibrium_perc95_in0.csv`
  5. `./statistical_analyses/Output/explanatory_vars.csv`
  6. `./statistical_analyses/Output/Facilitation_ClimaticDisequilibrium_chisq.csv`
  7. `./statistical_analyses/Output/change_abundances.csv`
  8. `./statistical_analyses/Output/disequilibrium_centroid_univariate.csv`
  9. `./statistical_analyses/Output/disequilibrium_mode_univariate.csv`
  10. `./statistical_analyses/Output/disequilibrium_perc_90_in0_univariate.csv`
  11. `./statistical_analyses/Output/disequilibrium_perc_95_in0_univariate.csv`

## Running the Scripts

Each script can be run independently, provided the required files are in place. It is advised to follow the script order as listed for a smooth workflow. Scripts 1 to 3 can be slow to run. Scripts 5 and 6 include the files and code required for the final analyses and allow to obtain the figures and results of the manuscript.

## License

This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0). This license allows others to share, copy and redistribute the material in any medium or format, and adapt, remix, transform, and build upon the material for any purpose. 

You can view a copy of this license at [http://creativecommons.org/licenses/by/4.0/](http://creativecommons.org/licenses/by/4.0/).

### How to Cite This Work

If you use this work, or any part of it, in your research or project, please provide appropriate credit to the authors and a link to the original source and the license. Here is a suggested format for citation:

- Perez-Navarro MA, Lloret F, Molina-Venegas R, Alcantara JM, Verdu M. (2024). Plant canopies promote climatic disequilibrium in Mediterranean recruit communities.





