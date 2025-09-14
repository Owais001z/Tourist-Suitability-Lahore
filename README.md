
# Tourist Suitability Analysis — Lahore, Pakistan

This project identifies candidate locations for a new tourist facility in Lahore (15 km buffer around city center), using a weighted-overlay suitability analysis.



## How to run
- Open `scripts/suitability_analysis_lahore.R` in R/RStudio.
- Install required packages: `sf`, `terra`, `osmdata`, `ggplot2`, `ggspatial`, `dplyr`, `scales`.
- Run the script to generate outputs.

## Data
All spatial features (roads, water, parks) are retrieved from OpenStreetMap using the `osmdata` package.
