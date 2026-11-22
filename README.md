# Harvest and Discards in the Recreational Gulf of Mexico Gag Fishery
 
Data and R code to reproduce the analysis in:
 
> Hyman, A. C., Chagaris, D., & Frazer, T. K. (2025). Influence of temporal regulations on harvest and discards in the recreational Gulf of Mexico gag (*Mycteroperca microlepis*) fishery. *Fisheries Research*, 285, 107332. https://doi.org/10.1016/j.fishres.2025.107332
 
A preprint of this work circulated as *Take it or leave it: influence of environmental variables and temporal regulations on harvest and discards in a multispecies fishery; a case study with gag (Mycteroperca microlepis)* ([SSRN](https://doi.org/10.2139/ssrn.5070788)). Some file headers in this repository still use that title.
 
## Overview
 
This study modeled monthly recreational gag harvest and discards along the west coast of Florida (Panhandle and Peninsula regions) as functions of temporal management, relative abundance, seasonality, environmental conditions, and angler effort. The fitted models were then combined with the effort model from the companion paper to simulate how different gag season start dates would affect harvest, landings, discards, and dead discards under a fixed annual catch target (ACT).
 
**Companion paper (effort model):** Hyman, A. C., Chagaris, D., Drexler, M., & Frazer, T. K. (2024). Modeling effort in a multispecies recreational fishery; influence of species-specific temporal closures, relative abundance, and seasonality on monthly angler-trips. *Fisheries Research*, 279, 107136. https://doi.org/10.1016/j.fishres.2024.107136 ([code](https://github.com/ChallenHymanPhD/Hyman-et-at-2024-Effort))
 
## Models
 
Monthly MRIP estimates of gag harvest ($H_t$) and discards ($D_t$), 2004 to 2023, are continuous, non-negative, and often zero when the season is closed. Both were modeled with Bayesian hurdle-Gamma regressions in [`brms`](https://paul-buerkner.github.io/brms/). For harvest:
 
$$
\begin{aligned}
H_t &\sim \text{HurdleGamma}(\lambda_t,\, \phi_t,\, \delta_t)\\
\log(\lambda_t) &= \omega_0 + \textstyle\sum_{u} e_{t,u}\,\omega_u && \text{mean of nonzero harvest}\\
\log(\phi_t) &= \nu_0 + \textstyle\sum_{j} w_{t,j}\,\nu_j && \text{Gamma shape}\\
\text{logit}(\delta_t) &= \gamma_0 + \textstyle\sum_{l} b_{t,l}\,\gamma_l && \text{probability of zero harvest}
\end{aligned}
$$
 
The discard model has the same structure. All terms interact with region, and all coefficients have $\text{Normal}(0, 10)$ priors.
 
| Component | Predictors |
|---|---|
| Mean (harvest and discards) | Prior-year gag CPUE index (log), juvenile gag index (log), annual and semi-annual harmonics, sea-surface temperature, fraction of month open to gag, fraction of month open to red snapper, effort (log angler-trips) |
| Shape (harvest and discards) | Harmonic terms, prior-year gag CPUE index (log) |
| Hurdle, harvest | Whether gag was open at any point in the month |
| Hurdle, discards | Prior-year gag CPUE index (log) |
 
Fifty randomly selected month-region observations were withheld for out-of-sample validation.
 
## Season simulations
 
The simulations propagate posterior uncertainty through three linked models: effort, harvest, and discards.
 
1. Draw expected effort from the [effort model](https://doi.org/10.1016/j.fishres.2024.107136) for a counterfactual design matrix.
2. Append those effort draws to the counterfactual data.
3. Draw expected harvest and discards from the hurdle-Gamma models.
4. Convert harvest to landings with mean monthly gag weight, and discards to dead discards with month- and region-specific release mortality from [Sauls (2014)](https://doi.org/10.1016/j.fishres.2013.10.008).
This is repeated for 1,000 posterior draws. Four 2024 scenarios were simulated: no gag season (baseline), and seasons opening June 1, September 1, or November 1. For each start date, season length starts at 60 days and is shortened one day at a time until median cumulative landings fall at or below an ACT of 635,364 lbs gutted weight. Non-management predictors are fixed at 2021–2023 region-month averages (juvenile index at its 2022 value), and red snapper and red grouper seasons are fixed at their 2024 Florida dates.
 
## Repository structure
 
| Branch | Contents |
|---|---|
| [`main`](https://github.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/tree/main) | This README |
| [`R-Files`](https://github.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/tree/R-Files) | Analysis code and the annotated supplement |
| [`Data-Files`](https://github.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/tree/Data-Files) | Model data, mortality estimates, and flow-chart figures |
 
### Code (`R-Files`)
 
| File | Description |
|---|---|
| `Hyman et al 2024 Supplemental Code github.Rmd` | Annotated supplement with full model equations, simulation procedure, and code |
| `Hyman-et-al-2024-Supplemental-Code.html` | Knitted version of the supplement ([view](https://htmlpreview.github.io/?https://github.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/blob/R-Files/Hyman-et-al-2024-Supplemental-Code.html)) |
| `Hyman et al 2024 Harvest All Analyses.R` | Standalone script that fits both models, runs the simulations, and produces the manuscript tables and figures |
| `Hyman et al 2024 Harvest and Discard helper functions.R` | Helper functions (season-to-month conversion, parameter tables, plot themes); sourced automatically by the scripts above |
 
### Data (`Data-Files`)
 
| File | Description |
|---|---|
| `MRIP Effort and Catch data.csv` | Monthly region-level model dataset (harvest, discards, effort, and all predictors). See the [effort repository](https://github.com/ChallenHymanPhD/Hyman-et-at-2024-Effort) for how it was built from raw MRIP data |
| `Sauls 2014 Mortality.csv` | Monthly gag release mortality by region, with confidence limits |
| `Hyman Flow Chart 1–3.png` | Conceptual diagrams of the simulation procedure |
 
The fitted effort model (`Effort_Model.rds`) is downloaded automatically from the effort repository.
 
## Reproducing the analysis
 
### Requirements
 
- R (≥ 4.1) with a working C++ toolchain for Stan
- R packages:
```r
install.packages(c("dplyr", "tidyr", "tidyverse", "lubridate", "ggplot2",
                   "ggrepel", "ggpubr", "RColorBrewer", "xtable", "kableExtra",
                   "readxl", "brms", "forecast"))
```
 
### Running
 
All data, helper functions, and the effort model are read directly from GitHub, so no manual downloads are needed.
 
1. Download `Hyman et al 2024 Harvest All Analyses.R` (or the `.Rmd`) from the `R-Files` branch.
2. Replace the `setwd()` path near the top with a local folder; downloaded files are written there.
3. Run the script, or knit the `.Rmd`.
Both models are fit with four chains on every run, and the season simulations loop over many candidate season lengths, so a full run takes some time.
 
## Citation
 
```bibtex
@article{Hyman2025Harvest,
  author  = {Hyman, A. Challen and Chagaris, David and Frazer, Thomas K.},
  title   = {Influence of temporal regulations on harvest and discards in the
             recreational {Gulf of Mexico} gag ({Mycteroperca microlepis}) fishery},
  journal = {Fisheries Research},
  year    = {2025},
  volume  = {285},
  pages   = {107332},
  doi     = {10.1016/j.fishres.2025.107332}
}
```
 
## Contact
 
A. Challen Hyman, PhD
University of South Florida, College of Marine Science / FWC Fish and Wildlife Research Institute
 
