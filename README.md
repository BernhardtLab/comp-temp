# Data and code for Davis et al., General predictions for the effects of warming on competition

Authors: Kaleigh E. Davis, Tess N. Grainger, Po-Ju Ke, Patrick L. Thompson, Mary I. O'Connor and Joey R. Bernhardt

Description: In this study, we incorporated temperature sensitivity into a MacArthur consumer-resource model with two consumers and two resources. We incorporated temperature sensitivity using Arrhenius-style temperature effects, and then we added empirically-derived temperature-sensitivity estimates for each mechanistic process in the consumer-resource model. We then simulated warming on species with randomly drawn thermal traits (temperature sensitivities) to investigate the effects of warming on competition. Raw temperature sensitivity data are located in the repository's [data/](https://github.com/BernhardtLab/comp-temp/tree/main/data) folder. The workflow for the analysis includes 7 scripts, meant to be run in order. Scripts are located in the repository's [R-scripts/](https://github.com/BernhardtLab/comp-temp/tree/main/R-scripts) folder and figures produced by each script are stored in the [figures/](https://github.com/BernhardtLab/comp-temp/tree/main/figures) folder.

## Scripts

**01-param-dists.R**: This script processes synthesized published estimates of temperature sensitivities for the processes underling competition in the MacArthur consumer-resource model. The script requires param-eas.csv as an input and then generates posterior distributions for each parameter and generates the output param_post_dists.csv, stored in the [data/processed-data/](https://github.com/BernhardtLab/comp-temp/tree/main/data/processed-data) folder. This .csv file is a required input for scripts 03, 04, 05, 06, and 07. This script also generates main text figure 2.

**02-temp-dep-macarthur.R**: This script provides the function that incorporates temperature sensitivity into the MacArthur consumer resource model. This function is used as an input for scripts 04-07.

**03-arrhenius.R**: This script provides the function that defines how model parameters respond to temperature. This function is used as an input for scripts 04-07.

**04-var-ea-one-by-one.R**: This script analyses the effects of temperature sensitivity in each parameter on niche and fitness differences. Each parameter is investigated on its own, while all other parameters are held invariant with temperature. This script requires inputs from scripts 02 and 03, and produces main text figures 3 and 4, as well as supplementary figures S6 and S9.

**05-full-temp-var-analysis.R**: This script analyses the effects of simultaneous temperature sensitivity in all parameters on niche and fitness differences. This script requires inputs from scripts 02 and 03, and produces main text figure 5, as well as supplementary figures S5, S10-S12.

**06-start-point-exercise.R**: This script analyses how different simulation start points, imposed here by differences in the resource use preferences of competing species pairs, influences the warming trajectories of species. This script requires inputs from scripts 02 and 03, and produces supplementary figures S2-4 and S8.

**07-different-CR-models.R**: This script analyses how warming affects competition when we consider consumer-resource models with more resource diversity (8 species) and unimodal temperature responses. This script requires inputs from scripts 02 and 03, as well as supplemental scripts 07S-temp-dep-macarthur-8r.R and 07S-temp-dep-macarthur-unimodal.R. This script produces Figures S1 and S13.

For scripts 02 - 07, temperature sensitivities are defined in each simulation, with the naming convention: "{parameter_EAik}", where ik captures the relevant consumer, resource, or both. The parameter value at the reference (ambient) temperature in each simulation is given following the naming convention: "{parameter-ik_b}". Consumer species are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively. The set of parameters at the ambient temperature define the model's starting conditions, and thus the position in ND--FD space of the all species pairs for a given simulation prior to warming.

Scripts only used for supplementary analyses are contained in the folder R-scripts \> supplemental-analysis and are arranged based on the main script that they are related to. Script 01S-param-dists-univsmulti.R explores the representation of unicellular and multicellular organisms in our synthesized dataset and generates posterior distributions for each paramter based on cellularity, where possible. Scripts 07S-temp-dep-macarthur-8r.R and 07S-temp-dep-macarthur-unimodal.R define the temperature-dependent competition model and calculations of niche and fitness differences for a consumer-resource model with 8 resources and a two-resource consumer-resource model where parameters respond to warming with unimodal curves, defined by the Johnson-Lewin model (Johnson and Lewin 1946, Journal of Cellular and Comparative Physiology).

## Data

The [data/raw-data/](https://github.com/BernhardtLab/comp-temp/tree/main/data/raw-data) folder contains two files:

-   **param-eas.csv**: Synthesized published estimates of temperature sensitivities (activation energies) for the processes underlying competition in the MacArthur consumer-resource model. This is the primary input for script 01.
-   **params-ea-metadata.csv**: Column-level metadata for param-eas.csv, describing each field including parameter groupings, specific metrics, literature references, experiment type, taxonomic groupings, and the method used to estimate activation energies.

## Interactive Shiny App

An interactive version of the MacArthur consumer-resource model is available as a Shiny app, allowing users to explore model behaviour without running the R scripts directly. The app is hosted at:

<https://joeybernhardt.shinyapps.io/macarthur-consumer-resource-model/>

The source code for the app is located in the [R-scripts/MacArthur-shiny/](https://github.com/BernhardtLab/comp-temp/tree/main/R-scripts/MacArthur-shiny) folder.

## Requirements

The analysis scripts require R and the following packages (versions used in the manuscript are noted):

-   `tidyverse` (1.3.2)
-   `MCMCpack` (1.7.0)
-   `MCMCvis` (0.16.3)
-   `bayesplot` (1.11.1)
-   `patchwork` (1.3.2)
-   `cowplot` (1.1.1)
-   `colorspace` (2.1.0)
-   `viridis` (0.6.5)
-   `janitor` (2.1.0)
-   `purrr` (1.0.2)
-   `car` (3.1.2)
-   `see` (0.13.0)
-   `visreg` (2.7.0)
-   `beepr` (2.0)

Packages can be installed via `install.packages()`. The analysis was developed and tested in R; users are advised to use a recent version of R (≥ 4.0) for compatibility.

## Reproducibility

To reproduce the exact results from the manuscript, package versions matter — particularly for the MCMC-based script (01-param-dists.R), where differences in `MCMCpack` or `MCMCvis` versions can affect posterior distributions. To check your current package environment, run the following in R:

``` r
sessionInfo()
```

If you encounter differences in results, please check your package versions against those listed in `sessionInfo()` output and consider opening an issue on the repository.

## Citation

If you use this code or data, please cite the associated manuscript:

> Davis, K.E., Grainger, T.N., Ke, P.-J., Thompson, P.L., O'Connor, M.I., and Bernhardt, J.R. General predictions for the effects of warming on competition. Ecology Letters, 2026. DOI: [DOI]

## License

This repository is released under the [MIT License](LICENSE).

## Contact

If you have any questions about this repository, please direct them to the manuscript's corresponding author, Kaleigh Davis, at [kaleigh.davis\@uoguelph.ca](mailto:kaleigh.davis@uoguelph.ca){.email}
