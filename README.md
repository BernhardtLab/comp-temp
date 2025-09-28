---
editor_options: 
  markdown: 
    wrap: 72
---

# Data and code for Davis et al., General predictions for the effects of warming on competition

Description: In this study, we incorporated temperature sensitivity into
a MacArthur consumer-resource model with two consumers and two
resources. We incorporated temperature sensitivity using Arrhenius-style
temperature effects, and then we added empirically-derived
temperature-sensitivity estimates for each mechanistic process in the
consumer-resource model. We then simulated warming on species with
randomly drawn thermal traits (temperature sensitivities) to investigate
the effects of warming on competition. Raw temperature sensitivity data
are located in the repository's
[data/](https://github.com/BernhardtLab/comp-temp/tree/main/data)
folder. The workflow for the analysis includes 5 scripts, meant to be
run in order. Scripts are located in the repository's
[R-scripts/](https://github.com/BernhardtLab/comp-temp/tree/main/R-scripts)
folder and figures produced by each script are stored in the
[figures/](https://github.com/BernhardtLab/comp-temp/tree/main/figures)
folder.

**01-param-scripts.R**: This script processes synthesized published
estimates of temperature sensitivities for the processes underling
competition in the MacArthur consumer-resource model. The script
requires param-eas.csv as an input and then generates posterior
distributions for each parameter and generates the output
param_post_dists.csv, stored in the
[data/processed-data/](https://github.com/BernhardtLab/comp-temp/tree/main/data/processed-data)
folder. This script also generates main text figure 2.

**02-temp-dep-macarthur.R**: This script provides the function that
incorporates temperature sensitivity into the MacArthur consumer
resource model, as well as the Arrhenius function. This function is used
as an input for scripts 03, 04, and 05.

**03-var-ea-one-by-one.R**: This script analyses the effects of
temperature sensitivity in each parameter on niche and fitness
differences. Each parameter is investigated on its own, while all other
parameters are held invariant with temperature. This script requires
inputs from scripts 01 and 02, and produces main text figures 3 and 4,
as well as supplementary figures S2, S4, and S8.

**04-full-temp-var-analysis.R**: This script analyses the effects of
simultaneous temperature sensitivity in all parameter on niche and
fitness differences. This script requires inputs from scripts 01 and 02,
and produces main text figure 5, as well as supplementary figures S1,
S5, and S7.

**05-start-point-exercise.R**: This script analyses how different
simulation start points, imposed here by differences in the resource use
preferences of competing species pairs, influences the warming
trajectories of species. This script requires inputs from scripts 01 and
02, and produces supplementary figures S9 and S10.

If you have any questions about this repository, please direct them to
the manuscript's corresponding author, Kaleigh Davis, at
[kaleigh.davis\@uoguelph.ca](mailto:kaleigh.davis@uoguelph.ca){.email}
