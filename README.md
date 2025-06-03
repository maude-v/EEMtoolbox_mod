# EEMtoolbox modification and application
This package provides a set of functions that propose modifications of the EEMtoolbox package to make its use more flexible and realistic for ecological and evolutionary modeling.
The modifications include the possibility to simulate systems with more than one non-native species, to set bounds on the systems species equilibrium abundances,
and to model systems where species introductions can lead to extinction of native species.
It includes functions to set abundance bounds on the native system EEMs, and provide the possibility to model systems where species introductions can lead to extinction of native species.
Artificial recruitment as well as species control and removal can also be simulated.

## Features
- **Multiple non-native species**: Simulate systems with more than one non-native species.
- **Equilibrium bounds**: Set bounds on the equilibrium abundances of species in the system.
- **Native species extinction**: Model systems where the introduction of non-native species can lead to the extinction of native species.
- **Artificial recruitment**: Simulate artificial recruitment of species.
- **Species control and removal**: Model species control and removal scenarios.

## Adapted functions
- `EEMtoolbox::adapted_dicrepancy_continuous_sum()`: Adapted from `EEMtoolbox::dicrepancy_continuous_sum()` to set bounds on equilibrium abundances.
- `add_species_names()`: Add names to the species in the system.
- `select_EEM_outputs()`: Check if the system' is at equilibrium and set bounds on the equilibrium abundances.'s equilibrium abundances are correct, and their disprition within the bounds.
- `merge_introductions()`: Allow for multiple non-native species introductions.
- `add_introduced_species()`: Add non-native species to the system, with no stability or feasability expectations. Some species may go extinct.
- `adapted_calculate_projections()`: Adapted from `EEMtoolbox::calculate_projections()`to allow for the artificial control of introduced and/or native species.
- `adapted_ode_solve()`: Adapted from `EEMtoolbox::ode_solve()` which is called inside `adapted_calculate_projections()`.
- `adapted_plot_projections()`: Slight adaptation of `EEMtoolbox::plot_projections()` to plot the results of the projections.

## Use
Users need to make sure that their version of R is at least 4.3.1. We recommend running the following code in R to preinstall all the necessary packages:
``` r
packages_to_install <- c("deSolve", "doParallel", "doSNOW","dplyr", "foreach", "ggplot2", "magrittr", "MASS", "nleqslv", "parallel", "parallelly","stats","tidyr")

# Install packages if not already installed
install_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Apply the function to install packages
invisible(lapply(packages_to_install, install_if_not_installed))
```

To install the modified EEMtoolbox package, download the repository as a zip file.
Rename the file as EEMtoolbox_mod.zip, then run the following code
``` r
devtools::install_local("path_to_file/EEMtoolbox_mod.zip", repos = NULL, type = "win.binary")
```

After installation, load the package with:
``` r
library(EEMtoolbox_mod)
```

Then, source the adapted functions that are not in the package:
``` r
source("path_to_file/adapted functions/add_species_names.R")
source("path_to_file/adapted functions/select_EEM_outputs.R")
source("path_to_file/adapted functions/merge_introductions.R")
source("path_to_file/adapted functions/add_introduced_species.R")
source("path_to_file/adapted functions/adapted_calculate_projections.R")
source("path_to_file/adapted functions/adapted_ode_solve.R")
source("path_to_file/adapted functions/adapted_plot_projections.R")
```

## Example usage
Open and run the following RMarkdown files:
``` r
file.edit("path_to_file/Application/pre_introduction/initial abundances/initial abundances.Rmd")
file.edit("path_to_file/Application/pre_introduction/pre_introduction_system.Rmd")
file.edit("path_to_file/Application/palm_control/palm_control.Rmd")
file.edit("path_to_file/Application/post_introduction/post_introduction_system.Rmd")
```

# Getting help
To seek support, report bugs or give suggestions, contact Maude Vernet at: maude.vernet@unibe.ch
