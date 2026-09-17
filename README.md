# AccretR: A planetary accretion and composition code in R

**Current release: v2.0.0 "Banana"**

AccretR is a planetary body accretion Monte Carlo code that calculates mass, radius and bulk composition along a specified growth track, for orderly/hierarchical, runaway, and random particle accretion models. Elements in the model include concentrations of: H, C, N, O, Na, Mg, Al, Si, S, Cl, K, Ca, and Fe. Maximal water is also computed, assuming all H goes into forming water. Accretional heat is also calculated.

Currently, the code is set up to build Jupiter's moon Europa, Saturn's moons Titan and Enceladus, the dwarf planet Ceres, and generic ocean-world/telluric-world super-Earths, from CI, CM, CR, CK, CO and CV carbonaceous chondrite meteorites, cometary material (using comet 67P/Churyumov-Gerasimenko), and pure water ice. Literature sources for each of these compositions are found in the code.

## Requirements

R packages: `foreach`, `iterators`, `parallel`, `doParallel`, `plyr`, `ggplot2`, `scales`.

```r
install.packages(c("foreach", "iterators", "doParallel", "plyr", "ggplot2", "scales"))
```

## Usage

As of v2.0.0, `AccretR()` is a parameterized function — no source editing is required to change scenarios (earlier versions required commenting/uncommenting blocks of code to configure a run). Source `AccretR.R` in R (or RStudio) and call it with the arguments you want:

```r
source("AccretR.R")

# Default scenario: generic 50/50 water-ice/CI-chondrite "OceanWorld", 100 bootstrap runs.
result <- AccretR()

# Europa built from a custom mix of CI + CM chondrites + comet 67P material,
# with mixture probabilities from Desch et al. (2018), ApJ.
result <- AccretR(body = "Europa",
                   material_list = list(list("CI", CI_composition), list("CM", CM_composition), list("Comet", Comet_67P)),
                   material_prob = c(6.94e-3, 1.73, runif(1, 1.37e-3, 6.94e-3)),
                   bootstrap_n = 100)
```

`AccretR()` returns a named list of summary statistics (median/mean/standard deviation for elemental wt. %, maximal H2O wt. %, particle count, body radius/mass/density, accretion energy, and surface temperature), and produces two ggplot2 figures (composition histograms and body property histograms) as a side effect. By default it also writes a text summary of the result to `AccretR_output_<body>.txt`.

### Arguments

| Argument | Description | Default |
|---|---|---|
| `body` | Preset name: one of `"Europa"`, `"Titan"`, `"Enceladus"`, `"Ceres"`, `"OceanWorld"`, `"TelluricWorld"`. Each preset bundles a growth track and a termination radius so the two can't drift out of sync. | `"OceanWorld"` |
| `material_list` | List of `list(name, composition_vector)` building blocks to accrete from (e.g. `CI_composition`, `CM_composition`, `Comet_67P`, `Water_ice` — see `AccretR.R` for the full set and their literature sources). | 50/50 water ice + CI chondrite |
| `material_prob` | Sampling probability for each entry of `material_list`. Must be the same length as `material_list`; normalized automatically. | `c(0.5, 0.5)` |
| `seed_radius_m` | Starting ("seed") body radius, in meters. | `500` |
| `T_nebula` | Nebula/circumplanetary disk temperature, in Kelvin, used by the Lunine & Stevenson (1982) surface temperature model. | `30` |
| `particle_radius_fn` | `function(total_body_radius)` returning one accreting particle's radius in meters. Default is "Type D" (100 km embryos ± 50%, appropriate for super-Earths). Other options given in `AccretR.R`: Type A (fastest, runaway growth proportional to body radius), Type B (Galilean-satellite-scale impactors), Type C (very slow, dehydrated Europa-building pebbles — best run on an HPC). | Type D |
| `bootstrap_n` | Number of full-body Monte Carlo realizations to run for statistics. | `100` |
| `cores` | Number of parallel workers to use. | all detected cores minus 1 |
| `out_file` | Path to write a text summary of the result list. `NULL` auto-names it from `body`; `NA` skips writing a file entirely. | `NULL` |

Passing an unknown `body` name, or a `material_prob` whose length doesn't match `material_list`, raises an immediate error naming the problem.

The code can run on a HPC, and has been run successfully on the Texas Advanced Computing Center (TACC) Stampede2 HPC (https://www.tacc.utexas.edu/systems/stampede2).

To do: 
1) Timing of accretion needs to be tweaked to allow a change in radiative heat dissipation rates (i.e., shedding accretional heat).
2) Evaluation and quantification of melting and vaporization of maximal water reservoir after each particle impact (computationally intensive).
3) Rewrite as a R notebook.
5) Change the multi-plotting functions to tidyverse plots.
6) Adding new constraints on accretion material hydration and temperature derived from literature (e.g. hydrodynamic escape constraints, Bierson & Nimmo (2020). ApJ, 897(2), L43. https://doi.org/10.3847/2041-8213/aba11a).

Desiderata:
1) Changing isotope reservoirs of stable isotopes.
2) A pretty GUI.
