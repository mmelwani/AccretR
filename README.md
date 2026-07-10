# AccretR: A planetary accretion and composition code in R

AccretR is a planetary body accretion Monte Carlo code that calculates mass, radius and bulk composition along a specified growth track, for orderly/hierarchical, runaway, and random particle accretion models. Elements in the model include concentrations of: H, C, N, O, Na, Mg, Al, Si, S, Cl, K, Ca, and Fe. Maximal water is also computed, assuming all H goes into forming water. Accretional heat is also calculated.

Currently, the code is set up to build Jupiter's moon Europa, Saturn's moons Titan and Enceladus, the dwarf planet Ceres, and generic ocean-world/telluric-world super-Earths, from CI, CM, CR, CK, CO and CV carbonaceous chondrite meteorites, cometary material (using comet 67P/Churyumov-Gerasimenko), and pure water ice. Literature sources for each of these compositions are found in the code.

`AccretR()` is a parameterized function — no source editing is required to change scenarios. Source `AccretR.R` in R (or RStudio) and call it with the arguments you want, e.g.:

```r
source("AccretR.R")
result <- AccretR(body = "Europa",
                   material_list = list(list("CI", CI_composition), list("CM", CM_composition), list("Comet", Comet_67P)),
                   material_prob = c(6.94e-3, 1.73, runif(1, 1.37e-3, 6.94e-3)),
                   bootstrap_n = 100)
```

Calling `AccretR()` with no arguments builds the default "OceanWorld" scenario. See the comment block above the `AccretR` function definition in `AccretR.R` for the full list of arguments (body preset, material mixture and probabilities, particle-size function, seed radius, nebula temperature, bootstrap count, core count, output file). The code can run on a HPC, and has been run successfully on the Texas Advanced Computing Center (TACC) Stampede2 HPC (https://www.tacc.utexas.edu/systems/stampede2).

To do: 
1) Timing of accretion needs to be tweaked to allow a change in radiative heat dissipation rates (i.e., shedding accretional heat).
2) Evaluation and quantification of melting and vaporization of maximal water reservoir after each particle impact (computationally intensive).
3) Rewrite as a R notebook.
5) Change the multi-plotting functions to tidyverse plots.
6) Adding new constraints on accretion material hydration and temperature derived from literature (e.g. hydrodynamic escape constraints, Bierson & Nimmo (2020). ApJ, 897(2), L43. https://doi.org/10.3847/2041-8213/aba11a).

Desiderata:
1) Changing isotope reservoirs of stable isotopes.
2) A pretty GUI.
