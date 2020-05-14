# AccretR: A planetary accretion and composition code in R

AccretR is a planetary body accretion Monte Carlo code that calculates mass, radius and bulk composition along a specified growth track, for orderly/hierarchical, runaway, and random particle accretion models. Elements in the model include concentrations of: H, C, N, O, Na, Mg, Al, Si, S, Cl, K, Ca, and Fe. Maximal water is also computed, assuming all H goes into forming water. Accretional heat is also calculated.

Currently, the code is optimized to build Jupiter's moon Europa, and Saturn's moons Titan and Enceladus, from CI, CM, CR, CK, CO and CV carbonaceous chondrite meteorites, cometary material (using comet 67P/Churyumov-Gerasimenko), and pure water ice. Literature sources for each of these compositions are found in the code.

Some care must be taken when selecting to accrete either Europa, Titan or Enceladus - simply comment out the unnecessary lines, and run in R. RStudio may be used.

To do: 
1) Timing of accretion needs to be tweaked to allow a change in radiative heat dissipation rates (i.e., shedding accretional heat).
2) Evaluation and quantification of melting and vaporization of maximal water reservoir after each particle impact (computationally intensive).
3) Various quality of life improvements, such as making the selection of bulk compositions and growth tracks more intuitive.

Desiderata:
1) Changing isotope reservoirs of stable isotopes.
2) A pretty GUI.
