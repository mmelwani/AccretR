# AccretR, a program to calculate the composition of an icy ocean world using the building blocks of our solar system.
# Copyright (C) 2022 Mohit Melwani Daswani

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Author contact details:
# Mohit Melwani Daswani
# melwani.mohit@gmail.com
# Latest update: 31st October 2022.

# Invariant physical constants:
Earth_mass_kg <- 5.972e+24
Earth_radius_m <- 6371e+3
G <- 6.67408e-11
H2O_LHV <- 2260e+3 # J/kg Latent heat of vaporization of water (Datt et al., 2014)
H2O_LHF <- 334e+3 # J/kg Latent heat of fusion of water (Datt et al., 2014)
H2O_LHS <- 2838e+3 # J/kg Latent heat of sublimation of water (Datt et al., 2014)
H2O_Cp <- 4184 # J/Kg/K Specific heat capacity of water
H_molar_mass <- 1.00794 # g/mol
O_molar_mass <- 15.999 # g/mol
H2O_molar_mass <- 2*H_molar_mass + O_molar_mass

# Compositions of the building materials in weight percent. Compositions from Lodders and Fegley (1998) for CM, CV, CO, CK and CR chondrites, Palme et al. (2014) for CI chondrites, Clay et al. (2017) for chlorine in all meteorites (average CV chlorine used for CK chondrites); all normalized to 100 wt. %. H was not reported in Lodders and Fegley (1998) for CK and CR chondrites, so taken from CV chondrites for CK chondrites. For CR chondrites, H is estimated from DTA/DTG for 10 CR chondrites in Garenne et al (2014), where  mass lost between 200 - 400 C is assumed to be H2O from oxyhydroxides and mass lost from 400 - 770 is OH- from phyllosilicates. Bulk densities from Flynn et al. (2018). Heat capacities approximated from Ostrowski and Bryson (2019) at 200 K. Comet 67P/Churyumov-Gerasimenko composition is a synthesis of Pätzold et al. (2016), Dhooghe et al. (2017), Le Roy et al. (2015), Bardyn et al. (2017). Heat capacity of comet 67P/C-G adopted from Hu et al. (2017). Comet 67P/C-G density and dust-to-ice ratio of 4 (by mass) (already taken into account in composition) is from Pätzold et al. (2016). In parentheses: (H wt. %, C wt. %, Mg wt. %, Al wt. %, Si wt. %, S wt. %, Ca wt. %, Fe wt. %, O wt. %, Na wt. %, K wt. %, Cl wt. %, N wt. %, density in kg/m3, heat capacity in J/(kg*K^-1)). Pass any combination of these into AccretR()'s material_list argument to build a body out of them.
CI_composition <- c(2.006,3.543,9.714,0.855,10.895,5.448,0.928,19.000,46.737,0.505,0.057,0.012,0.300,1570,500)
CM_composition <- c(1.428,2.244,11.733,1.153,12.957,2.755,1.316,21.731,44.073,0.398,0.038,0.019,0.155,2270,500)
CV_composition <- c(0.287,0.544,14.679,1.725,16.116,2.258,1.889,24.123,37.980,0.349,0.037,0.005,0.008,2970,500)
Comet_67P <- c(11.283,26.814,0.985,0.178,10.434,1.789,0.082,6.035,40.652,0.696,0.031,0.001,1.021,533,1000)
CO_composition <- c(0.071,0.447,14.727,1.422,16.047,2.234,1.605,25.391,37.578,0.427,0.037,0.006,0.009,3100,500)
CK_composition <- c(0.285,0.224,14.942,1.494,16.060,1.728,1.728,23.379,39.802,0.315,0.029,0.005,0.008,2900,500)
CR_composition <- c(0.319,2.038,13.961,1.172,15.286,1.936,1.315,24.254,39.280,0.336,0.032,0.007,0.063,3110,500)
Water_ice <- c(11.191,0,0,0,0,0,0,0,88.809,0,0,0,0,916.9,1800)

# Growth-track and termination-radius presets, keyed by body name (Sotin et al. 2007, Icarus, modified to fit specific bodies; termination radii from the JPL SSD Planetary Satellites Physical Parameters site). Pass one of these names as AccretR()'s `body` argument: growth track and termination radius are looked up together so they can't get out of sync with each other.
AccretR_body_presets <- list(
  "Europa"        = list(growth_track = function(mass_ratio) 1.0720*mass_ratio^0.306, termination_radius_m = 1560800),
  "Titan"         = list(growth_track = function(mass_ratio) 1.2705*mass_ratio^0.302, termination_radius_m = 2574730),
  "Enceladus"     = list(growth_track = function(mass_ratio) 1.0709*mass_ratio^0.302, termination_radius_m = 252100),
  "Ceres"         = list(growth_track = function(mass_ratio) 1.0754*mass_ratio^0.306, termination_radius_m = 469600),
  "OceanWorld"    = list(growth_track = function(mass_ratio) 1.2620*mass_ratio^0.275, termination_radius_m = Earth_radius_m*2), # 50 % water, M/M_Earth = 1 to 10.
  "TelluricWorld" = list(growth_track = function(mass_ratio) 1.0000*mass_ratio^0.274, termination_radius_m = Earth_radius_m*2)  # 0.05 % water, M/M_Earth = 1 to 10.
)

# AccretR() builds a body particle-by-particle via Monte Carlo accretion, bootstrapping `bootstrap_n` full-body realizations to report summary statistics and histograms.
#
# Arguments:
#   body               Name of a preset in AccretR_body_presets (growth track + termination radius, kept in sync automatically). One of "Europa", "Titan", "Enceladus", "Ceres", "OceanWorld", "TelluricWorld". Default "OceanWorld".
#   material_list      List of list(name, composition_vector) building blocks to accrete from. Default: 50/50 water ice + CI chondrite.
#                       Example (Europa built from CI + CM chondrites + comet 67P):
#                         material_list = list(list("CI", CI_composition), list("CM", CM_composition), list("Comet", Comet_67P))
#   material_prob      Sampling probability for each entry of material_list (must be the same length as material_list; normalized automatically by sample()).
#                       Example (Desch et al. 2018, ApJ, distance-based probabilities for the material_list above):
#                         material_prob = c(6.94e-3, 1.73, runif(1, 1.37e-3, 6.94e-3))
#   seed_radius_m      Starting ("seed") body radius in meters.
#   T_nebula           Nebula/circumplanetary disk temperature in Kelvin, used by the Lunine & Stevenson (1982) surface temperature model.
#   particle_radius_fn function(total_body_radius) returning one particle radius in meters for the next accretion step. Default is Type D (100 km embryos +- 50 %, appropriate for super-Earths). Other options:
#                         - Type A (fastest to run): runaway growth proportional to body radius: function(r) runif(1, 0.1*r, 0.5*r)
#                         - Type B (slow): Galilean-satellite-scale impactors (Barr & Canup 2008; Canup & Ward 2009): function(r) runif(1, 1, 1000)
#                         - Type C (very slow, best run on an HPC): dehydrated Europa-building pebbles (Ronnet et al. 2017): function(r) runif(1, 0.005, 0.5)
#   bootstrap_n        Number of full-body Monte Carlo realizations to run for statistics.
#   cores              Number of parallel workers to use.
#   out_file           Path to write a text summary of the result list. NULL (default) auto-names it from `body`; NA skips writing a file entirely.
AccretR <- function(body = "OceanWorld",
                     material_list = list(list("Ice", Water_ice), list("CI", CI_composition)),
                     material_prob = c(0.5, 0.5),
                     seed_radius_m = 500,
                     T_nebula = 30,
                     particle_radius_fn = function(total_body_radius) runif(n = 1, min = 0.5, max = 1.5)*100000,
                     bootstrap_n = 100,
                     cores = max(1, parallel::detectCores() - 1),
                     out_file = NULL){

  # Load the package to parallelize
  library(foreach)
  library(iterators)
  library(parallel)
  library(doParallel)
  library(plyr)
  library(ggplot2)
  library(scales)

  # Start the timer
  time_stamp <- proc.time()

  if (!body %in% names(AccretR_body_presets)) {
    stop("Unknown body '", body, "'. Choose one of: ", paste(names(AccretR_body_presets), collapse = ", "))
  }
  if (length(material_prob) != length(material_list)) {
    stop("material_prob must have the same length as material_list (", length(material_list), " entries).")
  }
  body_preset <- AccretR_body_presets[[body]]

  # Number of workers/cores used for parallelization.
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  AccretR_main_subroutine <- function(){
    # Initial system values
    H_mass_kg <- 0
    C_mass_kg <- 0
    Mg_mass_kg <- 0
    Al_mass_kg <- 0
    Si_mass_kg <- 0
    S_mass_kg <- 0
    Ca_mass_kg <- 0
    Fe_mass_kg <- 0
    O_mass_kg <- 0
    Na_mass_kg <- 0
    K_mass_kg <- 0
    Cl_mass_kg <- 0
    N_mass_kg <- 0
    total_body_mass <- 0
    N_particles <- 0
    E_accretion <- 0
    Delta_T <- 0
    total_body_heat_capacity <- 0
    total_body_radius <- seed_radius_m #Start with a seed. For an icy moon, 10 m is appropriate. For a planet, <273 K temperatures are maintained until planet radius is = 500 m (Solomon 1979, Toksoz et al. 1978, Hanks & Anderson 1969, McSween et al. 2002).
    Temperature_diff <- 0
    Temperature_LS1982 <- T_nebula #Temperature of the body surface from Lunine and Stevenson (1982), including latent heat loss from Lunine and Stevenson (1982).

    repeat{
      # Radius of the accreting particle, in meters.
      particle_radius <- particle_radius_fn(total_body_radius)

      # Volume of the accreting particle, in m^3.
      particle_volume <- 4/3*pi*(particle_radius^3)

	  #########################################

      # Select material class based on material_prob. Normalization of the probabilities to 1 is done automatically by the "sample" function.
      select_material <- sample(material_list, 1, replace = TRUE, prob = material_prob)

	    unlist_vector <- unlist(select_material)

	  #########################################

      # Obtain mass of the accreting particle, in kg, by multiplying particle volume by material class density:
      particle_mass <- (particle_volume)*(as.numeric(unlist_vector[[15]]))

      # Obtain total mass of elements:
      H_mass_kg <- H_mass_kg + (particle_mass*as.numeric(unlist_vector[[2]]))
      C_mass_kg <- C_mass_kg + (particle_mass*as.numeric(unlist_vector[[3]]))
      Mg_mass_kg <- Mg_mass_kg + (particle_mass*as.numeric(unlist_vector[[4]]))
      Al_mass_kg <- Al_mass_kg + (particle_mass*as.numeric(unlist_vector[[5]]))
      Si_mass_kg <- Si_mass_kg + (particle_mass*as.numeric(unlist_vector[[6]]))
      S_mass_kg <- S_mass_kg + (particle_mass*as.numeric(unlist_vector[[7]]))
      Ca_mass_kg <- Ca_mass_kg + (particle_mass*as.numeric(unlist_vector[[8]]))
      Fe_mass_kg <- Fe_mass_kg + (particle_mass*as.numeric(unlist_vector[[9]]))
      O_mass_kg <- O_mass_kg + (particle_mass*as.numeric(unlist_vector[[10]]))
      Na_mass_kg <- Na_mass_kg + (particle_mass*as.numeric(unlist_vector[[11]]))
      K_mass_kg <- K_mass_kg + (particle_mass*as.numeric(unlist_vector[[12]]))
      Cl_mass_kg <- Cl_mass_kg + (particle_mass*as.numeric(unlist_vector[[13]]))
      N_mass_kg <- N_mass_kg + (particle_mass*as.numeric(unlist_vector[[14]]))

      # Normalize composition to 100 wt. %
      Sum_composition <- sum(H_mass_kg, C_mass_kg, Mg_mass_kg, Al_mass_kg, Si_mass_kg, S_mass_kg, Ca_mass_kg, Fe_mass_kg, O_mass_kg, Na_mass_kg, K_mass_kg, Cl_mass_kg, N_mass_kg)
      H_wt_perc <- 100*H_mass_kg/Sum_composition
      C_wt_perc <- 100*C_mass_kg/Sum_composition
      Mg_wt_perc <- 100*Mg_mass_kg/Sum_composition
      Al_wt_perc <- 100*Al_mass_kg/Sum_composition
      Si_wt_perc <- 100*Si_mass_kg/Sum_composition
      S_wt_perc <- 100*S_mass_kg/Sum_composition
      Ca_wt_perc <- 100*Ca_mass_kg/Sum_composition
      Fe_wt_perc <- 100*Fe_mass_kg/Sum_composition
      O_wt_perc <- 100*O_mass_kg/Sum_composition
      Na_wt_perc <- 100*Na_mass_kg/Sum_composition
      K_wt_perc <- 100*K_mass_kg/Sum_composition
      Cl_wt_perc <- 100*Cl_mass_kg/Sum_composition
      N_wt_perc <- 100*N_mass_kg/Sum_composition

	  #########################################

      # Heat capacity equation (must come BEFORE mass, radius and growth track equations)
      particle_heat_capacity <- as.numeric(unlist_vector[[16]])
      total_body_heat_capacity <- (total_body_heat_capacity*total_body_mass + particle_heat_capacity*particle_mass)/(total_body_mass + particle_mass)

      # Calculate total body mass and radius according to the selected body's growth track (given in Earth radii and masses, from Sotin et al. 2007 (Icarus), modified to fit specific bodies).
      total_body_mass <- total_body_mass + particle_mass
      total_body_radius <- body_preset$growth_track(total_body_mass/Earth_mass_kg) * Earth_radius_m
      # Density of the body:
      total_body_bulk_density <- total_body_mass/(4/3*pi*(total_body_radius^3))

      # Energy and temperature equations (must come AFTER mass, radius and growth track equations). E_accretion and Temperature_diff is the energy of accretion, and temperature difference between the body's surface and the accretion disk, modified from Stevenson et al. (1986) Origins of Satellites chapter in Satellites book. Also in other places (e.g., Johansen et al. 2022, A&A Eq. 19).
      E_accretion <- E_accretion + G*total_body_mass/total_body_radius
      Temperature_diff <- Temperature_diff + (G*total_body_mass)/(total_body_radius*total_body_heat_capacity)

      #An alternative calculation of the surface temperature of the accreting body can be found in Bierson and Nimmo (2020), based on Kuramoto and Matsui (1994).
      #Temperature_BN2020 <-

      #Yet another alternative calculation of the surface temperature, from Lunine and Stevenson (1982), Equation 25, with evaporative and sublimation cooling.
      Particle_H2O_mass_frac <- (((particle_mass*as.numeric(unlist_vector[[2]])/100)*(1/H_molar_mass)*(1/2)*(H2O_molar_mass))/total_body_mass)
      if (Temperature_LS1982>=273.15) {
        Temperature_LS1982 <- Temperature_LS1982 + (T_nebula + 2*pi/3*G*total_body_bulk_density*(total_body_radius^2)/total_body_heat_capacity) - (Particle_H2O_mass_frac*H2O_LHV/H2O_Cp)
     } else {
        Temperature_LS1982 <- max(Temperature_LS1982, T_nebula) + (2*pi/3*G*total_body_bulk_density*(total_body_radius^2)/total_body_heat_capacity) - (Particle_H2O_mass_frac*H2O_LHS/H2O_Cp)
     }


      # Total number of accreted particles
      N_particles=N_particles+1

      # Exit once the body has grown to the target radius for the selected body preset.
      if (total_body_radius>=body_preset$termination_radius_m) break
    }
    return(list("H wt. %" = H_wt_perc, "C wt. %" = C_wt_perc, "Mg wt. %" = Mg_wt_perc, "Al wt. %" = Al_wt_perc, "Si wt. %" = Si_wt_perc, "S wt. %" = S_wt_perc, "Ca wt. %" = Ca_wt_perc, "Fe wt. %" = Fe_wt_perc, "O wt. %" = O_wt_perc, "Na wt. %" = Na_wt_perc, "K wt. %" = K_wt_perc, "Cl wt. %" = Cl_wt_perc, "N wt. %" = N_wt_perc, "body radius (m)" = total_body_radius, "body mass (kg)" = total_body_mass, "body bulk density (kg/m^3)" = total_body_bulk_density, "Energy of accretion (J/(kg*K^-1))" = E_accretion, "Accretion disk T (K)" = Temperature_diff, "Body surface T, with evaporative cooling (K)" = Temperature_LS1982, "Number of particles"=N_particles))
  }

  # Number of full bodies to build for statistics. Repeat the main subroutine bootstrap_n times and output all results into a transposed dataframe.
  Total_bootstrap_run <- foreach(i=1:bootstrap_n) %dopar% {AccretR_main_subroutine()}
  Total_bootstrap_run_frame <- as.data.frame((Total_bootstrap_run))

  # Output medians and standard deviations of the data gathered in the data frame
  N_particles_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Number.of.particles", colnames(Total_bootstrap_run_frame))]))
  N_particles_median <- median(N_particles_all_runs)
  N_particles_mean <- mean(N_particles_all_runs)
  N_particles_sd <- sd(N_particles_all_runs)
  H_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("H.wt.", colnames(Total_bootstrap_run_frame))]))
  H_wt_perc_median <- median(H_wt_perc_all_runs)
  H_wt_perc_mean <- mean(H_wt_perc_all_runs)
  H_wt_perc_sd <- sd(H_wt_perc_all_runs)
  C_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("C.wt.", colnames(Total_bootstrap_run_frame))]))
  C_wt_perc_median <- median(C_wt_perc_all_runs)
  C_wt_perc_mean <- mean(C_wt_perc_all_runs)
  C_wt_perc_sd <- sd(C_wt_perc_all_runs)
  Mg_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Mg.wt.", colnames(Total_bootstrap_run_frame))]))
  Mg_wt_perc_median <- median(Mg_wt_perc_all_runs)
  Mg_wt_perc_mean <- mean(Mg_wt_perc_all_runs)
  Mg_wt_perc_sd <- sd(Mg_wt_perc_all_runs)
  Al_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Al.wt.", colnames(Total_bootstrap_run_frame))]))
  Al_wt_perc_median <- median(Al_wt_perc_all_runs)
  Al_wt_perc_mean <- mean(Al_wt_perc_all_runs)
  Al_wt_perc_sd <- sd(Al_wt_perc_all_runs)
  Si_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Si.wt.", colnames(Total_bootstrap_run_frame))]))
  Si_wt_perc_median <- median(Si_wt_perc_all_runs)
  Si_wt_perc_mean <- mean(Si_wt_perc_all_runs)
  Si_wt_perc_sd <- sd(Si_wt_perc_all_runs)
  S_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("S.wt.", colnames(Total_bootstrap_run_frame))]))
  S_wt_perc_median <- median(S_wt_perc_all_runs)
  S_wt_perc_mean <- mean(S_wt_perc_all_runs)
  S_wt_perc_sd <- sd(S_wt_perc_all_runs)
  Ca_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Ca.wt.", colnames(Total_bootstrap_run_frame))]))
  Ca_wt_perc_median <- median(Ca_wt_perc_all_runs)
  Ca_wt_perc_mean <- mean(Ca_wt_perc_all_runs)
  Ca_wt_perc_sd <- sd(Ca_wt_perc_all_runs)
  Fe_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Fe.wt.", colnames(Total_bootstrap_run_frame))]))
  Fe_wt_perc_median <- median(Fe_wt_perc_all_runs)
  Fe_wt_perc_mean <- mean(Fe_wt_perc_all_runs)
  Fe_wt_perc_sd <- sd(Fe_wt_perc_all_runs)
  O_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("O.wt.", colnames(Total_bootstrap_run_frame))]))
  O_wt_perc_median <- median(O_wt_perc_all_runs)
  O_wt_perc_mean <- mean(O_wt_perc_all_runs)
  O_wt_perc_sd <- sd(O_wt_perc_all_runs)
  Na_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Na.wt.", colnames(Total_bootstrap_run_frame))]))
  Na_wt_perc_median <- median(Na_wt_perc_all_runs)
  Na_wt_perc_mean <- mean(Na_wt_perc_all_runs)
  Na_wt_perc_sd <- sd(Na_wt_perc_all_runs)
  K_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("K.wt.", colnames(Total_bootstrap_run_frame))]))
  K_wt_perc_median <- median(K_wt_perc_all_runs)
  K_wt_perc_mean <- mean(K_wt_perc_all_runs)
  K_wt_perc_sd <- sd(K_wt_perc_all_runs)
  Cl_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Cl.wt.", colnames(Total_bootstrap_run_frame))]))
  Cl_wt_perc_median <- median(Cl_wt_perc_all_runs)
  Cl_wt_perc_mean <- mean(Cl_wt_perc_all_runs)
  Cl_wt_perc_sd <- sd(Cl_wt_perc_all_runs)
  N_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("N.wt.", colnames(Total_bootstrap_run_frame))]))
  N_wt_perc_median <- median(N_wt_perc_all_runs)
  N_wt_perc_mean <- mean(N_wt_perc_all_runs)
  N_wt_perc_sd <- sd(N_wt_perc_all_runs)
  Total_body_radius_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("body.radius", colnames(Total_bootstrap_run_frame))]))
  Total_body_radius_median <- median(Total_body_radius_all_runs)
  Total_body_radius_mean <- mean(Total_body_radius_all_runs)
  Total_body_radius_sd <- sd(Total_body_radius_all_runs)
  Total_body_density_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("body.bulk.density", colnames(Total_bootstrap_run_frame))]))
  Total_body_density_median <- median(Total_body_density_all_runs)
  Total_body_density_mean <- mean(Total_body_density_all_runs)
  Total_body_density_sd <- sd(Total_body_density_all_runs)
  Total_body_mass_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("body.mass", colnames(Total_bootstrap_run_frame))]))
  Total_body_mass_median <- median(Total_body_mass_all_runs)
  Total_body_mass_mean <- mean(Total_body_mass_all_runs)
  Total_body_mass_sd <- sd(Total_body_mass_all_runs)
  Accretion_energy_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Energy.of.accretion", colnames(Total_bootstrap_run_frame))]))
  Accretion_energy_median <- median(Accretion_energy_all_runs)
  Accretion_energy_mean <- mean(Accretion_energy_all_runs)
  Accretion_energy_sd <- sd(Accretion_energy_all_runs)
  Temperature_diff_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Accretion.disk", colnames(Total_bootstrap_run_frame))]))
  Temperature_diff_median <- median(Temperature_diff_all_runs)
  Temperature_diff_mean <- mean(Temperature_diff_all_runs)
  Temperature_diff_sd <- sd(Temperature_diff_all_runs)
  Temperature_LS1982_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("with.evaporative", colnames(Total_bootstrap_run_frame))]))
  Temperature_LS1982_median <- median(Temperature_LS1982_all_runs)
  Temperature_LS1982_mean <- mean(Temperature_LS1982_all_runs)
  Temperature_LS1982_sd <- sd(Temperature_LS1982_all_runs)
  #Assuming all H is in water, calculate all water available
  Max_water_wt_perc_all_runs <- H_wt_perc_all_runs*(1/H_molar_mass)*(1/2)*(H2O_molar_mass)
  Max_water_wt_perc_median <- median(Max_water_wt_perc_all_runs)
  Max_water_wt_perc_mean <- mean(Max_water_wt_perc_all_runs)
  Max_water_wt_perc_sd <- sd(Max_water_wt_perc_all_runs)

  # Output the type of parallelization used, and the number of workers
  Workers_used <- c(getDoParName(), getDoParWorkers())

  # Plot histograms

  # Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
      print(plots[[1]])
      return(invisible(plots[[1]]))

    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }

      # Capture the assembled page as a grob so it can be reprinted or saved
      # (e.g. via ggsave) after this function returns, instead of being lost
      # as soon as the graphics device closes.
      return(invisible(grid.grab()))
    }
  }

  # AccretR histogram plots sensu strictu
  N_particles_all_runs_frame <- ldply(N_particles_all_runs, data.frame)
  N_particles_all_runs_plot <<- ggplot(N_particles_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    scale_x_continuous(labels = scientific) +
    xlab("Number of accreting planetesimals") +
    theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Total_body_radius_all_runs_frame <- ldply(Total_body_radius_all_runs, data.frame)
  Total_body_radius_all_runs_plot <<- ggplot(Total_body_radius_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    scale_x_continuous(labels = scientific) + xlab("Total body radius (m)") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Total_body_mass_all_runs_frame <- ldply(Total_body_mass_all_runs, data.frame)
  Total_body_mass_all_runs_plot <<- ggplot(Total_body_mass_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    scale_x_continuous(labels = scientific) + xlab("Total body mass (kg)") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Total_body_density_all_runs_frame <- ldply(Total_body_density_all_runs, data.frame)
  Total_body_density_all_runs_plot <<- ggplot(Total_body_density_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab((expression(paste("Bulk body densities (",kg/m^3,")")))) + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Accretion_energy_all_runs_frame <- ldply(Accretion_energy_all_runs, data.frame)
  Accretion_energy_all_runs_plot <<- ggplot(Accretion_energy_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    scale_x_continuous(labels = scientific) + xlab("Accretion energies (J)") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Temperature_diff_all_runs_frame <- ldply(Temperature_diff_all_runs, data.frame)
  Temperature_diff_all_runs_plot <<- ggplot(Temperature_diff_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    scale_x_continuous(labels = scientific) +
    xlab((expression(paste("Maximum surface-disk ", Delta,"T (K)")))) +
    theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Temperature_LS1982_all_runs_frame <- ldply(Temperature_LS1982_all_runs, data.frame)
  Temperature_LS1982_all_runs_plot <<- ggplot(Temperature_LS1982_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    scale_x_continuous(labels = scientific) +
    xlab(expression(paste("Surface T, ", H[2],O~italic("L"), " cooling (K)"))) +
    theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  AccretR_body_result_plot <<- multiplot(N_particles_all_runs_plot, Total_body_radius_all_runs_plot, Total_body_mass_all_runs_plot, Total_body_density_all_runs_plot, Accretion_energy_all_runs_plot, Temperature_diff_all_runs_plot, Temperature_LS1982_all_runs_plot, cols=2)
  ggsave("AccretR_body_result_plot.png", plot = AccretR_body_result_plot, width = 10, height = 12, dpi = 300)

  H_wt_perc_all_runs_frame <- ldply(H_wt_perc_all_runs, data.frame)
  H_wt_perc_all_runs_plot <<- ggplot(H_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk H wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  C_wt_perc_all_runs_frame <- ldply(C_wt_perc_all_runs, data.frame)
  C_wt_perc_all_runs_plot <<- ggplot(C_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk C wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Mg_wt_perc_all_runs_frame <- ldply(Mg_wt_perc_all_runs, data.frame)
  Mg_wt_perc_all_runs_plot <<- ggplot(Mg_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Mg wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Al_wt_perc_all_runs_frame <- ldply(Al_wt_perc_all_runs, data.frame)
  Al_wt_perc_all_runs_plot <<- ggplot(Al_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Al wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Si_wt_perc_all_runs_frame <- ldply(Si_wt_perc_all_runs, data.frame)
  Si_wt_perc_all_runs_plot <<- ggplot(Si_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Si wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  S_wt_perc_all_runs_frame <- ldply(S_wt_perc_all_runs, data.frame)
  S_wt_perc_all_runs_plot <<- ggplot(S_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +
    geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk S wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Ca_wt_perc_all_runs_frame <- ldply(Ca_wt_perc_all_runs, data.frame)
  Ca_wt_perc_all_runs_plot <<- ggplot(Ca_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Ca wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Fe_wt_perc_all_runs_frame <- ldply(Fe_wt_perc_all_runs, data.frame)
  Fe_wt_perc_all_runs_plot <<- ggplot(Fe_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Fe wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  O_wt_perc_all_runs_frame <- ldply(O_wt_perc_all_runs, data.frame)
  O_wt_perc_all_runs_plot <<- ggplot(O_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk O wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Na_wt_perc_all_runs_frame <- ldply(Na_wt_perc_all_runs, data.frame)
  Na_wt_perc_all_runs_plot <<- ggplot(Na_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Na wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  K_wt_perc_all_runs_frame <- ldply(K_wt_perc_all_runs, data.frame)
  K_wt_perc_all_runs_plot <<- ggplot(K_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk K wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Cl_wt_perc_all_runs_frame <- ldply(Cl_wt_perc_all_runs, data.frame)
  Cl_wt_perc_all_runs_plot <<- ggplot(Cl_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk Cl wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  N_wt_perc_all_runs_frame <- ldply(N_wt_perc_all_runs, data.frame)
  N_wt_perc_all_runs_plot <<- ggplot(N_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab("Bulk N wt. %") + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  Max_water_wt_perc_all_runs_frame <- ldply(Max_water_wt_perc_all_runs, data.frame)
  Max_water_wt_perc_all_runs_plot <<- ggplot(Max_water_wt_perc_all_runs_frame, aes(x=X..i..)) +
    geom_histogram(bins = 15) +  geom_vline(aes(xintercept=median(X..i..)), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=mean(X..i..)), color="blue", linetype="dotdash") +
    xlab((expression(paste("Maximal ", H[2],"O wt. %")))) + theme_bw() + labs(caption = "Red dashed = median; blue dot-dash = mean")

  AccretR_composition_result_plot <<- multiplot(H_wt_perc_all_runs_plot, C_wt_perc_all_runs_plot, Mg_wt_perc_all_runs_plot, Al_wt_perc_all_runs_plot, Si_wt_perc_all_runs_plot, S_wt_perc_all_runs_plot, Ca_wt_perc_all_runs_plot, Fe_wt_perc_all_runs_plot, O_wt_perc_all_runs_plot, Na_wt_perc_all_runs_plot, K_wt_perc_all_runs_plot, Cl_wt_perc_all_runs_plot, N_wt_perc_all_runs_plot,Max_water_wt_perc_all_runs_plot, cols=3)
  ggsave("AccretR_composition_result_plot.png", plot = AccretR_composition_result_plot, width = 14, height = 16, dpi = 300)

  # Free up the cores in the cluster
  stopCluster(cl)

  # Return results
  AccretR_result <<- (list("Median H wt. %" = H_wt_perc_median, "Mean H wt. %" = H_wt_perc_mean, "Standard deviation H wt. %" = H_wt_perc_sd, "Median C wt. %" = C_wt_perc_median, "Mean C wt. %" = C_wt_perc_mean, "Standard deviation C wt. %" = C_wt_perc_sd, "Median Mg wt. %" = Mg_wt_perc_median, "Mean Mg wt. %" = Mg_wt_perc_mean, "Standard deviation Mg wt. %" = Mg_wt_perc_sd, "Median Al wt. %" = Al_wt_perc_median, "Mean Al wt. %" = Al_wt_perc_mean, "Standard deviation Al wt. %" = Al_wt_perc_sd, "Median Si wt. %" = Si_wt_perc_median, "Mean Si wt. %" = Si_wt_perc_mean, "Standard deviation Si wt. %" = Si_wt_perc_sd, "Median S wt. %" = S_wt_perc_median, "Mean S wt. %" = S_wt_perc_mean, "Standard deviation S wt. %" = S_wt_perc_sd, "Median Ca wt. %" = Ca_wt_perc_median, "Mean Ca wt. %" = Ca_wt_perc_mean, "Standard deviation Ca wt. %" = Ca_wt_perc_sd, "Median Fe wt. %" = Fe_wt_perc_median, "Mean Fe wt. %" = Fe_wt_perc_mean, "Standard deviation Fe wt. %" = Fe_wt_perc_sd, "Median O wt. %" = O_wt_perc_median, "Mean O wt. %" = O_wt_perc_mean, "Standard deviation O wt. %" = O_wt_perc_sd, "Median Na wt. %" = Na_wt_perc_median, "Mean Na wt. %" = Na_wt_perc_mean, "Standard deviation Na wt. %" = Na_wt_perc_sd, "Median K wt. %" = K_wt_perc_median, "Mean K wt. %" = K_wt_perc_mean, "Standard deviation K wt. %" = K_wt_perc_sd, "Median Cl wt. %" = Cl_wt_perc_median, "Mean Cl wt. %" = Cl_wt_perc_mean, "Standard deviation Cl wt. %" = Cl_wt_perc_sd, "Median N wt. %" = N_wt_perc_median, "Mean N wt. %" = N_wt_perc_mean, "Standard deviation N wt. %" = N_wt_perc_sd, "Median maximum H2O wt. %" = Max_water_wt_perc_median, "Mean maximum H2O wt. %" = Max_water_wt_perc_mean, "Standard deviation maximum H2O wt. %" = Max_water_wt_perc_sd, "Median number of particles" = N_particles_median, "Mean number of particles" = N_particles_mean, "Standard deviation number of particles" = N_particles_sd, "Median total body radius (m)" = Total_body_radius_median, "Mean total body radius (m)" = Total_body_radius_mean, "Standard deviation total body radius (m)" = Total_body_radius_sd, "Median total body mass (kg)" = Total_body_mass_median, "Mean total body mass (kg)" = Total_body_mass_mean, "Standard deviation total body mass (kg)" = Total_body_mass_sd, "Median body bulk density" = Total_body_density_median, "Mean body bulk density" = Total_body_density_mean, "Standard deviation body bulk density" = Total_body_density_sd, "Median accretion energy (J)" = Accretion_energy_median, "Mean accretion energy (J)" = Accretion_energy_mean, "Standard deviation accretion energy (J)" = Accretion_energy_sd, "Median temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_median, "Mean temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_mean, "Standard deviation temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_sd, "Median surface temperature with evaorative cooling (K)" = Temperature_LS1982_median, "Mean surface temperature with evaorative cooling (K)" = Temperature_LS1982_mean, "Standard deviation surface temperature with evaorative cooling (K)" = Temperature_LS1982_sd, "Workers used" = Workers_used, "Elapsed time" = proc.time()-time_stamp))

  # Write the results to the output file (auto-named from `body` unless out_file is given; pass out_file = NA to skip writing a file), then free up the print out sink.
  if (is.null(out_file)) out_file <- paste0("AccretR_output_", body, ".txt")
  if (!identical(out_file, NA)) {
    sink(out_file)
    print(AccretR_result)
    sink()
  }

  return(AccretR_result)
}
