	# AccretR, a program to calculate the composition of an icy ocean world using the building blocks of our solar system.
	# Copyright (C) 2019 Mohit Melwani Daswani

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

AccretR <- function(){
    
    # Load the package to parallelize
    library(foreach)
    library(iterators)
    library(parallel)
    library(doParallel)
    library(plyr)
    library(ggplot2)
	library(scales)
    
    #Set the output file
    
    # Start the timer
    time_stamp <- proc.time()
    
    # Input number of workers/cores to be used. Function detectCores gives the total number of cores available on your system. Change if you will need fewer cores (a good idea if you are running other tasks, and AccretR is using too much memory). 
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    
    # Invariant values
    Earth_mass_kg <- 5.972e+24
    Earth_radius_m <- 6371000
    G <- 6.67408e-11
    
    # Compositions of the building materials in weight percent. Compositions from Lodders and Fegley (1998) for CM and CV chondrites, Lodders (2010) for CI chondrites, CLay et al. (2017) for chlorine in all meteorites; all normalized to 100 wt. %. Bulk densities from Flynn et al. (2018). Heat capacities from Ostrowski and Bryson (2019).In parentheses: (H wt. %, C wt. %, Mg wt. %, Al wt. %, Si wt. %, S wt. %, Ca wt. %, Fe wt. %, O wt. %, Na wt. %, K wt. %, Cl wt. %, density in kg/m3, heat capacity in J/(kg*K^-1)). Comment out the building blocks you want to leave out.
    CI_composition <- c(2.014,3.558,9.794,0.869,10.939,5.469,0.943,18.913,46.924,0.510,0.057,0.012,1570,500)
    CM_composition <- c(1.431,2.248,11.751,1.155,12.977,2.759,1.318,21.764,44.142,0.399,0.038,0.020,2270,500)
    CV_composition <- c(0.287,0.544,14.680,1.725,16.117,2.258,1.889,24.125,37.983,0.349,0.037,0.005,2970,500)
    Water_ice <- c(11.19,0,0,0,0,0,0,0,88.808,0,0,0,916.9,1800)
    material_list <- list(list("CI",CI_composition),list("CM",CM_composition),list("CV",CV_composition), list("Ice",Water_ice))
    #material_list <- list(list("CI",CI_composition),list("CM",CM_composition),list("CV",CV_composition))
    
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
        total_body_mass <- 0
        N_particles <- 0
        E_accretion <- 0
        Delta_T <- 0
        total_body_heat_capacity <- 0
        total_body_radius <- 10
        Temperature_diff <- 0
        
        repeat{
            # Radius of the accreting particles, in m. Here, it is a random uniform distribution proportional to the total body radius (a type of runaway growth).
            particle_radius <- runif(1,0.1*total_body_radius,0.5*total_body_radius)
            # Volume of the accreting particle, in m^3.
            particle_volume <- 4/3*pi*(particle_radius^3)
            
            # Randomly select material class:
            select_material <- sample(material_list,1)
            unlist_vector <- unlist(select_material)
            # Obtain mass of the accreting particle, in kg, by multiplying particle volume by material class density:
            particle_mass <- (particle_volume)*(as.numeric(unlist_vector[[14]]))
            
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
            
            # Normalize composition to 100 wt. %
            Sum_composition <- sum(H_mass_kg, C_mass_kg, Mg_mass_kg, Al_mass_kg, Si_mass_kg, S_mass_kg, Ca_mass_kg, Fe_mass_kg, O_mass_kg, Na_mass_kg, K_mass_kg, Cl_mass_kg)
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
            
            # Heat capacity equation (must come BEFORE mass, radius and growth track equations)
            particle_heat_capacity <- as.numeric(unlist_vector[[15]])
            total_body_heat_capacity <- (total_body_heat_capacity*total_body_mass + particle_heat_capacity*particle_mass)/(total_body_mass + particle_mass)
            
            # Calculate total body mass and radius according to specified growth track (given in Earth radii and masses, from Sotin et al. 2007 Icarus paper, here modified to fit Europa specifically)
            total_body_mass <- total_body_mass + particle_mass
            total_body_radius <- ((1.072*(total_body_mass/Earth_mass_kg)^0.306)*Earth_radius_m)
            total_body_bulk_density <- total_body_mass/(4/3*pi*(total_body_radius^3))
            
            # Energy and temperature equations (must come AFTER mass, radius and growth track equations). E_accretion and Temperature_diff is the energy of accretion, and temperature difference between the body's surface and the accretion disk, modified from Stevenson et al. (1986) Origins of Satellites chapter in Satellites book.
            E_accretion <- E_accretion + G*total_body_mass/total_body_radius
            Temperature_diff <- Temperature_diff + (G*total_body_mass)/(total_body_radius*total_body_heat_capacity)
            
            # Total number of accreted particles
            N_particles=N_particles+1
            
            # Exit if radius = specified radius in meters
            if (total_body_radius>=1560800) break
        }
        return(list("H wt. %" = H_wt_perc, "C wt. %" = C_wt_perc, "Mg wt. %" = Mg_wt_perc, "Al wt. %" = Al_wt_perc, "Si wt. %" = Si_wt_perc, "S wt. %" = S_wt_perc, "Ca wt. %" = Ca_wt_perc, "Fe wt. %" = Fe_wt_perc, "O wt. %" = O_wt_perc, "Na wt. %" = Na_wt_perc, "K wt. %" = K_wt_perc, "Cl wt. %" = Cl_wt_perc, "body radius (m)" = total_body_radius, "body mass (kg)" = total_body_mass, "body bulk density (kg/m^3)" = total_body_bulk_density, "Energy of accretion (J/(kg*K^-1))" = E_accretion, "Temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff, "Number of particles"=N_particles))
    }
    
    # Number of full bodies to build for statistics. Repeat the main subroutine X times and output all results into a transposed dataframe.
    Total_bootstrap_run <- foreach(i=1:10000) %dopar% {AccretR_main_subroutine()}
    Total_bootstrap_run_frame <- as.data.frame((Total_bootstrap_run))
    
    # Output means and standard deviations of the data gathered in the data frame
    N_particles_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Number.of.particles", colnames(Total_bootstrap_run_frame))]))
    N_particles_mean <- mean(N_particles_all_runs)
    N_particles_sd <- sd(N_particles_all_runs)
    H_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("H.wt.", colnames(Total_bootstrap_run_frame))]))
    H_wt_perc_mean <- mean(H_wt_perc_all_runs)
    H_wt_perc_sd <- sd(H_wt_perc_all_runs) 	
    C_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("C.wt.", colnames(Total_bootstrap_run_frame))]))
    C_wt_perc_mean <- mean(C_wt_perc_all_runs)
    C_wt_perc_sd <- sd(C_wt_perc_all_runs) 	
    Mg_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Mg.wt.", colnames(Total_bootstrap_run_frame))]))
    Mg_wt_perc_mean <- mean(Mg_wt_perc_all_runs)
    Mg_wt_perc_sd <- sd(Mg_wt_perc_all_runs) 	
    Al_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Al.wt.", colnames(Total_bootstrap_run_frame))]))
    Al_wt_perc_mean <- mean(Al_wt_perc_all_runs)
    Al_wt_perc_sd <- sd(Al_wt_perc_all_runs) 	
    Si_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Si.wt.", colnames(Total_bootstrap_run_frame))]))
    Si_wt_perc_mean <- mean(Si_wt_perc_all_runs)
    Si_wt_perc_sd <- sd(Si_wt_perc_all_runs) 	
    S_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("S.wt.", colnames(Total_bootstrap_run_frame))]))
    S_wt_perc_mean <- mean(S_wt_perc_all_runs)
    S_wt_perc_sd <- sd(S_wt_perc_all_runs)
    Ca_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Ca.wt.", colnames(Total_bootstrap_run_frame))]))
    Ca_wt_perc_mean <- mean(Ca_wt_perc_all_runs)
    Ca_wt_perc_sd <- sd(Ca_wt_perc_all_runs)
    Fe_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Fe.wt.", colnames(Total_bootstrap_run_frame))]))
    Fe_wt_perc_mean <- mean(Fe_wt_perc_all_runs)
    Fe_wt_perc_sd <- sd(Fe_wt_perc_all_runs) 	
    O_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("O.wt.", colnames(Total_bootstrap_run_frame))]))
    O_wt_perc_mean <- mean(O_wt_perc_all_runs)
    O_wt_perc_sd <- sd(O_wt_perc_all_runs)
	Na_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Na.wt.", colnames(Total_bootstrap_run_frame))]))
    Na_wt_perc_mean <- mean(Na_wt_perc_all_runs)
    Na_wt_perc_sd <- sd(Na_wt_perc_all_runs)
	K_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("K.wt.", colnames(Total_bootstrap_run_frame))]))
    K_wt_perc_mean <- mean(K_wt_perc_all_runs)
    K_wt_perc_sd <- sd(K_wt_perc_all_runs)
	Cl_wt_perc_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Cl.wt.", colnames(Total_bootstrap_run_frame))]))
    Cl_wt_perc_mean <- mean(Cl_wt_perc_all_runs)
    Cl_wt_perc_sd <- sd(Cl_wt_perc_all_runs)
    Total_body_radius_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("body.radius", colnames(Total_bootstrap_run_frame))]))
    Total_body_radius_mean <- mean(Total_body_radius_all_runs)
    Total_body_radius_sd <- sd(Total_body_radius_all_runs)
    Total_body_density_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("body.bulk.density", colnames(Total_bootstrap_run_frame))]))
    Total_body_density_mean <- mean(Total_body_density_all_runs)
    Total_body_density_sd <- sd(Total_body_density_all_runs)
    Total_body_mass_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("body.mass", colnames(Total_bootstrap_run_frame))]))
    Total_body_mass_mean <- mean(Total_body_mass_all_runs)
    Total_body_mass_sd <- sd(Total_body_mass_all_runs)
    Accretion_energy_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Energy.of.accretion", colnames(Total_bootstrap_run_frame))]))
    Accretion_energy_mean <- mean(Accretion_energy_all_runs)
    Accretion_energy_sd <- sd(Accretion_energy_all_runs)
    Temperature_diff_all_runs <- as.numeric(unlist(Total_bootstrap_run_frame[, grep("Temperature.difference", colnames(Total_bootstrap_run_frame))]))
    Temperature_diff_mean <- mean(Temperature_diff_all_runs)
    Temperature_diff_sd <- sd(Temperature_diff_all_runs)
    #Assuming all H is in water, calculate all water available
    Max_water_wt_perc_all_runs <- H_wt_perc_all_runs*(1/1.00794)*(1/2)*((1.00794*2+15.999)/1)
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
        }
    }
	
	# Fancy scientific labels in the ggplot graphs. Choose labels = fancy_scientific" for 10^y notation, or "scientific" for +e^y notation.
	fancy_scientific <- function(l) {
        # turn in to character string in scientific notation
        l <- format(l, scientific = TRUE)
        # Print "0 x 10\u{207a}" as 0
        l <- gsub("0e\\+00","0",l)
        # quote the part before the exponent to keep all the digits
        l <- gsub("^(.*)e", "'\\1'e", l)
        # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
        l <- gsub("e\\+","e",l)
        # turn the 'e+' into plotmath format
        l <- gsub("e", "%*%10^", l)
        # convert 1x10^ or 1.000x10^ -> 10^ 
        l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
        # return this as an expression
        parse(text=l)
    }
	
    # AccretR histogram plots sensu strictu
    N_particles_all_runs_frame <- ldply(N_particles_all_runs, data.frame)
    N_particles_all_runs_plot <<- ggplot(N_particles_all_runs_frame, aes(x=X..i..)) + geom_histogram() + geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Number of accreting planetesimals") + theme_bw()
    
    Total_body_radius_all_runs_frame <- ldply(Total_body_radius_all_runs, data.frame)
    Total_body_radius_all_runs_plot <<- ggplot(Total_body_radius_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + scale_x_continuous(labels = scientific) + xlab("Total body radius (m)") + theme_bw()
    
    Total_body_mass_all_runs_frame <- ldply(Total_body_mass_all_runs, data.frame)
    Total_body_mass_all_runs_plot <<- ggplot(Total_body_mass_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + scale_x_continuous(labels = scientific) + xlab("Total body mass (kg)") + theme_bw()
    
    Total_body_density_all_runs_frame <- ldply(Total_body_density_all_runs, data.frame)
    Total_body_density_all_runs_plot <<- ggplot(Total_body_density_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab((expression(paste("Bulk body densities (",kg/m^3,")")))) + theme_bw()
    
    Accretion_energy_all_runs_frame <- ldply(Accretion_energy_all_runs, data.frame)
    Accretion_energy_all_runs_plot <<- ggplot(Accretion_energy_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + scale_x_continuous(labels = scientific) + xlab("Accretion energies (J)") + theme_bw()
    
    Temperature_diff_all_runs_frame <- ldply(Temperature_diff_all_runs, data.frame)
    Temperature_diff_all_runs_plot <<- ggplot(Temperature_diff_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + scale_x_continuous(labels = scientific) + xlab("Temperature difference between body surface and accretion disk (K)") + theme_bw()
    
    AccretR_body_result_plot <<- multiplot(N_particles_all_runs_plot, Total_body_radius_all_runs_plot, Total_body_mass_all_runs_plot, Total_body_density_all_runs_plot, Accretion_energy_all_runs_plot, Temperature_diff_all_runs_plot, cols=2)
    
    H_wt_perc_all_runs_frame <- ldply(H_wt_perc_all_runs, data.frame)
    H_wt_perc_all_runs_plot <<- ggplot(H_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk H wt. %") + theme_bw()
    
    C_wt_perc_all_runs_frame <- ldply(C_wt_perc_all_runs, data.frame)
    C_wt_perc_all_runs_plot <<- ggplot(C_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk C wt. %") + theme_bw()
    
    Mg_wt_perc_all_runs_frame <- ldply(Mg_wt_perc_all_runs, data.frame)
    Mg_wt_perc_all_runs_plot <<- ggplot(Mg_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Mg wt. %") + theme_bw()
    
    Al_wt_perc_all_runs_frame <- ldply(Al_wt_perc_all_runs, data.frame)
    Al_wt_perc_all_runs_plot <<- ggplot(Al_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Al wt. %") + theme_bw()
    
    Si_wt_perc_all_runs_frame <- ldply(Si_wt_perc_all_runs, data.frame)
    Si_wt_perc_all_runs_plot <<- ggplot(Si_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Si wt. %") + theme_bw()
    
    S_wt_perc_all_runs_frame <- ldply(S_wt_perc_all_runs, data.frame)
    S_wt_perc_all_runs_plot <<- ggplot(S_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk S wt. %") + theme_bw()
    
    Ca_wt_perc_all_runs_frame <- ldply(Ca_wt_perc_all_runs, data.frame)
    Ca_wt_perc_all_runs_plot <<- ggplot(Ca_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Ca wt. %") + theme_bw()
    
    Fe_wt_perc_all_runs_frame <- ldply(Fe_wt_perc_all_runs, data.frame)
    Fe_wt_perc_all_runs_plot <<- ggplot(Fe_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Fe wt. %") + theme_bw()
    
    O_wt_perc_all_runs_frame <- ldply(O_wt_perc_all_runs, data.frame)
    O_wt_perc_all_runs_plot <<- ggplot(O_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk O wt. %") + theme_bw()
	
	Na_wt_perc_all_runs_frame <- ldply(Na_wt_perc_all_runs, data.frame)
    Na_wt_perc_all_runs_plot <<- ggplot(Na_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Na wt. %") + theme_bw()
	
	K_wt_perc_all_runs_frame <- ldply(K_wt_perc_all_runs, data.frame)
    K_wt_perc_all_runs_plot <<- ggplot(K_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk K wt. %") + theme_bw()
	
	Cl_wt_perc_all_runs_frame <- ldply(Cl_wt_perc_all_runs, data.frame)
    Cl_wt_perc_all_runs_plot <<- ggplot(Cl_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Bulk Cl wt. %") + theme_bw()
    
    Max_water_wt_perc_all_runs_frame <- ldply(Max_water_wt_perc_all_runs, data.frame)
    Max_water_wt_perc_all_runs_plot <<- ggplot(Max_water_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab((expression(paste("Maximal ", H[2],"O wt. %")))) + theme_bw()
    
    AccretR_composition_result_plot <<- multiplot(H_wt_perc_all_runs_plot, C_wt_perc_all_runs_plot, Mg_wt_perc_all_runs_plot, Al_wt_perc_all_runs_plot, Si_wt_perc_all_runs_plot, S_wt_perc_all_runs_plot, Ca_wt_perc_all_runs_plot, Fe_wt_perc_all_runs_plot, O_wt_perc_all_runs_plot, Na_wt_perc_all_runs_plot, K_wt_perc_all_runs_plot, Cl_wt_perc_all_runs_plot, Max_water_wt_perc_all_runs_plot, cols=3)
    
    # Free up the cores in the cluster
    stopCluster(cl)
    
    # Return results
    AccretR_result <<- (list("Mean H wt. %" = H_wt_perc_mean, "Standard deviation H wt. %" = H_wt_perc_sd, "Mean C wt. %" = C_wt_perc_mean, "Standard deviation C wt. %" = C_wt_perc_sd, "Mean Mg wt. %" = Mg_wt_perc_mean, "Standard deviation Mg wt. %" = Mg_wt_perc_sd, "Mean Al wt. %" = Al_wt_perc_mean, "Standard deviation Al wt. %" = Al_wt_perc_sd, "Mean Si wt. %" = Si_wt_perc_mean, "Standard deviation Si wt. %" = Si_wt_perc_sd, "Mean S wt. %" = S_wt_perc_mean, "Standard deviation S wt. %" = S_wt_perc_sd, "Mean Ca wt. %" = Ca_wt_perc_mean, "Standard deviation Ca wt. %" = Ca_wt_perc_sd, "Mean Fe wt. %" = Fe_wt_perc_mean, "Standard deviation Fe wt. %" = Fe_wt_perc_sd, "Mean O wt. %" = O_wt_perc_mean, "Standard deviation O wt. %" = O_wt_perc_sd, "Mean Na wt. %" = Na_wt_perc_mean, "Standard deviation Na wt. %" = Na_wt_perc_sd, "Mean K wt. %" = K_wt_perc_mean, "Standard deviation K wt. %" = K_wt_perc_sd, "Mean Cl wt. %" = Cl_wt_perc_mean, "Standard deviation Cl wt. %" = Cl_wt_perc_sd, "Mean maximum H2O wt. %" = Max_water_wt_perc_mean, "Standard deviation maximum H2O wt. %" = Max_water_wt_perc_sd, "Mean number of particles" = N_particles_mean, "Standard deviation number of particles" = N_particles_sd, "Mean total body radius (m)" = Total_body_radius_mean, "Standard deviation total body radius (m)" = Total_body_radius_sd, "Mean total body mass (kg)" = Total_body_mass_mean, "Standard deviation total body mass (kg)" = Total_body_mass_sd, "Mean body bulk density" = Total_body_density_mean, "Standard deviation body bulk density" = Total_body_density_sd, "Mean accretion energy (J)" = Accretion_energy_mean, "Standard deviation accretion energy (J)" = Accretion_energy_sd, "Mean temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_mean, "Standard deviation temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_sd, "Workers used" = Workers_used, "Elapsed time" = proc.time()-time_stamp))
    return(AccretR_result)
}
