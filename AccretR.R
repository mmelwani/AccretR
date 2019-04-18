AccretR_hierarchy <- function(){
    
    # Load the package to parallelize
    library(foreach)
    library(iterators)
    library(parallel)
    library(doParallel)
    library(plyr)
    library(ggplot2)
    
    #Set the output file
    
    # Start the timer
    time_stamp <- proc.time()
    
    # Number of workers to be used
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    
    # Invariant values
    Earth_mass_kg <- 5.972e+24
    Earth_radius_m <- 6371000
    G <- 6.67408e-11
    
    # Compositions of the building materials in weight percent. Composition from Lodders and Fegley (1998), normalized to 100 wt. %. Bulk densities from Flynn et al. (2018). Heat capacities from Ostrowski and Bryson (2019).In parentheses: (H wt. %, C wt. %, Mg wt. %, Al wt. %, Si wt. %, S wt. %, Ca wt. %, Fe wt. %, O wt. %, density in kg/m3, heat capacity in J/(kg*K^-1))
    CI_composition <- c(2.07,3.530,9.94,0.89,10.9,5.54,0.95,18.65,47.54,1570,500)
    CM_composition <- c(1.44,2.26,11.8,1.16,13.04,2.77,1.32,21.86,44.34,2270,500)
    CV_composition <- c(0.29,0.55,14.74,1.73,16.18,2.27,1.9,24.22,38.13,2970,500)
    Water_ice <- c(11.19,0,0,0,0,0,0,0,88.808,916.9,1800)
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
        total_body_mass <- 0
        N_particles <- 0
        E_accretion <- 0
        Delta_T <- 0
        total_body_heat_capacity <- 0
        total_body_radius <- 0
        Temperature_diff <- 0
        particle_radius <- 1
        
        repeat{
            # Radius of the accreting particle, in m.
            particle_radius <- runif(1,0.5*particle_radius,2*particle_radius)
            # Volume of the accreting particle, in m^3.
            particle_volume <- 4/3*pi*(particle_radius^3)
            
            # Randomly select material class:
            select_material <- sample(material_list,1)
            unlist_vector <- unlist(select_material)
            # Obtain mass of the accreting particle, in kg, by multiplying particle volume by material class density:
            particle_mass <- (particle_volume)*(as.numeric(unlist_vector[[11]]))
            
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
            
            # Normalize composition to 100 wt. %
            Sum_composition <- sum(H_mass_kg, C_mass_kg, Mg_mass_kg, Al_mass_kg, Si_mass_kg, S_mass_kg, Ca_mass_kg, Fe_mass_kg, O_mass_kg)
            H_wt_perc <- 100*H_mass_kg/Sum_composition
            C_wt_perc <- 100*C_mass_kg/Sum_composition
            Mg_wt_perc <- 100*Mg_mass_kg/Sum_composition
            Al_wt_perc <- 100*Al_mass_kg/Sum_composition
            Si_wt_perc <- 100*Si_mass_kg/Sum_composition
            S_wt_perc <- 100*S_mass_kg/Sum_composition
            Ca_wt_perc <- 100*Ca_mass_kg/Sum_composition
            Fe_wt_perc <- 100*Fe_mass_kg/Sum_composition
            O_wt_perc <- 100*O_mass_kg/Sum_composition
            
            # Heat capacity equation (must come BEFORE mass, radius and growth track equations)
            particle_heat_capacity <- as.numeric(unlist_vector[[12]])
            total_body_heat_capacity <- (total_body_heat_capacity*total_body_mass + particle_heat_capacity*particle_mass)/(total_body_mass + particle_mass)
            
            # Calculate total body mass and radius according to specified growth track (given in Earth radii and masses, from Sotin et al. 2007 Icarus paper, modified to fit Europa)
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
        return(list("H wt. %" = H_wt_perc, "C wt. %" = C_wt_perc, "Mg wt. %" = Mg_wt_perc, "Al wt. %" = Al_wt_perc, "Si wt. %" = Si_wt_perc, "S wt. %" = S_wt_perc, "Ca wt. %" = Ca_wt_perc, "Fe wt. %" = Fe_wt_perc, "O wt. %" = O_wt_perc, "body radius (m)" = total_body_radius, "body mass (kg)" = total_body_mass, "body bulk density (kg/m^3)" = total_body_bulk_density, "Energy of accretion (J/(kg*K^-1))" = E_accretion, "Temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff, "Number of particles"=N_particles))
    }
    
    # Repeat the main subroutine X times and output all results into a transposed dataframe
    Total_bootstrap_run <- foreach(i=1:10000) %dopar% {AccretR_main_subroutine()}
    Total_bootstrap_run_frame <- as.data.frame((Total_bootstrap_run))
    # Repeat the main subroutine X times and output all results into a transposed dataframe
    #Total_bootstrap_run_frame <- as.data.frame(t(replicate(3,{AccretR_main_subroutine()})))
    
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
    # AccretR histogram plots sensu strictu
    N_particles_all_runs_frame <- ldply(N_particles_all_runs, data.frame)
    N_particles_all_runs_plot <<- ggplot(N_particles_all_runs_frame, aes(x=X..i..)) + geom_histogram() + geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Number of accreting particles") + theme_bw()
    
    Total_body_radius_all_runs_frame <- ldply(Total_body_radius_all_runs, data.frame)
    Total_body_radius_all_runs_plot <<- ggplot(Total_body_radius_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Total body radius (m)") + theme_bw()
    
    Total_body_mass_all_runs_frame <- ldply(Total_body_mass_all_runs, data.frame)
    Total_body_mass_all_runs_plot <<- ggplot(Total_body_mass_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Total body mass (kg)") + theme_bw()
    
    Total_body_density_all_runs_frame <- ldply(Total_body_density_all_runs, data.frame)
    Total_body_density_all_runs_plot <<- ggplot(Total_body_density_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab((expression(paste("Bulk body densities ",kg/m^3)))) + theme_bw()
    
    Accretion_energy_all_runs_frame <- ldply(Accretion_energy_all_runs, data.frame)
    Accretion_energy_all_runs_plot <<- ggplot(Accretion_energy_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Accretion energies (J)") + theme_bw()
    
    Temperature_diff_all_runs_frame <- ldply(Temperature_diff_all_runs, data.frame)
    Temperature_diff_all_runs_plot <<- ggplot(Temperature_diff_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab("Temperature difference between body surface and accretion disk (K)") + theme_bw()
    
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
    
    Max_water_wt_perc_all_runs_frame <- ldply(Max_water_wt_perc_all_runs, data.frame)
    Max_water_wt_perc_all_runs_plot <<- ggplot(Max_water_wt_perc_all_runs_frame, aes(x=X..i..)) + geom_histogram() +  geom_vline(aes(xintercept=mean(X..i..)), color="red", linetype="dashed") + xlab((expression(paste(H[2],"O wt. %")))) + theme_bw()
    
    AccretR_composition_result_plot <<- multiplot(H_wt_perc_all_runs_plot, C_wt_perc_all_runs_plot, Mg_wt_perc_all_runs_plot, Al_wt_perc_all_runs_plot, Si_wt_perc_all_runs_plot, S_wt_perc_all_runs_plot, Ca_wt_perc_all_runs_plot, Fe_wt_perc_all_runs_plot, O_wt_perc_all_runs_plot, Max_water_wt_perc_all_runs_plot, cols=2)
    
    # Free up the cluster
    stopCluster(cl)
    
    # Return results
    AccretR_result <<- (list("Mean H wt. %" = H_wt_perc_mean, "Standard deviation H wt. %" = H_wt_perc_sd, "Mean C wt. %" = C_wt_perc_mean, "Standard deviation C wt. %" = C_wt_perc_sd, "Mean Mg wt. %" = Mg_wt_perc_mean, "Standard deviation Mg wt. %" = Mg_wt_perc_sd, "Mean Al wt. %" = Al_wt_perc_mean, "Standard deviation Al wt. %" = Al_wt_perc_sd, "Mean Si wt. %" = Si_wt_perc_mean, "Standard deviation Si wt. %" = Si_wt_perc_sd, "Mean S wt. %" = S_wt_perc_mean, "Standard deviation S wt. %" = S_wt_perc_sd, "Mean Ca wt. %" = Ca_wt_perc_mean, "Standard deviation Ca wt. %" = Ca_wt_perc_sd, "Mean Fe wt. %" = Fe_wt_perc_mean, "Standard deviation Fe wt. %" = Fe_wt_perc_sd, "Mean O wt. %" = O_wt_perc_mean, "Standard deviation O wt. %" = O_wt_perc_sd, "Mean maximum H2O wt. %" = Max_water_wt_perc_mean, "Standard deviation maximum H2O wt. %" = Max_water_wt_perc_sd, "Mean number of particles" = N_particles_mean, "Standard deviation number of particles" = N_particles_sd, "Mean total body radius (m)" = Total_body_radius_mean, "Standard deviation total body radius (m)" = Total_body_radius_sd, "Mean total body mass (kg)" = Total_body_mass_mean, "Standard deviation total body mass (kg)" = Total_body_mass_sd, "Mean body bulk density" = Total_body_density_mean, "Standard deviation body bulk density" = Total_body_density_sd, "Mean accretion energy (J)" = Accretion_energy_mean, "Standard deviation accretion energy (J)" = Accretion_energy_sd, "Mean temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_mean, "Standard deviation temperature difference at body surface over accretion disk temperature (K)" = Temperature_diff_sd, "Workers used" = Workers_used, "Elapsed time" = proc.time()-time_stamp))
    return(AccretR_result)
}

