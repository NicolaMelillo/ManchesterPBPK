# Server PBPK model in R
# Nicola Melillo, Hitesh Mistry, 11/05/2021

setwd("C:/Users/nicol/Documents/POSTDOC/projects/systemsforcasting/PBPK/codes/2021_05_12_PBPK_app_v05p1")

# libraries
library(shiny)
library(readxl)
library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# my functions
source("./functions/import_param.R")
source("./functions/PBPK_model.R")
source("./functions/functions_simulations.R")
source("./functions/functions_plot2.R")



shinyServer(function(input, output, session) {
  
  global.system.out <- reactive(list())
  
  runmodel <- eventReactive(input$runButton, {
    
    mw     <- input$MW
    type   <- as.numeric(input$type)
    logPow <- input$logPow
    fup    <- input$fup
    BP     <- input$BP
    pKa    <- input$pKa
    
    # derive the other molecular parameters (*1)
    # Handerson Hasselback equation
    fut     <- 1/( 1 + 0.5*(1-fup)/fup )
    logDvow <- 1.115 * logPow - 1.35
    pH.tiss <- 7.4                         # HP, see (*1)
    if(type==1){
      logDvow_s <- logDvow - log10(1 + 10^(pH.tiss - pKa))
    }else if(type==2){
      logDvow_s <- logDvow - log10(1 + 10^(-pH.tiss + pKa))
    }else{
      logDvow_s <- logDvow
    }
    
    ### formulation related parameters
    # r      [um]      radius of the particle size of the formulation
    # rho    [g/L]     density of the formulation
    # Csint  [mg/L]    intrinsic solubility
    # Peff   [cm/h]    effective permeability across gut wall
    r     <- input$r
    rho   <- input$rho
    Csint <- input$Csint
    
    # if you have the water solubility and the pH of the solvent, you can derive the intrinsic solubility
    #Csw <- 100
    #pHw <- 6
    #if(type==1){
    #  Csint_w <- Csw / (1 + 10^(pHw - pKa))
    #}else if(type==2){
    #  Csint_w <- Csw / (1 + 10^(-pHw + pKa))
    #}else{
    #  Csint_w <- logDvow
    #}
    # Csint <- Csint_w  # uncomment if you have water solubility!
    
    
    
    # FOR humans: from Papp to Peff, regression!! must be in 10^-4 cm/s (*2)
    # FOR mice:   consider to use directly the caco2 permeability... must be in 10^-4 cm/s 
    Peff_caco2  <- input$Peff_caco2  # [10^-4 cm/s]
    #logPeff     <- 0.4926 * log10(Peff_caco2) - 0.1454
    Peff        <- Peff_caco2 # 10^(logPeff) * 10^-4 * 3600
    
    ### clearances
    # CLh    [L/h]     intrinsic hepatic clearance
    # CLr    [L/h]     intrinsic renal clearance
    # CLent  [L/h]     enterocyte clerance
    CLh   <- input$Clh
    CLr   <- input$Clr
    CLent <- input$Clent
    
    ### choose partition coefficients
    # "PT"    - Poulin & Theil
    # "bere"  - Berezhkovsky
    type_part_coeff <- input$PCM
    
    param_drug = c(Pow    = 10^(logPow),
                   Dvow   = 10^(logDvow),
                   Dvow_s = 10^(logDvow_s),
                   fup    = fup,
                   fut    = fut,
                   BP     = BP,
                   CLh    = CLh,
                   CLr    = CLr,
                   CLent  = CLent,
                   Peff   = Peff,
                   r      = r,
                   mw     = mw,
                   rho    = rho,
                   Csint  = Csint,
                   pKa    = pKa,
                   type   = type
    )
    # load PBPK parameters
    specie          <- input$Species
    if(specie=="human"){
      filename <- "./data/PBPK_parameters/2021_02_27_pbpk_parameters_human.xlsx"
    }else if(specie=="mouse"){
      filename <- "./data/PBPK_parameters/2021_03_23_pbpk_parameters_mices.xlsx"
    }else if(specie=="beagle"){
      filename <- "./data/PBPK_parameters/2021_04_16_pbpk_parameters_beagles.xlsx"
    }else if(specie=="dog"){
      filename <- "./data/PBPK_parameters/2021_04_16_pbpk_parameters_dogs.xlsx"
    }
    
    
    param.PBPK      <- getPBPKParam(filename, param_drug, type_part_coeff, specie)
    
    comp.names      <- c(param.PBPK$comp_PBPK_names, param.PBPK$comp_ACAT_names)
    lo              <- length(comp.names)
    func <- PBPK.ACAT
    
    # define schedule
    if (input$Schedule==1){
      n_rep <- 1 * input$days
      time_interval <- 24
      if (input$days==1){
        time_interval_end <- 0
      }else{
        time_interval_end <- 24
      }
    }else if(input$Schedule==2){
      n_rep <- 2 * input$days
      time_interval <- 12
      time_interval_end <- 12
    }
    
    # select compartment for initial conditions
    if (input$Route==1){
      comp_dose = "venous_blood"
    }else if(input$Route==2){
      comp_dose <- "stomach_s"
    }else if(input$Route==3){
      comp_dose <- "stomach_d"
    }
    
    dose              <- input$dose #* param.PBPK$general_p["weight"] # [mg/kg] -> [mg]
    time_begin        <- 0
    #time_interval     <- 1
    #n_rep             <- 1
    
    
    
    schedule <- defineSchedule(dose, time_begin, time_interval, n_rep, time_interval_end)
    delta.t <- 0.01
    
    sim1 <- list(func=func, schedule=schedule, delta.t=delta.t, comp_dose=comp_dose, param.PBPK=param.PBPK)
    list.sim <- list(sim1)
    l.param.set <- length(list.sim)
    
    names.sim <- c("sim1")
    
    ### simulate the model & plot the system ---------------------------------------------------------------------------
    system.out.list <- list()
    for(ii in 1:l.param.set){
      simi <- list.sim[[ii]]
      system.out.list[[ii]] <- simulation(simi$func, 
                                         simi$schedule, 
                                         simi$delta.t, 
                                         simi$comp_dose, 
                                         simi$param.PBPK )
    }

    # plot
    p.tot <- plotPBPK(system.out.list, list.sim, names.sim = names.sim)
    
    # plot all the organs and tissues (except absorbed and sink enterocytes compartments)
    #do.call("grid.arrange", c(p.tot$p.pbpk, ncol=5, nrow=4))      # plot all organs mass PK
    #do.call("grid.arrange", c(p.tot$p.acat.1, ncol=6, nrow=4))    # plot ACAT compartments mass PK
    #p.tot$p.f.excr    # plot fraction excreted
    #p.tot$p.f.abs     # plot fraction absorbed
    #p.tot$p.plasma    # plot plasma concentration PK
    
    ptlist <- list(p.tot$p.plasma, p.tot$p.f.excr, p.tot$p.f.abs)
    do.call("grid.arrange", c(ptlist, ncol=1, nrow=3)) 
    #return(system.out.list)
  })
  
  
  output$PK<-renderPlot({
    
    runmodel()
    #t<-runmodel$time
    #c<-runmodel$venous_blood
    
    #plot(t,c,type="l",xlim=c(0,24),
    #     xlab="Time (hours)",ylab="Plasma Conc. (mg/L)")

    })
  
  
  
  
})