# Server PBPK model in R
# Nicola Melillo, Hitesh Mistry, 16/06/2021

#setwd("")

# libraries
library(shiny)
library(shinyjs)
library(readxl)
library(RxODE)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# my functions
source("./functions/import_param.R")
source("./functions/PBPK_model_rxode.R")
source("./functions/functions_plot3.R")



shinyServer(function(input, output, session) {
  
  
  values <- reactiveValues(system.out.list = list(),
                           list.sim = list(),
                           names.PBPK = c(),
                           names.ACAT = c(),
                           names.sim = c(),
                           count = 1)  
  flag.clear <- reactiveValues(flag=0)
  flag.simulated <- reactiveValues(flag=0)
  param.drug.library <- reactiveValues(param.drug = list(), 
                                       PK.data = list(), 
                                       PK.data.sel = list())
  
  ### run the model ------------------------------------
  #runmodel <- eventReactive(input$runButton, {
  observeEvent(input$runButton, {
    withProgress(message = 'Simulating the model', value = 0, {
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
      Peff_caco2  <- input$Peff_caco2 * 3600 # [10^-4 cm/s] -> [10^-4 cm/h]
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
        if(input$Sex=="female"){
          filename <- "./data/PBPK_parameters/2021_02_27_pbpk_parameters_human_female.xlsx"
        }else{
          filename <- "./data/PBPK_parameters/2021_02_27_pbpk_parameters_human_male.xlsx"
        }
      }else if(specie=="mouse"){
        filename <- "./data/PBPK_parameters/2021_03_23_pbpk_parameters_mices.xlsx"
      }else if(specie=="beagle"){
        filename <- "./data/PBPK_parameters/2021_04_16_pbpk_parameters_beagles.xlsx"
      }else if(specie=="dog"){
        filename <- "./data/PBPK_parameters/2021_04_16_pbpk_parameters_dogs.xlsx"
      }
      
      param.PBPK       <- getPBPKParam(filename, param_drug, type_part_coeff, specie)
      param.PBPK.rxode <- reorganizeParam.rxode(param.PBPK)
      comp.names       <- c(param.PBPK$comp_PBPK_names, param.PBPK$comp_ACAT_names)
      names.PBPK       <- param.PBPK$comp_PBPK_names
      names.ACAT       <- param.PBPK$comp_ACAT_names
      lo               <- length(comp.names)
      
      param.rxode <- c(param.PBPK.rxode, param_drug)
      
      
      # define schedule
      n_rep <- input$daily_admin * input$days
      time_interval <- 24 / input$daily_admin
      time_interval_end <- 24
      
      
      # select compartment for initial conditions
      if (input$Route==1 || input$Route == 4){
        comp_dose = "venous_blood"
      }else if(input$Route==2){
        comp_dose <- "stomach_s"
      }else if(input$Route==3){
        comp_dose <- "stomach_d"
      }
      
      if(input$Route == 4){
        ev <- eventTable(amount.units="mg", time.units="hr") %>%
          add.dosing(dose=input$dose, dosing.to=comp_dose, nbr.doses=n_rep, dosing.interval=time_interval, dur=input$inf_dur) %>%
          add.sampling(seq(0,n_rep * time_interval + time_interval_end,by=0.01))
      }else{
        ev <- eventTable(amount.units="mg", time.units="hr") %>%
          add.dosing(dose=input$dose, dosing.to=comp_dose, nbr.doses=n_rep, dosing.interval=time_interval) %>%
          add.sampling(seq(0,n_rep * time_interval + time_interval_end,by=0.01))
      }
      
      inits <- c()
      sim1 <- list(param.rxode = param.rxode, ev = ev, inits = inits)
      list.sim <- list(sim1)
      l.param.set <- length(list.sim)
      
      #incProgress(1/2)
      
      ### simulate the model & plot the system 
      system.out.list <- list()
      system.out.list[[1]] <- PBPK.ACAT %>% rxSolve(sim1$param.rxode, sim1$ev, sim1$inits)
      
      # plot
      
      if(input$keepPlots && flag.clear$flag==0){
        names.sim <- paste("sim",values$count,sep="")
        values$system.out.list <- c(values$system.out.list, system.out.list)
        values$list.sim        <- c(values$list.sim, list.sim)
        values$names.PBPK      <- names.PBPK
        values$names.ACAT      <- names.ACAT
        values$names.sim       <- c(values$names.sim, names.sim)
        values$count           <- values$count + 1
      }else{
        names.sim <- paste("sim","1",sep="")
        values$system.out.list <- c(system.out.list)
        values$list.sim        <- c(list.sim)
        values$names.PBPK      <- names.PBPK
        values$names.ACAT      <- names.ACAT
        values$names.sim       <- c(names.sim)
        values$count           <- 2
        flag.clear$flag <- 0
      }
      
      list.out <- list(system.out.list, list.sim, names.sim, names.PBPK, names.ACAT)
      
      flag.simulated$flag = 1
      incProgress(1/2)
      

    })
    
  })
  
  
  ### plot --------------------------------------------------------------------------
  output$PK<-renderPlot({
    withProgress(message = 'Plotting plasma PK', value = 0, {
      if(input$logscale){
        logscale <- 1
      }else{
        logscale <- 0
      }
      
      if(flag.clear$flag==0 && flag.simulated$flag == 1){
        p.tot <- plotPBPK(values$system.out.list, values$list.sim, names.sim = values$names.sim, values$names.PBPK, values$names.ACAT,logscale,param.drug.library$PK.data.sel)
        incProgress(1/4)
        ptlist <- list(p.tot$p.plasma, p.tot$p.f.excr, p.tot$p.f.abs)
        do.call("grid.arrange", c(ptlist, ncol=1, nrow=3))
        flag.simulated$flag = 1
        incProgress(3/4)
      }
    })
  })
  
  output$PK_comp_PBPK<-renderPlot({
    
    if(flag.clear$flag==0 && flag.simulated$flag==1 && input$plotOrgansPK){
      withProgress(message = 'Plotting organs PK', value = 0, {
        if(input$logscale){
          logscale <- 1
        }else{
          logscale <- 0
        }
        p.tot <- plotPBPK(values$system.out.list, values$list.sim, names.sim = values$names.sim, values$names.PBPK, values$names.ACAT,logscale,param.drug.library$PK.data.sel)
        incProgress(1/4)
        do.call("grid.arrange", c(c(p.tot$p.pbpk,p.tot$p.acat.1), ncol=3, nrow=13))
        incProgress(3/4)
      })
    }
    
  })
  
  
  ### functions for activating/deactivating options according to given events -------------------------------------------
  observeEvent(input$Route, {
    if(input$Route == 4){
      shinyjs::enable("inf_dur")
    }else{
      shinyjs::disable("inf_dur")
    }
  })
  
  observeEvent(input$Species, {
    if(input$Species == "human"){
      shinyjs::enable("Sex")
    }else{
      shinyjs::disable("Sex")
    }
  })
  
  
  observeEvent(input$clearPlot, {
    #hide("PK")
    flag.clear$flag <- 1
    flag.simulated$flag <- 0
    param.drug.library$PK.data.sel <- list()
  })
  
  
  ### upload parameters values & PK from library ------------------------------------------------------
  output$libraryDrugs <- renderUI({
    param_drug <- read_excel("./data/library_drugs/drugsParam.xlsx", sheet="drug_param")
    param.drug.library$param.drug <- param_drug
    mydata <- param_drug$drug
    selectInput('selectedDrug', 'Select drug', c(Choose='', mydata), selectize=FALSE)
  })
  
  # upload the parameters from library
  observeEvent(input$uploadDrugParam, {
    
    # remove plots & old stored data for other drugs
    param.drug.library$PK.data.sel <- list()
    flag.clear$flag <- 1
    flag.simulated$flag <- 0
    
    if(input$selectedDrug!=""){
      
      # select drug parameters
      library.param <- param.drug.library$param.drug[param.drug.library$param.drug$drug==input$selectedDrug,]
      
      # molecular related parameters
      updateNumericInput(session,"MW", value = library.param$mw)
      if(library.param$type=="neutral"){
        updateRadioButtons(session,"type",selected=0)
        Csint <- library.param$Cs_w
      }else if(library.param$type=="acid"){
        updateRadioButtons(session,"type",selected=1)
        updateNumericInput(session,"pKa", value = library.param$pKa)
        Csint <- library.param$Cs_w / (1 + 10^(-library.param$pKa + library.param$pH_ref))
      }else{
        updateRadioButtons(session,"type",selected=2)
        updateNumericInput(session,"pKa", value = library.param$pKa)
        Csint <- library.param$Cs_w / (1 + 10^( library.param$pKa - library.param$pH_ref))
      }
      updateNumericInput(session,"logPow", value = library.param$logPow)
      updateNumericInput(session,"fup", value = library.param$fup)
      updateNumericInput(session,"BP", value = library.param$BP)
      
      # dissolution and absorption parameters
      updateNumericInput(session,"r", value = library.param$r)
      updateNumericInput(session,"rho", value = library.param$rho)
      updateNumericInput(session,"Csint", value = Csint)
      updateNumericInput(session,"Peff_caco2", value = library.param$Peff)
      
      # clearance parameters
      updateNumericInput(session,"Clh", value = library.param$Clh)
      updateNumericInput(session,"Clr", value = library.param$Clr)
      updateNumericInput(session,"Clent", value = library.param$Clent)
    }
    
  })
  
  # define selectInput with drug specific data
  output$PKData <- renderUI({
    req(input$selectedDrug)
    PK.data <- read_excel("./data/library_drugs/dataPK.xlsx", sheet=input$selectedDrug)
    param.drug.library$PK.data <- PK.data
    mydata <- unique(PK.data$schedule_ID)
    selectInput('selectPKData', 'Select PK data', c(Choose='', mydata), selectize=FALSE)
  })
  
  # show action button only if a drug is chosen
  output$UploadPKData <- renderUI({
    req(input$selectedDrug)
    actionButton("uploadDrugPK", "Upload PK", width='100px')
  })
  
  # upload the chosen PK data
  observeEvent(input$uploadDrugPK,{
    param.drug.library$PK.data.sel <- param.drug.library$PK.data[param.drug.library$PK.data$schedule_ID==input$selectPKData,]
    
    # update schedule
    updateNumericInput(session,"dose", value = param.drug.library$PK.data.sel$dose[1])
    updateRadioButtons(session,"Species",selected="human")
    updateRadioButtons(session,"Sex",selected="female")
    if(param.drug.library$PK.data.sel$route[1]=="PO"){
      updateRadioButtons(session,"Route",selected=2)
    }else if(param.drug.library$PK.data.sel$route[1]=="IV"){
      updateRadioButtons(session,"Route",selected=4)
      updateNumericInput(session,"inf_dur", value = 0.1) # default to 0.1 to allow better visual comparison
    }
  })
  
  # reset button
  observeEvent(input$resetButton,{
    
    # remove plots & old stored data for other drugs
    param.drug.library$PK.data <- list()
    param.drug.library$PK.data.sel <- list()
    flag.clear$flag <- 1
    flag.simulated$flag <- 0
    updateSelectInput(session, "selectedDrug",selected="")
    
    ### default parameters
    
    # schedule
    updateRadioButtons(session,"Route",selected=2) 
    updateNumericInput(session,"inf_dur", value = 0.5)
    updateNumericInput(session,"dose", value = 10)
    updateNumericInput(session,"days", value = 1)
    updateSliderInput(session,"daily_admin", value=1)
    
    # molecular related parameters
    updateNumericInput(session,"MW", value = 500)
    updateRadioButtons(session,"type",selected=0)
    updateNumericInput(session,"logPow", value = 2)
    updateNumericInput(session,"fup", value = 0.8)
    updateNumericInput(session,"BP", value = 0.8)
    
    # dissolution and absorption parameters
    updateNumericInput(session,"r", value = 25)
    updateNumericInput(session,"rho", value = 1000)
    updateNumericInput(session,"Csint", value = 100)
    updateNumericInput(session,"Peff_caco2", value = 2)
    
    # clearance parameters
    updateNumericInput(session,"Clh", value = 10)
    updateNumericInput(session,"Clr", value = 10)
    updateNumericInput(session,"Clent", value = 0)
    updateRadioButtons(session,"PCM",selected="PT")
    
    # specie and sex
    updateRadioButtons(session,"Species",selected="human")
    updateRadioButtons(session,"Sex",selected="female")
    
  })
  
  
})





