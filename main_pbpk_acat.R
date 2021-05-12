### MAIN - simulated PBPK-ACAT coupled model

### units
#   mass    [mg]
#   volumes [L]
#   time    [h]

### notes
# - the structure is hard coded, therefore, it is not straightforward changing compartments and compartment order...


setwd("C:/Users/nicol/Documents/POSTDOC/projects/systemsforcasting/PBPK/codes/2021_04_19_PBPK_model_r_v03")

### load libraries & functions -------------------------------------------------------------------------------------

# libraries
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


### define drug parameters -----------------------------------------------------------------------------------------
# notes:
#     - diffusion coefficient (Diff) and effective thickness of the hydrodynamic diffusion layer are caculated in getPBPKParam function
#     - types: 0 neutral | 1 acid | 2 base

### parameters related to the molecule
# mw     [g/mol]   molecular weight
# type   flag      0 neutral | 1 acid | 2 base
# Pow    [adim]    octanole to water partition coefficient
# Dvow   [adim]    olive oil to water partition coefficient      - generally this parameter is derived from the logP
# Dvow_s [adim]    corrected oil to water partition coefficient  - generally it is derived from logDov
# fup    [adim]    fraction unbound in plasma
# fut    [adim]    fraction unbound in tissues                   - used to calculate the partition coefficients, derived from fup
# BP     [adim]    blood to plasma ratio
# pKa    [adim]
mw     <- 300
type   <- 2
logPow <- 1
fup    <- 1
BP     <- 0.5
pKa    <- 7

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
r     <- 2.5
rho   <- 1000
Csint <- 100

# if you have the water solubility and the pH of the solvent, you can derive the intrinsic solubility
Csw <- 100
pHw <- 6
if(type==1){
  Csint_w <- Csw / (1 + 10^(pHw - pKa))
}else if(type==2){
  Csint_w <- Csw / (1 + 10^(-pHw + pKa))
}else{
  Csint_w <- logDvow
}
# Csint <- Csint_w  # uncomment if you have water solubility!



# FOR humans: from Papp to Peff, regression!! must be in 10^-4 cm/s (*2)
# FOR mice:   consider to use directly the caco2 permeability... must be in 10^-4 cm/s 
Peff_caco2  <- 1000  # [10^-4 cm/s]
logPeff     <- 0.4926 * log10(Peff_caco2) - 0.1454
Peff        <- 10^(logPeff) * 10^-4 * 3600

### clearances
# CLh    [L/h]     intrinsic hepatic clearance
# CLr    [L/h]     intrinsic renal clearance
# CLent  [L/h]     enterocyte clerance
CLh   <- 10
CLr   <- 10
CLent <- 0

### choose partition coefficients
# "PT"    - Poulin & Theil
# "bere"  - Berezhkovsky
type_part_coeff <- "PT"

### build parameters vector
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
specie          <- "beagle"
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


### define schedule -----------------------------------------------------------------------------------------

func <- PBPK.ACAT

# define schedule
comp_dose         <- "stomach_s"
dose              <- 10 * param.PBPK$general_p["weight"] # [mg/kg] -> [mg]
time_begin        <- 0
time_interval     <- 3
n_rep             <- 3
time_interval_end <- 24

schedule <- defineSchedule(dose, time_begin, time_interval, n_rep, time_interval_end)
delta.t <- 0.01


### define changes for multiple simulations ------------------------------------------------------------------------

# if you want to change some param_drug you need to re-run function for param.PBPK! Especially for formulation properties and parameters used to derive the partition coefficients
# param.PBPK  <- getPBPKParam(filename, param_drug, type_part_coeff)
sim1 <- list(func=func, schedule=schedule, delta.t=delta.t, comp_dose=comp_dose, param.PBPK=param.PBPK)


sim2 <- sim1
sim2$param.PBPK$param_drug["Csint"] <- sim1$param.PBPK$param_drug["Csint"]/1000

sim3 <- sim1
sim3$param.PBPK$param_drug["Peff"] <- sim1$param.PBPK$param_drug["Peff"]*10

sim4 <- sim1
sim4$param.PBPK$param_drug["Peff"] <- sim1$param.PBPK$param_drug["Peff"]*0

sim5 <- sim1
sim5$param.PBPK$param_drug["CLh"] <- sim1$param.PBPK$param_drug["CLh"]*0

sim6 <- sim1
sim6$param.PBPK$param_drug["CLr"] <- sim1$param.PBPK$param_drug["CLr"]*0

sim7 <- sim1
sim7$param.PBPK$param_drug["CLh"] <- sim1$param.PBPK$param_drug["CLh"]*0
sim7$param.PBPK$param_drug["CLr"] <- sim1$param.PBPK$param_drug["CLr"]*0

# here you can add any simulation param that you want...
list.sim <- list(sim1, sim2, sim3)
l.param.set <- length(list.sim)

names.sim <- c("baseline", "Csint/100", "Peff*10", "sim4", "sim5", "sim6", "sim7")

### simulate the model & plot the system ---------------------------------------------------------------------------
system.out.list <- list()
for(i in 1:l.param.set){
  print(i)
  simi <- list.sim[[i]]
  system.out.list[[i]] <- simulation(simi$func, 
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
p.tot$p.f.excr    # plot fraction excreted
p.tot$p.f.abs     # plot fraction absorbed
p.tot$p.plasma    # plot plasma concentration PK

### some references ------------------------------------------------------------------------------------------------
# (*1) Poulin & Theil 2001 https://doi.org/10.1002/jps.10005
# (*2) Sun 2002 https://doi.org/10.1023/a:1020483911355




