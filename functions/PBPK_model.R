### functions describing the PBPK model

### Units
# x        - [mg]
# volumes  - [L]
# time     - [h]

### PBPK organs for which we have the various parameters (in order)

# PBPK distribution
# 1  lungs
# 2  brain
# 3  heart
# 4  kidneys
# 5  bone
# 6  muscle
# 7  stomach
# 8  spleen
# 9  liver
# 10 gut
# 11 pancreas
# 12 skin
# 13 fat
# 14 arterial_blood
# 15 venous_blood

# PBPK sink
# 16 sink liver
# 17 sink kidney

# ACAT solid
# 18 stomach_s
# 19 duodenum_s
# 20 jejunum1_s
# 21 jejunum2_s
# 22 ileum1_s
# 23 ileum2_s
# 24 ileum3_s

# ACAT dissolved
# 25 stomach_d
# 26 duodenum_d
# 27 jejunum1_d
# 28 jejunum2_d
# 29 ileum1_d
# 30 ileum2_d
# 31 ileum3_d

# ACAT enterocytes
# 32 duodenum_e
# 33 jejunum1_e
# 34 jejunum2_e
# 35 ileum1_e
# 36 ileum2_e
# 37 ileum3_e

# ACAT_absorbed
# 38 duodenum_ea
# 39 jejunum1_ea
# 40 jejunum2_ea
# 41 ileum1_ea
# 42 ileum2_ea
# 43 ileum3_ea

# ACAT cleared
# 44 duodenum_ea
# 45 jejunum1_ea
# 46 jejunum2_ea
# 47 ileum1_ea
# 48 ileum2_ea
# 49 ileum3_ea

# ACAT sink
# 50 sink_s
# 51 sink_d


PBPK.ACAT <- function(t, x, parms) {
  
  # define some support variable (indices hard coded)
  x_pbpk      <- x[1:17]
  x_acat      <- x[18:51]
  idx_liver   <- 9
  idx_ent_abs <- 21:26

  # derive the differential equations
  dx_pbpk   <- PBPK_distribution(t, x_pbpk, parms)
  dx_acat   <- ACAT_model(t, x_acat, parms)
  dx_pbpk_v <- dx_pbpk[[1]]
  dx_acat_v <- dx_acat[[1]]
  
  # couple the models
  dx_pbpk_v[idx_liver] <- dx_pbpk_v[idx_liver] + sum(dx_acat_v[idx_ent_abs])
  dx <- list(c(dx_pbpk_v, dx_acat_v))
  
  # output
  return(dx)
  
}

### PBPK distribution function --------------------------------------------------------------------
PBPK_distribution <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    # some design parameters
    lo <- length(comp_PBPK_names)
    map_c <- structure(1:lo, names=comp_PBPK_names)
    
    # define type of tissues and useful indices
    parallel_tissues  <- c("brain","heart","kidneys","bone","muscle","skin","fat")
    splanchnic_organs <- c("stomach", "spleen", "pancreas", "gut")
    lpt       <- length(parallel_tissues)
    lso       <- length(splanchnic_organs)
    idx_art   <- map_c["arterial_blood"]
    idx_ven   <- map_c["venous_blood"]
    idx_lung  <- map_c["lungs"]
    idx_liv   <- map_c["liver"]
    idx_sclh  <- map_c["sink_CLh"]
    idx_sclr  <- map_c["sink_CLr"]
    idx_kid   <- map_c["kidneys"]
    
    # preallocation of some vectors and support variable...
    dx  <- numeric(lo)
    bf_tot      <- 0
    outflow_tot <- 0
    outflow_spl <- 0
    
    # loop for defining parallel tissues equations...
    for(i in 1:lpt){
      idx_i         <- map_c[parallel_tissues[i]]
      org_i_outflow <- organ_bf[idx_i] * (x[idx_i]/organ_v[idx_i])/(part_coeff[idx_i]/param_drug["BP"])
      dx[idx_i]     <- organ_bf[idx_i] * x[idx_art]/organ_v[idx_art] - org_i_outflow
      
      bf_tot        <- bf_tot + organ_bf[idx_i]
      outflow_tot   <- outflow_tot + org_i_outflow
    }
    
    # add elimination to the kidney compartment
    clear_kid   <- param_drug["CLr"] * x[idx_kid]/organ_v[idx_kid]
    dx[idx_kid] <- dx[idx_kid] - clear_kid
    
    # loop for defining splanchnic organs equations
    for(i in 1:lso){
      idx_i         <- map_c[splanchnic_organs[i]]
      org_i_outflow <- organ_bf[idx_i] * (x[idx_i]/organ_v[idx_i])/(part_coeff[idx_i]/param_drug["BP"])
      dx[idx_i]     <- organ_bf[idx_i] * x[idx_art]/organ_v[idx_art] - org_i_outflow
      
      bf_tot        <- bf_tot + organ_bf[idx_i]
      outflow_spl   <- outflow_spl + org_i_outflow
    }
    
    # liver equations
    outflow_liv  <- organ_bf[idx_liv] * (x[idx_liv]/organ_v[idx_liv])/(part_coeff[idx_liv]/param_drug["BP"])
    clear_liv    <- param_drug["CLh"] * x[idx_liv]/organ_v[idx_liv]
    dx[idx_liv]  <- outflow_spl + organ_bf[idx_liv] * x[idx_art]/organ_v[idx_art] - outflow_liv - clear_liv
    bf_tot       <- bf_tot + organ_bf[idx_liv]
    outflow_tot  <- outflow_tot + outflow_liv
    
    # venous blood equation
    outflow_ven <- bf_tot * x[idx_ven]/organ_v[idx_ven]
    dx[idx_ven] <- outflow_tot - outflow_ven
    
    # lung equation
    outflow_lung <- bf_tot * ( x[idx_lung]/organ_v[idx_lung] ) / ( part_coeff[idx_lung]/param_drug["BP"] )
    dx[idx_lung] <- outflow_ven - outflow_lung
    
    # arterial blood equation
    dx[idx_art] <- outflow_lung - bf_tot * x[idx_art]/organ_v[idx_art]
    
    # sink compartments
    dx[idx_sclh] <- clear_liv # liver sink
    dx[idx_sclr] <- clear_kid # kidney sink
    
    res  <- structure(dx, names=comp_PBPK_names)
    list(res)
  })
}


### ACAT model function ---------------------------------------------------------------------------
ACAT_model <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    ### define some useful parameters
    
    # map idx name
    lo <- length(comp_ACAT_names)
    map_c <- structure(1:lo, names=comp_ACAT_names)
    section_names <- names(ACAT_l)
    
    # preallocate dx vector
    dx  <- numeric(lo)
    
    # compartment names
    comp_solid <- comp_ACAT_names[1:7]    # solid drug in lumen
    comp_diss  <- comp_ACAT_names[8:14]   # dissolved drug in lume
    comp_ent   <- comp_ACAT_names[15:20]  # drug in enterocytes
    comp_abs   <- comp_ACAT_names[21:26]  # drug absorbed in a given section
    comp_cl    <- comp_ACAT_names[27:32]  # drug cleared in a given section
    
    # term in front of the dissolution equation
    Kd1 <- 10^5*(3*param_drug["Diff"]/(param_drug["rho"]*param_drug["ht"]*param_drug["r"])); # [L/(h*mg)] 
    
    ### stomach equations
    
    # define parameters and useful indices
    idx_ss <- map_c["stomach_s"]
    idx_sd <- map_c["stomach_d"]
    if(param_drug["type"]==0){           # neutral
       Y <- 1
    }else if(param_drug["type"]==1){     # acid
       Y <- 1 + 10^( ACAT_pH["stomach"] - param_drug["pKa"] )
    } else if (param_drug["type"]==2){   # base
       Y <- 1 + 10^( -ACAT_pH["stomach"] + param_drug["pKa"] ) 
    }
    Cs_st <- param_drug["Csint"] * Y
    Kd_st <- Kd1 * ( Cs_st - x[idx_sd]/(ACAT_v_lum["stomach"] + general_p["water_po"]))
    kst   <- 1/general_p["GET"] # [1/h]
    
    # define differential equations
    output_s   <- kst * x[idx_ss]
    output_d   <- kst * x[idx_sd]
    dx[idx_ss] <- -output_s - Kd_st * x[idx_ss]
    dx[idx_sd] <- -output_d + Kd_st * x[idx_ss]
    
    
    ### equations for all the small intestine transit compartments
    
    lsi <- length(comp_ent)
    length_small_int <- sum(ACAT_l[2:length(ACAT_l)])
    CO <- organ_bf["venous_blood"]  # cardiac output [L/h]
    
    for(i in 1:lsi){
      
      # define some indices
      idx_is <- map_c[comp_solid[i+1]]  # solid
      idx_id <- map_c[comp_diss[i+1]]   # dissolved
      idx_ie <- map_c[comp_ent[i]]      # enterocytes
      idx_ia <- map_c[comp_abs[i]]      # absorbed
      idx_ic <- map_c[comp_cl[i]]       # cleared
      
      # define some useful param
      if(param_drug["type"]==0){           # neutral
        Y <- 1
      }else if(param_drug["type"]==1){     # acid
        Y <- 1 + 10^( ACAT_pH[i+1] - param_drug["pKa"] )
      } else if (param_drug["type"]==2){   # base
        Y <- 1 + 10^( -ACAT_pH[i+1] + param_drug["pKa"] ) 
      }
      Cs_it <- param_drug["Csint"] * Y
      Kd_it <- Kd1 * ( Cs_st - x[idx_id]/ACAT_v_lum[i+1])
      kit   <- 1/(general_p["SITT"] * ACAT_l[i+1]/length_small_int) # [1/h]
      ka    <- 2 * param_drug["Peff"] / ACAT_l[i+1]                 # [1/h]
      Qe_i  <- ACAT_f_CO[i+1] * CO                                  # [L/h]
      
      # define differential equations
      dx[idx_is] <- +output_s - Kd_it * x[idx_is] - kit * x[idx_is]
      dx[idx_id] <- +output_d + Kd_it * x[idx_is] - kit * x[idx_id] - ka * x[idx_id]
      dx[idx_ie] <- ka * x[idx_id] - Qe_i * x[idx_ie] / ACAT_v_ent[i+1] - param_drug["CLent"] * x[idx_ie]/  ACAT_v_ent[i+1]
      dx[idx_ia] <- Qe_i * x[idx_ie] / ACAT_v_ent[i+1]
      dx[idx_ic] <- param_drug["CLent"] * x[idx_ie]/  ACAT_v_ent[i+1]
      
      output_s   <- kit * x[idx_is]
      output_d   <- kit * x[idx_id]
      
    }
    
    ### equations for the sink compartments
    dx[map_c["sink_s"]] <- output_s
    dx[map_c["sink_d"]] <- output_d
    
    ### output
    res  <- structure(dx, names=comp_ACAT_names)
    list(res)
    
  })
  
}







