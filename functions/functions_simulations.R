### Functions used for simulating the set of ODE
# for now just single and multiple boluses in one compartment are supported
# consider in the future to do a general function, supporting single/multiple doses simultaneously in different compartments



### function for simple schedule ------------------------------------------------------------
# define a simple schedule of one dose administered each time_interval for n_rep 
# and the simulation finish at time_interval_end
defineSchedule <- function(dose, time_begin, time_interval, n_rep, time_interval_end){
  
  if(n_rep>1){
    dose_v   <- rep(dose, n_rep)
    time_v   <- seq(from=time_begin, to=time_begin + time_interval*(n_rep-1), by=time_interval)
    time_end <- time_v[length(time_v)] + time_interval_end
    dose_v   <- c(dose_v, 0)
    time_v   <- c(time_v, time_end)
  }else if(n_rep == 1){
    dose_v <- c(dose, 0)
    time_v <- c(0, time_interval + time_interval_end)
  }
  
  
  return(list(dose = dose_v, time = time_v))
  
}


### simulation ----------------------------------------------------------------------
# ode functions there is an "event"
simulation <- function(func, schedule, delta.t, comp_dose, param.PBPK){
  
  comp.names <- c(param.PBPK$comp_PBPK_names, param.PBPK$comp_ACAT_names)
  lo         <- length(comp.names)
  n_admin    <- length(schedule$dose) - 1  # the last element is supposed to be the end of the simulation (without dosing)
  ini.SA     <- structure(numeric(lo), names=comp.names)
  system.out.g <- c() 
  #t.g <- c()
  
  for(i in 1:n_admin){
    
    ini.SA[comp_dose] <- ini.SA[comp_dose] + schedule$dose[i]
    t <- seq(schedule$time[i], schedule$time[i+1], by=delta.t)
    
    system.out <- data.frame(lsodes(y = ini.SA, times = t, func = func, 
                                 parms = param.PBPK))
    
    # set the initial conditions for the next iteration
    idx.end <- length(system.out$time)
    ini.SA  <- unlist(as.vector(system.out[idx.end,2:(lo+1)]))
    
    # bind general solution
    system.out.g <- rbind(system.out.g, system.out[1:(idx.end-1),])
    #t.g <- rbind(t.g, t[1:(idx.end-1)])
    
  }
  
  return(system.out.g)
  
}
  







  
  
  
