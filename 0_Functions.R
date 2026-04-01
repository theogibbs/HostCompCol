library(tidyverse)
library(deSolve)
library(reshape2)
library(plotly)
library(gridExtra)
library(see)
library(cowplot)
library(doFuture)
library(profvis)
library(ggnewscale)
library(ggpattern)

# returns the derivatives from the current state and
# model parameters to integrate the ODEs
CoralMutPathDynamics <- function(time, state, pars) {
  dstatedt <- with(pars, {
    
    h <- state[1]
    p <- state[2]
    m <- state[3]
    
    dhdt <- ((h - m) * ch + chm * m) * (1 - h) - dh * h - dhp * p
    dpdt <- cp * (h - p) * p - (dp + dh + dhp) * p
    dmdt <- cm * (h - p - m) * m - (dm + cp * p + dh) * m
    
    dstatedt <- list(c(dhdt, dpdt, dmdt))
    return(dstatedt)
  })
  return(dstatedt)
}

# integrates the dynamics using deSolve
IntegrateDynamics <- function(inistate, pars, endtime, timestep, fn){
  times <- seq(0, endtime, by = timestep)
  timeseries <- as.data.frame(ode(inistate, times, fn, pars,
                                  method = "ode45"))  
  return(timeseries)
}

# Just returns the value of the ODE at some final time points
EndDynamics <- function(inistate, pars, endtime, timestep, timelength, fn){
  times <- c(0, seq(endtime - timestep*timelength, endtime, length.out = timelength+1))
  timeseries <- as.data.frame(ode(inistate, times, fn, pars,
                                  method = "ode45"))
  return(timeseries)
}

# returns the equation for the host so we can solve
# it numerically using the current parameters
PredEq <- function(h, pars) {
  ret <- with(pars, {
    X <- (dp + dh + dhp) / cp
    Y <- (dp + dhp - dm) / cm
    
    m <- X + Y - cp / cm * h
    p <- h - X
    ret <- ((h - m) * ch + chm * m) * (1 - h) - (dh * h + dhp * p)
    return(ret)
  })
  return(ret)
}

# returns the values of the pathogen and mutualist
# given a host value and the parameters
GetPM <- function(h, pars) {
  ret <- with(pars, {
    X <- (dp + dh + dhp) / cp
    Y <- (dp + dhp - dm) / cm
    
    p <- h - X
    m <- X + Y - cp / cm * h
    ret <- c(p, m)
    names(ret) <- c("p", "m")
    return(ret)
  })
  return(ret)
}

# builds the Jacobian of the model as in
# Appendix 1 of the main text
BuildJacobian <- function(h, p, m, pars) {
  J <- with(pars, {
    J <- matrix(0, 3, 3)
    
    J[1,1] <- ch - dh - (chm - ch) * m - 2 * ch * h
    J[1,2] <- - dhp
    J[1,3] <- (chm - ch) * (1 - h)

    J[2,1] <- cp * p
    J[2,2] <- - cp * p
    
    J[3,1] <- cm * m
    J[3,2] <- - (cm + cp) * m
    J[3,3] <- -cm * m
    
    return(J)
  })
  return(J)
}

# helper function to find the eigenvalues of the Jacobian
# and return the one with largest real part
GetEig <- function(J) {
  eigs <- eigen(J, only.values = TRUE)$values
  re_eigs <- Re(eigs)
  max_eig <- max(re_eigs)
  return(max_eig)
}

HM_LimitCycle <- function(h, p, m, pars) {
  ret <- with(pars, {
    LHS = -ch * h^2 + (ch - dh) * h + (chm - ch) * (1 - h) * m
    RHS = dhp * p
    # h_p0 = dhp * ((dh + dp + dhp) / cp)
    # h_p1 = dh - ch + dhp
    # h_p2 = ch
    # h_p = h_p0 + h_p1 * h + h_p2 * h^2
    # 
    # m_p0 = 0
    # m_p1 = (chm - ch) * (1 + ((dh + dm) / cp) + (1 + (cp / cm)) * ((dh + dp) / cp))
    # m_p2 = (chm - ch) * cm / cp
    # m_p = m_p0 + m_p1 * m + m_p2 * m^2
    
    # ret <- m_p - h_p
    ret <- LHS - RHS
    return(ret)
  })
  return(ret)
}

BuildKernel <- function(SpatialStructure, N, kern_length) {
  
  neighborhoods <- vector("list", N)
  
  if(SpatialStructure == "1D Lattice") {
    for(i in 0:(N-1)) {
      neighbors <- seq(i-kern_length, i+kern_length, by = 1)
      neighbors <- neighbors %% N
      neighbors <- unique(neighbors)
      neighborhoods[[i+1]] <- neighbors + 1
    }
  }
  return(neighborhoods)
}

GetInitialRates <- function(cur_host, cur_path, cur_mut, pars) {
  
  host_col_host <- pars$ch * (cur_host + cur_path)
  host_col_mut  <- pars$chm * cur_mut
  
  host_death <- pars$dh * (cur_host + cur_path + cur_mut)
  host_death_path <- pars$dhp * cur_path
  
  path_col <- pars$cp * cur_path
  path_death <- pars$dp * cur_path
  
  mut_col <- pars$cm * cur_mut
  mut_death <- pars$dm * cur_mut
  
  event_vec <- c("HostColHost", "HostColMut",
                 "HostDeath", "HostDeathPath",
                 "PathCol", "PathDeath",
                 "MutCol", "MutDeath")
  
  rate_vec <- c(host_col_host, host_col_mut,
                host_death, host_death_path,
                path_col, path_death,
                mut_col, mut_death)
  
  new_rates <- data.frame(Event = event_vec, Rate = rate_vec,
                          CumProb = cumsum(rate_vec) / sum(rate_vec))
  
  return(new_rates)
}

RunSpatialSim <- function(pars, ini_state, end_step, collect_spatial_steps, collect_freq_steps, print_time = TRUE) {
  
  cur_state <- ini_state
  
  ini_host <- sum(ini_state == "Host")
  ini_path <- sum(ini_state == "Pathogen")
  ini_mut <- sum(ini_state == "Mutualist")
  
  cur_host <- ini_host
  cur_path <- ini_path
  cur_mut <- ini_mut
  
  cur_time <- 0
  
  out_data <- data.frame()
  out_freq <- matrix(0, nrow = length(collect_freq_steps), ncol = 5)
  
  cur_ind <- 0
  start_time <- Sys.time()
  
  for(cur_step in 1:end_step) {
    
    cur_samp <- runif(1, min = 0, max = 1)
    if(cur_step == 1) {
      cur_rates <- GetInitialRates(cur_host, cur_path, cur_mut, pars)
    }
    
    cur_probs <- cur_rates %>%
      mutate(Success = CumProb > cur_samp)
    
    filt_probs <- cur_probs %>%
      filter(Success == TRUE)
    
    cur_event <- filt_probs[1,]
    cur_time <- cur_time + log(1 / cur_samp) / sum(cur_rates$Rate)
    
    if(cur_event$Event == "HostColHost") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Host", "Pathogen"))
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Empty") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Host"
        cur_host <- cur_host + 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] + pars$dh
        
      }
      
    } else if(cur_event$Event == "HostColMut") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Mutualist"))
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Empty") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Host"
        cur_host <- cur_host + 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] + pars$dh
        
      }
      
    } else if(cur_event$Event == "HostDeath") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Host", "Pathogen", "Mutualist"))
      cur_site <- sample(cur_sites$Site, size = 1)
      dead_type <- cur_state$Occ[cur_state$Site == cur_site]
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Empty"
      if(dead_type == "Host") {
        cur_host <- cur_host - 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
        
      } else if(dead_type == "Pathogen") {
        cur_path <- cur_path - 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
        cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] - pars$dhp
        cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] - pars$cp
        cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] - pars$dp
        
      } else if(dead_type == "Mutualist") {
        cur_mut <- cur_mut - 1
        cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] - pars$chm
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
        cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] - pars$cm
        cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] - pars$dm
        
      }
      
    } else if(cur_event$Event == "HostDeathPath") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Pathogen"))
      cur_site <- sample(cur_sites$Site, size = 1)
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Empty"
      cur_path <- cur_path - 1
      
      cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
      cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
      cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] - pars$dhp
      cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] - pars$cp
      cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] - pars$dp
      
    } else if(cur_event$Event == "PathCol") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Pathogen"))
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Host") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Pathogen"
        cur_path <- cur_path + 1
        cur_host <- cur_host - 1
        
        cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] + pars$dhp
        cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] + pars$cp
        cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] + pars$dp
        
      } else if(cur_state$Occ[cur_state$Site == cur_target] == "Mutualist") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Pathogen"
        cur_path <- cur_path + 1
        cur_mut <- cur_mut - 1
        
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] + pars$dhp
        cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] + pars$cp
        cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] + pars$dp
        
        cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] - pars$chm
        cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] - pars$cm
        cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] - pars$dm
        
      }
      
    } else if(cur_event$Event == "PathDeath") {
      
      cur_sites <- cur_state %>%
        filter(Occ == "Pathogen")
      cur_site <- ifelse(length(cur_sites$Site) > 1, sample(cur_sites$Site, size = 1), cur_sites$Site)
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Host"
      cur_path <- cur_path - 1
      cur_host <- cur_host + 1
      
      cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] - pars$dhp
      cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] - pars$cp
      cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] - pars$dp
      
    } else if(cur_event$Event == "MutCol") {
      
      cur_sites <- cur_state %>%
        filter(Occ == "Mutualist")
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Host") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Mutualist"
        cur_mut <- cur_mut + 1
        cur_host <- cur_host - 1
        
        cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] + pars$chm
        cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] + pars$cm
        cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] + pars$dm
        
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
        
      }
    } else if(cur_event$Event == "MutDeath") {

      cur_sites <- cur_state %>%
        filter(Occ == "Mutualist")
      
      cur_site <- ifelse(length(cur_sites$Site) > 1, sample(cur_sites$Site, size = 1), cur_sites$Site)
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Host"
      cur_mut <- cur_mut - 1
      cur_host <- cur_host + 1
      
      cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] - pars$chm
      cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] - pars$cm
      cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] - pars$dm
      
      cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
      
    }
    
    cur_rates$CumProb <- cumsum(cur_rates$Rate) / sum(cur_rates$Rate)
    
    if(cur_step %in% collect_spatial_steps) {
      cur_data <- cur_state
      cur_data$Time <- cur_time
      cur_data$Step <- cur_step
      out_data <- rbind(out_data, cur_data)
    }
    
    if(cur_step %in% collect_freq_steps) {
      cur_ind <- cur_ind+1
      out_freq[cur_ind,] <- c(cur_step, cur_time, (cur_host + cur_path + cur_mut) / pars$N, cur_path / pars$N, cur_mut / pars$N)
    }
    
    #check_host <- sum(cur_state == "Host")
    #check_path <- sum(cur_state == "Pathogen")
    #check_mut <- sum(cur_state == "Mutualist")
    #
    #if(!prod(check_host == cur_host, check_path == cur_path, check_mut == cur_mut)) {
    #  print("UH OH")
    #  print(check_host)
    #  print(cur_host)
    #  print(check_path)
    #  print(cur_path)
    #  print(check_mut)
    #  print(cur_mut)
    #  
    #  print(cur_event)
    #  print(cur_step)
    #  
    #  print(cur_site)
    #  print(cur_sites)
    #  print(cur_state)
    # break
    #}
    
    if((cur_host + cur_path + cur_mut) == 0) {
      print("       everyone died :(")
      cur_ind <- cur_ind+1
      out_freq[cur_ind,] <- c(cur_step, cur_time, (cur_host + cur_path + cur_mut) / pars$N, cur_path / pars$N, cur_mut / pars$N)
      break
    }
    
  }
  #})
  
  elaps_time <- Sys.time() - start_time
  if(print_time) print(elaps_time)
  
  out_freq <- as.data.frame(out_freq)
  colnames(out_freq) <- c("Step", "Time", "Host", "Pathogen", "Mutualist")
  
  return(list(freq = out_freq, spatial = out_data))
  
}

GetInitialRates <- function(cur_host, cur_path, cur_mut, pars) {
  
  host_col_host <- pars$ch * (cur_host + cur_path)
  host_col_mut  <- pars$chm * cur_mut
  
  host_death <- pars$dh * (cur_host + cur_path + cur_mut)
  host_death_path <- pars$dhp * cur_path
  
  path_col <- pars$cp * cur_path
  path_death <- pars$dp * cur_path
  
  mut_col <- pars$cm * cur_mut
  mut_death <- pars$dm * cur_mut
  
  event_vec <- c("HostColHost", "HostColMut",
                 "HostDeath", "HostDeathPath",
                 "PathCol", "PathDeath",
                 "MutCol", "MutDeath")
  
  rate_vec <- c(host_col_host, host_col_mut,
                host_death, host_death_path,
                path_col, path_death,
                mut_col, mut_death)
  
  new_rates <- data.frame(Event = event_vec, Rate = rate_vec,
                          CumProb = cumsum(rate_vec) / sum(rate_vec))
  
  return(new_rates)
}

RunSpatialSim <- function(pars, ini_state, end_step, collect_spatial_steps, collect_freq_steps, print_time = TRUE) {
  
  cur_state <- ini_state
  
  ini_host <- sum(ini_state == "Host")
  ini_path <- sum(ini_state == "Pathogen")
  ini_mut <- sum(ini_state == "Mutualist")
  
  cur_host <- ini_host
  cur_path <- ini_path
  cur_mut <- ini_mut
  
  cur_time <- 0
  
  out_data <- data.frame()
  out_freq <- matrix(0, nrow = length(collect_freq_steps), ncol = 5)
  
  cur_ind <- 0
  start_time <- Sys.time()
  
  for(cur_step in 1:end_step) {
    
    cur_samp <- runif(1, min = 0, max = 1)
    if(cur_step == 1) {
      cur_rates <- GetInitialRates(cur_host, cur_path, cur_mut, pars)
    }
    
    cur_probs <- cur_rates %>%
      mutate(Success = CumProb > cur_samp)
    
    filt_probs <- cur_probs %>%
      filter(Success == TRUE)
    
    cur_event <- filt_probs[1,]
    cur_time <- cur_time + log(1 / cur_samp) / sum(cur_rates$Rate)
    
    if(cur_event$Event == "HostColHost") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Host", "Pathogen"))
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Empty") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Host"
        cur_host <- cur_host + 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] + pars$dh
        
      }
      
    } else if(cur_event$Event == "HostColMut") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Mutualist"))
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Empty") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Host"
        cur_host <- cur_host + 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] + pars$dh
        
      }
      
    } else if(cur_event$Event == "HostDeath") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Host", "Pathogen", "Mutualist"))
      cur_site <- sample(cur_sites$Site, size = 1)
      dead_type <- cur_state$Occ[cur_state$Site == cur_site]
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Empty"
      if(dead_type == "Host") {
        cur_host <- cur_host - 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
        
      } else if(dead_type == "Pathogen") {
        cur_path <- cur_path - 1
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
        cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] - pars$dhp
        cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] - pars$cp
        cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] - pars$dp
        
      } else if(dead_type == "Mutualist") {
        cur_mut <- cur_mut - 1
        cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] - pars$chm
        cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
        cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] - pars$cm
        cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] - pars$dm
        
      }
      
    } else if(cur_event$Event == "HostDeathPath") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Pathogen"))
      cur_site <- sample(cur_sites$Site, size = 1)
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Empty"
      cur_path <- cur_path - 1
      
      cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
      cur_rates$Rate[cur_rates$Event == "HostDeath"] <- cur_rates$Rate[cur_rates$Event == "HostDeath"] - pars$dh
      cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] - pars$dhp
      cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] - pars$cp
      cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] - pars$dp
      
    } else if(cur_event$Event == "PathCol") {
      
      cur_sites <- cur_state %>%
        filter(Occ %in% c("Pathogen"))
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Host") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Pathogen"
        cur_path <- cur_path + 1
        cur_host <- cur_host - 1
        
        cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] + pars$dhp
        cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] + pars$cp
        cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] + pars$dp
        
      } else if(cur_state$Occ[cur_state$Site == cur_target] == "Mutualist") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Pathogen"
        cur_path <- cur_path + 1
        cur_mut <- cur_mut - 1
        
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
        cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] + pars$dhp
        cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] + pars$cp
        cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] + pars$dp
        
        cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] - pars$chm
        cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] - pars$cm
        cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] - pars$dm
        
      }
      
    } else if(cur_event$Event == "PathDeath") {
      
      cur_sites <- cur_state %>%
        filter(Occ == "Pathogen")
      cur_site <- ifelse(length(cur_sites$Site) > 1, sample(cur_sites$Site, size = 1), cur_sites$Site)
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Host"
      cur_path <- cur_path - 1
      cur_host <- cur_host + 1
      
      cur_rates$Rate[cur_rates$Event == "HostDeathPath"] <- cur_rates$Rate[cur_rates$Event == "HostDeathPath"] - pars$dhp
      cur_rates$Rate[cur_rates$Event == "PathCol"] <- cur_rates$Rate[cur_rates$Event == "PathCol"] - pars$cp
      cur_rates$Rate[cur_rates$Event == "PathDeath"] <- cur_rates$Rate[cur_rates$Event == "PathDeath"] - pars$dp
      
    } else if(cur_event$Event == "MutCol") {
      
      cur_sites <- cur_state %>%
        filter(Occ == "Mutualist")
      cur_site <- sample(cur_sites$Site, size = 1)
      neighbors <- pars$host_kern[[cur_site]]
      cur_target <- sample(neighbors, size = 1)
      
      if(cur_state$Occ[cur_state$Site == cur_target] == "Host") {
        cur_state$Occ[cur_state$Site == cur_target] <- "Mutualist"
        cur_mut <- cur_mut + 1
        cur_host <- cur_host - 1
        
        cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] + pars$chm
        cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] + pars$cm
        cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] + pars$dm
        
        cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] - pars$ch
        
      }
    } else if(cur_event$Event == "MutDeath") {

      cur_sites <- cur_state %>%
        filter(Occ == "Mutualist")
      
      cur_site <- ifelse(length(cur_sites$Site) > 1, sample(cur_sites$Site, size = 1), cur_sites$Site)
      
      cur_state$Occ[cur_state$Site == cur_site] <- "Host"
      cur_mut <- cur_mut - 1
      cur_host <- cur_host + 1
      
      cur_rates$Rate[cur_rates$Event == "HostColMut"] <- cur_rates$Rate[cur_rates$Event == "HostColMut"] - pars$chm
      cur_rates$Rate[cur_rates$Event == "MutCol"] <- cur_rates$Rate[cur_rates$Event == "MutCol"] - pars$cm
      cur_rates$Rate[cur_rates$Event == "MutDeath"] <- cur_rates$Rate[cur_rates$Event == "MutDeath"] - pars$dm
      
      cur_rates$Rate[cur_rates$Event == "HostColHost"] <- cur_rates$Rate[cur_rates$Event == "HostColHost"] + pars$ch
      
    }
    
    cur_rates$CumProb <- cumsum(cur_rates$Rate) / sum(cur_rates$Rate)
    
    if(cur_step %in% collect_spatial_steps) {
      cur_data <- cur_state
      cur_data$Time <- cur_time
      cur_data$Step <- cur_step
      out_data <- rbind(out_data, cur_data)
    }
    
    if(cur_step %in% collect_freq_steps) {
      cur_ind <- cur_ind+1
      out_freq[cur_ind,] <- c(cur_step, cur_time, (cur_host + cur_path + cur_mut) / pars$N, cur_path / pars$N, cur_mut / pars$N)
    }
    
    #check_host <- sum(cur_state == "Host")
    #check_path <- sum(cur_state == "Pathogen")
    #check_mut <- sum(cur_state == "Mutualist")
    #
    #if(!prod(check_host == cur_host, check_path == cur_path, check_mut == cur_mut)) {
    #  print("UH OH")
    #  print(check_host)
    #  print(cur_host)
    #  print(check_path)
    #  print(cur_path)
    #  print(check_mut)
    #  print(cur_mut)
    #  
    #  print(cur_event)
    #  print(cur_step)
    #  
    #  print(cur_site)
    #  print(cur_sites)
    #  print(cur_state)
    # break
    #}
    
    if((cur_host + cur_path + cur_mut) == 0) {
      print("       everyone died :(")
      cur_ind <- cur_ind+1
      out_freq[cur_ind,] <- c(cur_step, cur_time, (cur_host + cur_path + cur_mut) / pars$N, cur_path / pars$N, cur_mut / pars$N)
      break
    }
    
  }
  #})
  
  elaps_time <- Sys.time() - start_time
  if(print_time) print(elaps_time)
  
  out_freq <- as.data.frame(out_freq)
  colnames(out_freq) <- c("Step", "Time", "Host", "Pathogen", "Mutualist")
  
  return(list(freq = out_freq, spatial = out_data))
  
}

