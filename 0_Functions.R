library(tidyverse)
library(deSolve)
library(reshape2)
library(plotly)
library(gridExtra)
library(see)
library(cowplot)
library(doFuture)

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