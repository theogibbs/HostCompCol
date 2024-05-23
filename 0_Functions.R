library(tidyverse)
library(deSolve)
library(reshape2)
library(plotly)
library(gridExtra)

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
  timeseries <- as.data.frame(ode(inistate, times, fn, pars))  
  return(timeseries)
}

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

GetEig <- function(J) {
  eigs <- eigen(J, only.values = TRUE)$values
  re_eigs <- Re(eigs)
  max_eig <- max(re_eigs)
  return(max_eig)
}




