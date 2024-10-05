# ODE Solvers #
library(tidyverse)
library(deSolve)
library(Hmisc)
library("doParallel")
library("RColorBrewer")
#==============================================================================================#

#==========================================================#
#   SECTION4.1 & 4.2 Simulation STUDY: 1 OR 2 TRUE MODES   #
#==========================================================#

# The ODE models we gonna use #

# 4.1_Simple model with a unique mode(S.Wang et al, JCGS, 2020) #
ode_sim1 <- function(t, y, params){
  # t		: time
  # y		: species
  # v		: parameters (ODE coefficients) 
  dy1dt <- 72/(36 + y[2]) - params[1]
  dy2dt <- params[2]*y[1] -1
  dxdt <- c(dy1dt, dy2dt)
  
  result <- list(dxdt)
  return(result)
}

# 4.2 A model with two true modes (S.Wang et al, JCGS, 2020)
ode_simbimod <- function(t, y, params){
  # t		: time
  # y		: species
  # v		: parameters (ODE coefficients) 
  dy1dt <- 72/(36 + y[2]) - abs(params[1])
  dy2dt <- params[2]*y[1] -1
  dxdt <- c(dy1dt, dy2dt)
  
  result <- list(dxdt)
  return(result)
}

#===============#
#DATA Generation#
#===============#

Simulator_ODE <- function (model, pars, tp)
  # Function to generate the simulated data based on models
  # model: ODE model Specification
  # pars: a list of parameters according to the model
  # tp	: Vector of time points
{
  d <- length(pars$init) 
  tt <- length(tp) 
  
  #= ODE SOLUTION using numerical solvers =#
  X_test <- ode(func = model, y = pars$init, times = tp, parms = pars$theta) 
  X_test <- X_test[,2:(d+1)]
  
  #= Generate Data =#
  traj = vector('list', d)
  
  for (i in 1:d) {
    
    traj[[i]] = X_test[,i] + rnorm(n= tt, mean= 0, sd = pars$sigma[i])
    
  }
  
  return(traj)	
}

###########
# Example #
###########

# Specifying a list of true parameters
theta.true <- c(2,1)
init.true <- c(7,-10)
sigma.true <- c(1,3)
sigma2.true <- sigma.true^2

true.pars <- list()
true.pars$theta <- theta.true
true.pars$init <- init.true
true.pars$sigma <- sigma.true
true.pars$sigma2 <- sigma2.true

# Specifying the timepoints
fromto <- seq(from = 0, to=60, length.out = 121)

ode_sim1_dat <- Simulator_ODE(model = ode_sim1, pars = true.pars, tp = fromto)
ode_sim1_dat_bimod <- Simulator_ODE(model = ode_simbimod, pars = true.pars, tp = fromto)




# STORAGE OF DATA for the simulation studies#
Simulation_Data <- list()
num.rep <- 50

set.seed(20230922)
for (j in 1:num.rep) {
  
  Simulation_Data[[j]] <- Simulator_ODE(model = ode_sim1, pars = true.pars, tp = fromto)
  
}
Simulation_Data[[1]]
