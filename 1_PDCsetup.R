#==================#
#    ODE solver    #
#==================#
ODEf <- function(theta, initial, model, tp){
  Init <- c(X1=initial[1], X2=initial[2])
  parameters <- c(theta1 = theta[1], theta2 = theta[2])
  output <- ode(y=Init, times=tp, func=model, parms=parameters, method = "lsoda")
  return(list(X1 = output[,2], X2 = output[,3]))
}

#==================================#
#    Cloned Likelihood Function    #
#==================================#
clone_Lik <- function(data, theta, initial, sigma,
                      model=ode_sim1, clone, tp=fromto){
  y1 <- data[[1]]
  y2 <- data[[2]]
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  output <- ODEf(theta, initial, model,tp)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  logL <- -sum((y1 - X1)^2/sigma1^2 + (y2 - X2)^2/sigma2^2)-length(y1)*log(sigma1^2)-length(y2)*log(sigma2^2)
  return(logL*clone)
}

#=============================================#
#    Prior Function for Bayesian framework    #
#=============================================#
odesim_prior <- function(theta, initial, sigma){
  # return the prior values with log scale
  
  # theta's ~ Normal(5, 5^2) 
  thetaprior <- dnorm(theta, mean = 5, sd=5, log = T) 
  # initial's ~ Normal(2,4^2)
  initprior <- dnorm(initial, mean = 2, sd = 4, log = T)
  # Sigma^2 ~ Inverse Gamma(1,1)
  sigprior <- invgamma::dinvgamma(sigma^2, shape = 1, rate = 1, log = T)
  
  lprior <- thetaprior + initprior + sigprior
  
  return(lprior)
}

#==========================#
# Gibbs move for sigma     #
#==========================#
GBKernel <- function(theta, initial, data, clone, model, tp){
  # To generate new sigma values from the full conditional dist'n of sigma^2
  # clone: the number of clones
  # model: ode model specification
  # tp: timepoint vector
  # output: sigma values, not squared.
  output <- ODEf(theta, initial, model = model, tp=tp)
  y1 <- data[[1]]
  y2 <- data[[2]]
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  temp1 <- rgamma(1, shape = 1+ length(y1)/2*clone, rate = 1 + clone*sum((data[[1]] - X1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2*clone, rate = 1 + clone*sum((data[[2]] - X2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}


#========================#
# MH-kernel within Gibbs #
#========================#
MHGibb <- function(particle, theta_cov, init_cov, data, clone, model,tp, alpha){
  # For both theta and init, we don't have a closed form of full conditional dist'n
  # To propose new theta and initial values given the new sigma values from Gibbsstep
  # MH within Gibbs
  # @particle = a list of parameters ($theta, $initial,)
  # theta_cov/init_cov = covariance structure of proposal dist'n 
  
  theta_old <- particle$theta
  #sigma_old <- particle$sigma
  init_old <- particle$initial
  #Gibbs step
  sigma_new <- GBKernel(theta = theta_old, initial = init_old, 
                        data = data, clone=clone, model=model, tp=tp)
  # Metropolis-Hastings
  
  # Simple proposal from multivariate normal dist'n: sample new theta
  theta_new <- MASS::mvrnorm(1, theta_old, theta_cov)
  init_new <- MASS::mvrnorm(1, init_old, init_cov)
  # compute numerator and denominator of the MH ratio
  num <- (clone_Lik(data=data, theta=theta_new, initial=init_new, model=model, sigma=sigma_new, clone)) * alpha +  
    odesim_prior(theta=theta_new, initial = init_new, sigma = sigma_new)
  den <- (clone_Lik(data=data, theta=theta_old, initial=init_old, model=model, sigma=sigma_new, clone)) * alpha +  
    odesim_prior(theta=theta_new, initial = init_new, sigma = sigma_new)
  
  # compute ratio
  ratio <- min(1, exp(num - den))
  
  # accept/reject step
  if(runif(1) < ratio){
    return(list(theta = theta_new, sigma= sigma_new, initial=init_new, accept = T))
  } else{
    return(list(theta = theta_old, sigma= sigma_new, initial=init_old, accept = F))
  }
  
  
  # MH(particle$theta, theta_cov, data, likelihood, prior, alpha[[r]])
}


#========================#
#   Functions for ASMC   #
#========================#

# log-sum-exponential evaluation of log(sum(w))
logsum <- function(logw){
  logmax = max(logw)
  log(sum(exp(logw-logmax)))+logmax
}

# # relative conditional effective sample size
# rCESS <- function(W, u, a, phi) {
#   logw <- a*u          # weight update
#   exp(2*logsum(log(W)+logw) - logsum(log(W)+2*logw)) - phi
# }

# effective sample size
rESS <- function(logW){
  NP <- length(logW)
  logWmax <- max(logW)
  logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(NP)
  return(exp(logRESS))
}

# relative conditional effective sample size
rCESS <- function(W, u, a, phi) {
  logw <- a*u          # weight update
  exp(2*logsum(log(W)+logw) - logsum(log(W)+2*logw)) - phi
}


# bisection function
# - recursive implementation of bisection algorithm
bisection <- function(low, high, W, u, phi ){
  
  mid <- (low+high)/2
  f.low <- rCESS( W, u, low, phi )
  f.mid <- rCESS( W, u, mid, phi )
  f.high <- rCESS( W, u, high, phi )
  
  if( f.low*f.high>0 )
    stop('Invalid endpoint for bisection.')
  
  try({if( low>=high )
    stop('bisection overlap')
    
    if( (abs(f.mid)<1e-10)||((high-low)/2<1e-10) )
      return( mid )
    if( (f.low*f.mid)<0 )
      return( bisection( low, mid, W, u, phi ) )
    if( (f.high*f.mid)<0 )
      return( bisection( mid, high, W, u, phi ) )
  })
  
  stop('bisection flawed')
}



#=====================================================================================#
#                   Function for Particle data cloning                                #
#=====================================================================================#

GDCASMC <- function(data, NP=500, theta_cov, init_cov, adaptive=F, alpha_f,
                    resamplethres = 0.5, model, clone=1, tp, phi
                    #init,  likelihood, prior, reference
){
  # @NP= the number of particles
  # @adaptive = T (automatic selection of annealing path) / F (predefined path)
  # @alpha_f = only Specifying for adaptive =F
  # @phi = recommendation: 0.99 or 0.999
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  #alpha <- list(alpha_f)
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  particles <- list()
  W <- list()
  
  # check if adaptive annealing scheme is being used
  #adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1               # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(alpha_f)
  }
  
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/NP,NP)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,NP)
  logw <- rep(0,NP)
  
  # initialize particles
  
  particles[[r]] <- lapply(1:NP, function(k) list(theta = c(sapply(result_1_propad47$particles[[1]], `[[`, "theta")[1,k], 
                                                            sapply(result_1_propad47$particles[[2]], `[[`, "theta")[2,k]),
                                                  sigma = c(runif(1, 0.3, 2), runif(1, 0.3, 4)),
                                                  initial = c(runif(1, -10, 10), runif(1, -15, 15)),
                                                  accept = F))
  # list(theta = c(runif(1, -4, 4), runif(1, 0, 2)),
  #sigma = c(runif(1, 0.3, 2), runif(1, 0.3, 4)),
  #initial = c(runif(1, -10, 10), runif(1, -15, 15)),
  #accept = F)
  # list(theta = c(runif(1, 0, 4), runif(1, 0, 3)),
  #      sigma = c(runif(1, 0.3, 3), runif(1, 0.3, 5)),
  #      initial = c(runif(1, 0, 15), runif(1, -15, 1)),
  #      accept = F)
  
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    # evaluate the log-likelihood for updating alpha
    incw <- rep(0, NP)   # incremental log-importance weights
    incw <- sapply(1:NP, function(k){
      #clone_Lik <- function(data, theta, initial, sigma,
      #                      model=ode_sim1, clone, tp=fromto){
      
      logL <- clone_Lik(data, 
                        theta = particles[[r-1]][[k]]$theta, 
                        initial = particles[[r-1]][[k]]$initial,
                        sigma = particles[[r-1]][[k]]$sigma, 
                        model = model, clone = clone, tp=tp)
      return(logL)
    })
    
    if(adaptive){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u=incw, phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    
    particles[[r]] <- lapply(particles[[r-1]],  MHGibb, theta_cov=theta_cov, init_cov=init_cov,
                             data=data, clone= clone, model=model, tp=tp, alpha = alpha[[r]])
    # accept_rate <- mean(sapply(particles[[r]], `[[`, "accept"))
    # cat("acceptance rate: ", accept_rate, "\n")
    
    # compute the ESS
    log_incremental_w <- alphaDiff * incw
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<resamplethres )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- sample.int(NP, prob = W[[r]], replace = TRUE) #alternative: systematic_resample( W[[r]] )
      particles[[r]] <-  particles[[r]][ancestors]
      W[[r]] <- rep(1/NP,NP)
      logW <- log(W[[r]])
      w <- rep(1,NP)
      logw <- rep(0,NP)
    }
  }
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W))
  
}



#=============================#
#           EXAMPLE           #
#=============================#
nclone <- 30
#alpha_vec <- (0:1000)/1000

theta_cov_temp <- diag(0.04^2, 2,2)/500 #diag(0.1^2, 2,2)
init_cov_temp <- diag(0.5^2, 2,2)/c(500,250) #diag(0.2^2,2,2)
#k=10/ init: diag(0.5^2, 2,2)/c(50,25)/ theta: diag(0.04^2, 2,2)/50



cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
start_time = Sys.time()
dcasmc30_foreach999 <- foreach(i=1:50, .packages = c("deSolve")) %dopar% {
  
  GDCASMC(data = Simulation_Data[[i]], NP=500, theta_cov = theta_cov_temp, init_cov = init_cov_temp,
          adaptive = T, alpha_f = alpha_vec, resamplethres = 0.5, model = ode_sim1, clone=nclone,
          tp=fromto,phi = 0.999)  
  
}
asmc30_Time_odesim999 <- Sys.time()-start_time  

parallel::stopCluster(cl)