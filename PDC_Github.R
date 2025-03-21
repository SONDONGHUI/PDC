#==================#
#    ODE solver    #
#==================#

ODEf <- function(theta, initial, model, tp){
  Init <- initial 
  parameters <- theta
  output <- ode(y=Init, times=tp, func=model, parms=parameters, method = "lsoda")[,-1]
  result <- lapply(seq_len(ncol(output)), function(i) output[, i])
  return(result)
}

#==================================#
#    Cloned Likelihood Function    #
#==================================#

clone_Lik_loop_new <- function(data, sigma, sol,clone){
  #@ data: List of Data,
  #@ sigma: Variance parameters from particles,
  #@ sol: ODE solutions from particles,
  #@ clone: The number of clones,
  #@ It returns log-likelihood * clone value.
  
  output <- sol
  logL <- 0
  for (i in 1:length(sigma)) {
    y <- data[[i]]
    sigmai <- sigma[i]
    X <- as.vector(output[[i]])
    logL <-  -sum((y - X)^2/sigmai^2)-length(y)*log(sigmai^2) + logL
  }
  #logS <- sum(logL)
  return(logL*clone)
}


#=============================================#
#    Prior Function for Bayesian framework    #
#=============================================#

odesim_prior68 <- function(theta, initial, sigma){
  #@ theta: ODE parameters,
  #@ initial: Initial conditions,
  #@ sigma: Variance parameters,
  #@ It evaluates the log-priors
  #@ This function is customized to be applied to prey68 data.
  #@ It is using hyperparameters based on the previous research results of Cao (2012).
  
  thetaprior <- sum(dnorm(theta, mean = c(3.9, 1.97, 4.3, 15.7, 
                                          0.11, 0.01, 0.152), 
                          sd=c(0.5, 0.25, 1.95, 2,
                               0.02, 0.14, 0.07),
                          log = T)) 
  initprior <- sum(dnorm(initial, mean = c(30.5, 1.5, 2.8, 0.2),
                         sd = c(4,0.25,0.5,0.2), log = T))
  sigprior <- sum(invgamma::dinvgamma(sigma^2, shape = 1, rate = 1, log = T))
  
  lprior <- thetaprior + initprior + sigprior
  
  return(lprior)
}


#==========================#
# Gibbs move for sigma     #
#==========================#

GBKernel_r_new <- function(theta, initial, data, sol, d, alpha, clone){
  
  #@ pi_r_invariant Gibbs Kernel to generate new candidate parameters.
  #@ theta/initial: particles of ODE and Initial parameters from the previous iteration
  #@ data: given data
  #@ sol: ODE solutions from the previous iteration
  #@ d: dimension of observed ODE components.
  #@ alpha: annealing parameter
  #@ clone: the number of clones.
  
  tempc <- numeric(d)
  #output <- ODEf(theta, initial, model = model, tp=tp)
  for (i in 1:d) {
    y <- data[[i]]
    X <- as.vector(sol[[i]])
    temp <- rgamma(1, shape = 1+ (length(y)/2)*clone*alpha,
                   rate = 1 + clone*alpha*sum((y - X)^2)/2)
    tempc[i] <- sqrt(1/temp)
  }
  
  return(tempc)
}


#========================#
# MH-kernel within Gibbs #
#========================#

MHGibb68_r_ref <- function(particle, refmu, refcov,
                           theta_cov, init_cov, 
                           adap_prop, mix,
                           data, clone, model,tp, alpha, d){
  #@ MH within Gibbs scheme incorporating reference distribution in subsection 2.4.
  #@ particle: particle objects
  #@ refmu: Multivariate normal mean vector for reference distribution in Eq (11) 
  #@ refcov: Multivariate normal covariance structure for reference distribution in Eq (11)   
  #@ adap_prop: T, PDC will use adaptive MCMC proposal
  #@ mix: mixture proportion (lower)
  #@ model: ODE model
  #@ tp: time points
  
  theta_old <- particle$theta
  #sigma_old <- particle$sigma
  init_old <- particle$initial
  
  d_theta <- length(theta_old) 
  d_init <- length(init_old) 
  d_cov <- d_theta  + d_init
  
  odesol_old <- particle$odesol
  #Gibbs
  sigma_new <- GBKernel_r_new(theta = theta_old, initial = init_old, alpha = alpha,
                              d=d, data = data,  clone=clone, sol = odesol_old)
  
  # Metropolis-Hastings
  likelihood_old <- clone_Lik_loop_new(data=data, 
                                       sol=odesol_old,
                                       sigma=sigma_new, 
                                       clone)
  
  # Sample new theta using Mixture distribution in Eq (10) of our paper 
  if(adap_prop==F){
    theta_new <- MASS::mvrnorm(1, theta_old, theta_cov)
    init_new <- MASS::mvrnorm(1, init_old, init_cov)
  } else {
    component <- sample(1:2, 1, replace = TRUE, prob = c(1-mix,mix))
    covlist_th <- list(theta_cov, diag((0.1)^2/(d_cov), d_theta))
    covlist_init <- list(init_cov, diag((0.1)^2/(d_cov), d_init))
    theta_new <- MASS::mvrnorm(1, theta_old, covlist_th[[component]]) 
    init_new <- MASS::mvrnorm(1, init_old, covlist_init[[component]]) 
    
  }
  odesol_new <- ODEf(theta = theta_new, initial = init_new, model = model, tp=tp)
  
  likelihood_new <- clone_Lik_loop_new(data=data, 
                                       sigma=sigma_new,
                                       sol=odesol_new,
                                       clone)  
  # compute numerator and denominator of the MH ratio
  num <- (likelihood_new) * alpha + 
    #proposal_mix(x=theta_old, mn= theta_new, covlist = covlist_th, mix = mix)+ 
    #proposal_mix(x=init_old, mn= init_new, covlist = covlist_init, mix = mix)+  
    odesim_prior68(theta=theta_new, initial = init_new, sigma = sigma_new)*alpha +
    dmvn(c(theta_new,init_new), mu= refmu, sigma = refcov,log = T)*(1-alpha)
  den <- (likelihood_old) * alpha + 
    #proposal_mix(x=theta_new, mn= theta_old, covlist = covlist_th, mix = mix)+ 
    #proposal_mix(x=init_new, mn= init_old, covlist = covlist_init, mix = mix)+
    odesim_prior68(theta=theta_old, initial = init_old, sigma = sigma_new)*alpha+
    dmvn(c(theta_old,init_old), mu= refmu, sigma = refcov,log = T)*(1-alpha)
  
  # compute ratio
  ratio <- min(1, exp(num - den))
  
  # accept/reject step
  if(runif(1) < ratio){
    return(list(theta = theta_new, sigma= sigma_new, initial=init_new, odesol=odesol_new, accept = T))
  } else{
    return(list(theta = theta_old, sigma= sigma_new, initial=init_old, odesol=odesol_old, accept = F))
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

GDCASMC68_r_ref <- function(data, NP=500, theta_cov, init_cov, 
                            refmu, refcov,
                            ip_option=F, ip_data,
                            adap_prop =F, mix=0.05,
                            adaptive=F, alpha_f,
                            resamplethres = 0.5, model, clone=1, tp, phi,d
                            #init,  likelihood, prior, reference
){
  # @NP= the number of particles
  # @adaptive = T (automatic selection of annealing path) / F (predefined path)
  # @alpha_f = only Specifying for adaptive =F
  # @phi = recommendation: 0.99 or 0.999
  # @refmu/refcov: reference distribution hyperparameter obtained from previous PDC with a smaller number of clones.
  # @ip_option: initial particle option. If T, initial particles are set to be the particles from previous PDC results.
  # @ip_data: if ip_option==T, previous PDC results data.
  # @adap_prop: adaptive MCMC proposal, or not.
  
  
  
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
  
  if(ip_option==F){
    ## Initialize weights
    # normalized weights
    # Particles were sampled from the same reference distribution. They have the equal weights.
    W[[r]] <- rep(1/NP,NP)
    logW <- log(W[[r]])
    
    # unnormalized weights
    w <- rep(1,NP)
    logw <- rep(0,NP)
    
    # initialize particles
    set.seed(364)
    particles[[r]] <- lapply(1:NP, function(k) {
      theta =  c(runif(1, 2,5), runif(1, 1,3), runif(1,3.5,6), runif(1,12,18),
                 runif(1,0,1), runif(1,0,0.5), runif(1,0,1))
      initial = c(runif(1, 25,35), runif(1,1,3), runif(1,1,5), runif(1,0,1))
      odesol = ODEf(theta = theta, initial = initial, model = model, tp=tp)
      
      list(theta = theta,
           sigma = c(runif(1,1,3),runif(1,1,3)),
           initial=initial,
           odesol = odesol,
           accept = F)})
    
    
  }else{
    ## Initialize weights from previous resultdata
    # normalized weights
    laststep_pre <- length(ip_data$alpha)
    W[[r]] <- rep(1/NP,NP)
    logW <- log(W[[r]])
    
    # unnormalized weights
    w <- rep(1,NP)
    logw <- rep(0,NP)
    
    # initialize particles from previous resultdata
    
    particles[[r]] <- ip_data$particles[[laststep_pre]]
  }
  
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    incw <- rep(0, NP)   # incremental log-importance weights
    incw <- sapply(1:NP, function(k){
      #clone_Lik <- function(data, theta, initial, sigma,
      #                      model=ode_sim1, clone, tp=fromto){
      
      logL <- clone_Lik_loop_new(data, 
                                 sol = particles[[r-1]][[k]]$odesol,
                                 sigma = particles[[r-1]][[k]]$sigma, 
                                 clone = clone)
      logprior <- sum(dnorm(particles[[r-1]][[k]]$theta, mean = c(3.9, 1.97, 4.3, 15.7, 
                                                                  0.11, 0.01, 0.152), 
                            sd=c(0.5, 0.25, 1.95, 2,
                                 0.02, 0.14, 0.07),
                            log = T))+
        sum(dnorm(particles[[r-1]][[k]]$initial, mean = c(30.5, 1.5, 2.8, 0.2),
                  sd = c(4,0.25,0.5,0.2), log = T))
      
      ref <- dmvn(c(particles[[r-1]][[k]]$theta,particles[[r-1]][[k]]$initial),
                  mu=refmu, sigma = refcov, log = T)
      
      return(logL+logprior-ref)
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
    
    if(adap_prop==F){
      cov_theta = theta_cov
      cov_init = init_cov
    } else {
      theta_temp <- t(sapply(particles[[r-1]],'[[', "theta"))
      init_temp <- t(sapply(particles[[r-1]],'[[', "initial"))
      dcov <- ncol(theta_temp)+ncol(init_temp)
      weight_temp <- W[[r-1]]
      
      cov_theta = ((2.38)^2)*(cov.wt(theta_temp, wt = weight_temp)$cov)/(dcov)
      cov_init = ((2.38)^2)*(cov.wt(init_temp, wt = weight_temp)$cov)/(dcov)
    }
    particles[[r]] <- lapply(particles[[r-1]],  MHGibb68_r_ref, 
                             theta_cov=cov_theta, init_cov=cov_init, 
                             refmu=refmu, refcov=refcov,
                             adap_prop=adap_prop, mix=mix,
                             data=data, clone= clone, model=model, tp=tp, alpha = alpha[[r]],d=d)
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
#           Result            #
#=============================#

Summary_asmc <- function(resultlist){
  # @ It extract the weighted particles from the last Annealing iteration.
  
  R <- length(resultlist$particles)
  
  theta <- t(sapply(resultlist$particles[[R]],'[[', "theta"))
  
  sigma <- t(sapply(resultlist$particles[[R]],'[[', "sigma"))
  
  init <- t(sapply(resultlist$particles[[R]],'[[', "initial"))
  weight <- resultlist$W[[R]]
  
  summarydat <- as.data.frame(cbind(theta,sigma,init,weight))
  #colnames(summarydat) <- c('theta1', 'theta2', 'sig1', 'sig2', 'init1', 'init2', 'weight')
  return(summarydat)
} 

summary_stat_asmc <- function(data, clone){
  # @ Inference Results based on the weighted particles (weighted mean & variance/ upper and lower bounds).
  # @data: from Summary_asmc output.
  # @clone: the number of clones.
  
  weight <- data$weight
  data_reduced <- data[,-ncol(data)]
  wmean <- colSums(data_reduced * weight )
  wsd <- apply(data_reduced, MARGIN = 2, 
               FUN = function(x) sqrt(wtd.var((x), weights = weight, normwt = T)*clone))  
  
  LB <- wmean -1.96*wsd
  UB <- wmean +1.96*wsd
  #Cover <- (LB<= truepar) & (UB>=truepar)
  
  summarystat <- rbind(wmean, wsd, LB, UB)
  rownames(summarystat) <- c('wmean', 'wsd', 'LB','UB')
  
  return(summarystat)
}


#=============================#
#           EXAMPLE           #
#=============================#
mle_thex0_k5 <- c(summary_stat_asmc(Summary_asmc(result_pp68k5_r_999_500_ref), clone = 5)[1,-c(8,9)])
cov_thex0_k5 <- (cov.wt(Summary_asmc(result_pp68k5_r_999_500_ref)[,-c(8,9,14)],
                        wt = Summary_asmc(result_pp68k5_r_999_500_ref)$weight)$cov)


result_pp68k20_r_999_500_ref <- GDCASMC68_r_ref(data = preydtset68, NP=500, 
                                                refmu = mle_thex0_k5, refcov = cov_thex0_k5,
                                                theta_cov = theta_cov_temp68/100, init_cov = init_cov_temp68/100,
                                                ip_option = T, ip_data = result_pp68k5_r_999_500_ref,
                                                adap_prop = T,
                                                adaptive = T, alpha_f = alpha_vec, 
                                                resamplethres = 0.5, model = ppmodel68, clone=20,
                                                tp=preytp68, phi = 0.999, d=2)  

summary_stat_asmc(Summary_asmc(result_pp68k20_r_999_500_ref), clone = 20)
