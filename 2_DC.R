#======================================================================================#
#======================================================================================#
# Basic Data Cloning + Simple ODE: to compare GDCASMC and DC have same posterior dist'n#
#======================================================================================#
#======================================================================================#

DC <- function(data, iter=500000, theta_cov, init_cov, 
               #adaptive=F, alpha_f, resamplethres = 0.5, 
               model, clone=1, tp, alpha
               #init,  likelihood, prior, reference
){
  
  # initialize storage list
  particles_dc <- list()
  
  r <- 1               # SMC iteration
  
  
  # initialize particles
  set.seed(364)
  particles_dc[[r]] <- list(theta = c(runif(1, -4, 4), runif(1, 0, 4)),
                            sigma = c(runif(1, 0.3, 2), runif(1, 0.3, 4)),
                            initial = c(runif(1, -1, 10), runif(1, -15, 1)),
                            accept = F)
  
  for(r in 1:(iter-1))   # repeat this step if the tempering parameter is less than 1
  {
    #cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    
    #cat("annealing parameter:",alpha[[r]],"\n")
    
    
    particles_dc[[r]] <- MHGibb(particle = particles_dc[[r-1]], 
                                theta_cov = theta_cov, init_cov=init_cov,
                                data=data, clone=clone, model=model, tp=tp, alpha = alpha) 
    
    if(r%%10000 == 0) cat("iteration:",r,"\n")
  }
  
  return(particles_dc)
  
}


#===========================#
#          EXAMPLE          #
#===========================#

theta_cov_dc12 <- diag(0.04^2, 2,2)/500 #diag(0.1^2, 2,2)
init_cov_dc12 <- diag(0.5^2, 2,2)/c(500,250)


cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

start_time = Sys.time()
DC12_foreach_bimod <- foreach(i = 1:50, .packages = c("deSolve")) %dopar% {
  
  DC(data = Simulation_Data_bimod[[i]], iter = iter, 
     theta_cov = theta_cov_dc12, init_cov = init_cov_dc12,
     model = ode_simbimod, clone = 12, tp=fromto, alpha = 1)
  
}
DC12bimo_Time_odesim <- Sys.time()-start_time #4.1h
parallel::stopCluster(cl)