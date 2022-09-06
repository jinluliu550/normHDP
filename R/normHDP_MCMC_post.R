#' MCMC function to obtain posterior samples of mean expressions and dispersions
#' 
#' Due to the label switching property, we will need to run the following code to 
#' carry out posterior inferences on mean expressions and dispersions.
#' @param normHDP_output Output from the normHDP_mcmc function.
#' @param number_iter Total number of iterations.
#' @param thinning Default is to set thinning to 5; saves posterior estimates for every fifth iteration.
#' @param burn_in Default is to set to 3000; the first 3000 iterations are ignored.
#' @param Z The clustering estimate.
#' @param Y Input dataset.
#' @param partial.save.name Default is to set to NULL, an option to give a name to the output data if set auto.save = TRUE.
#' @param auto.save An option to save the output every 100 iterations. Default is to set to FALSE.
#' @param iter_update Controls for text, i.e.: if iter_update is set to 100, then the algorithm will output a message every 100 iterations.
#' @param num.cores Number of cores to run the function.
#' @return A list of posterior samples of unique parameters; acceptance probabilities and variance-covariance
#' for adaptive Metropolis Hastings; and all the fixed parameters required for the simulation of unique parameters.
#' @export
normHDP_mcmc_post <- function(normHDP_output,
                              burn_in,
                              thinning, 
                              number_iter, 
                              Z, 
                              auto.save = FALSE,
                              partial.save.name = NULL,
                              iter_update = 100,
                              Y,
                              num.cores = 4){
  
  ##--------------------------- Input parameters from previous MCMC -------------------------------------------
  
  ##-- Register cores
  num.cores <- normHDP_output$num.cores
  registerDoParallel(cl <- makeCluster(num.cores))
  
  
  G <- ncol(normHDP_output$mu_star_1_J_new)
  D <- normHDP_output$D
  C <- normHDP_output$C
  
  ##-- Only carry out MCMC for occupied clusters
  J <- length(unique(unlist(Z)))
  
  alpha_mu_2 <- normHDP_output$alpha_mu_2
  mu_star_1_J_initial <- normHDP_output$mu_star_1_J_initial
  phi_star_1_J_initial <- normHDP_output$phi_star_1_J_initial
  quadratic <- normHDP_output$quadratic
  adaptive_prop <- normHDP_output$adaptive_prop
  
  mu_star_1_J_output <- NULL
  phi_star_1_J_output <- NULL
  
  ##-- Posterior mean of capture efficiencies
  Beta_output <- normHDP_output$Beta_output
  Beta_posterior_mean <- lapply(1:D, function(d) rep(0, C[d]))
  
  for(d in 1:D){
    
    for(rep in 1:length(Beta_output)){
      Beta_posterior_mean[[d]] <- Beta_posterior_mean[[d]] + Beta_output[[rep]][[d]]
    }
    
    Beta_posterior_mean[[d]] <- Beta_posterior_mean[[d]]/length(Beta_output)
  }
  
  ##-- Posterior mean of alpha_phi_2
  alpha_phi_2 <- mean(normHDP_output$alpha_phi2_output)
  
  
  ##-- Posterior mean of b
  b_posterior_mean <- NULL
  b_trace <- normHDP_output$b_output
  for(i in 1:length(b_trace[[1]])){
    b_mean_each <- mean(unlist(lapply(b_trace, `[[`, i)))
    b_posterior_mean <- append(b_posterior_mean,b_mean_each)
  }
  
  ##--------------------- Measure acceptance probabilities ------------------------------------
  
  acceptance_prob_vec <- rep(0, number_iter-1)
  unique_count <- 0
  
  ##-- Initial values
  mu_star_1_J_new <- mu_star_1_J_initial
  phi_star_1_J_new <- phi_star_1_J_initial
  
  ##-- Variance-covariance structure of unique mean and disperison parameters
  covariance_unique_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  tilde_s_unique_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  mean_X_unique_new <- rep(list(rep(list(matrix(0,nrow=1,ncol=2)),G)),J)
  
  
  for(j in 1:J){
    
    tilde_s_unique_new[[j]] <- lapply(1:G, function(g){
      matrix(c(log(mu_star_1_J_new[j,g]),log(phi_star_1_J_new[j,g])),ncol=1)%*%
        matrix(c(log(mu_star_1_J_new[j,g]),log(phi_star_1_J_new[j,g])),nrow=1)}
    )
    
    mean_X_unique_new[[j]] <- lapply(1:G, function(g){
      matrix(c(log(mu_star_1_J_new[j,g]),log(phi_star_1_J_new[j,g])),nrow=1)
    })
    
  }
  
  ## Initial index
  output_index <- 0
  
  ##----------------------- MCMC of unique parameters ------------------------------------------
  
  ## Starting from iteration 2
  for(iter in 2:number_iter){
    
    # Starting value of the output index = 1
    # If the current iteration is greater than the burn in and divisible by the thinning index,
    # then the output index becomes output index + 1
    if(iter >= burn_in & iter%%thinning == 0){
      output_index <- output_index + 1
      update <- TRUE
    }else{
      update <- FALSE
    }
    
    if(iter %% iter_update == 0){
      print(paste("iter:", iter))
    }
    
    ## MCMC of unique mean and dispersion parameters
    unique_output_sim <- unique_parameters_mcmc(mu_star_1_J = mu_star_1_J_new, 
                                                phi_star_1_J = phi_star_1_J_new,
                                                mean_X_mu_phi = mean_X_unique_new,
                                                tilde_s_mu_phi = tilde_s_unique_new,
                                                Z = Z,
                                                b = b_posterior_mean, 
                                                alpha_phi_2 = alpha_phi_2,
                                                Beta = Beta_posterior_mean,
                                                alpha_mu_2 = alpha_mu_2,
                                                covariance = covariance_unique_new,
                                                Y = Y,
                                                iter_num = iter, 
                                                quadratic = quadratic,
                                                adaptive_prop = adaptive_prop)
    
    mu_star_1_J_new <- unique_output_sim$mu_star_1_J_new
    phi_star_1_J_new <- unique_output_sim$phi_star_1_J_new
    tilde_s_unique_new <- unique_output_sim$tilde_s_mu_phi_new
    mean_X_unique_new <- unique_output_sim$mean_X_mu_phi_new
    covariance_unique_new <- unique_output_sim$covariance_new
    unique_count <- unique_count + unique_output_sim$accept_count
    acceptance_prob_vec[iter-1] <- unique_count/((iter-1)*J*G)
    
    
    if(update == TRUE){
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      phi_star_1_J_output[[output_index]] <- phi_star_1_J_new
    }
    
    ## If running for long, can choose the option to output the final list every 100 iteration
    if(iter %% 100 == 0 && auto.save == TRUE){
      my_list <- list('mu_star_1_J_output' = mu_star_1_J_output, 
                      'phi_star_1_J_output' = phi_star_1_J_output, 
                      'mu_star_1_J_new' = mu_star_1_J_new,
                      'phi_star_1_J_new' = phi_star_1_J_new, 
                      'tilde_s_unique_new' = tilde_s_unique_new, 
                      'mean_X_unique_new' = mean_X_unique_new,
                      'covariance_unique_new' = covariance_unique_new,
                      'unique_count' = unique_count,
                      'acceptance_prob_vec' = acceptance_prob_vec,
                      'G' = G,
                      'J' = J,
                      'D' = D,
                      'C' = C,
                      'alpha_mu_2' = alpha_mu_2,
                      'quadratic'= quadratic,
                      'Beta_posterior_mean' = Beta_posterior_mean,
                      'b_posterior_mean' = b_posterior_mean,
                      'final_iter' = iter,
                      'output_index' = output_index,
                      'adaptive_prop' = adaptive_prop,
                      'auto.save' = auto.save,
                      'partial.save.name' = partial.save.name,
                      'iter_update' = iter_update,
                      'Y' = Y,
                      'Z' = Z,
                      'alpha_phi_2' = alpha_phi_2,
                      'num.cores' = num.cores)
      
      ##-- Save to workplace
      save(my_list, file=partial.save.name)
    }
  }
  
  # Return final output
  my_list <- list('mu_star_1_J_output' = mu_star_1_J_output, 
                  'phi_star_1_J_output' = phi_star_1_J_output, 
                  'mu_star_1_J_new' = mu_star_1_J_new,
                  'phi_star_1_J_new' = phi_star_1_J_new, 
                  'tilde_s_unique_new' = tilde_s_unique_new, 
                  'mean_X_unique_new' = mean_X_unique_new,
                  'covariance_unique_new' = covariance_unique_new,
                  'unique_count' = unique_count,
                  'acceptance_prob_vec' = acceptance_prob_vec,
                  'G' = G,
                  'J' = J,
                  'D' = D,
                  'C' = C,
                  'alpha_mu_2' = alpha_mu_2,
                  'quadratic'= quadratic,
                  'Beta_posterior_mean' = Beta_posterior_mean,
                  'b_posterior_mean' = b_posterior_mean,
                  'final_iter' = iter,
                  'output_index' = output_index,
                  'adaptive_prop' = adaptive_prop,
                  'auto.save' = auto.save,
                  'partial.save.name' = partial.save.name,
                  'iter_update' = iter_update,
                  'Y' = Y,
                  'Z' = Z,
                  'alpha_phi_2' = alpha_phi_2,
                  'num.cores' = num.cores)
  
  return(my_list)
}

