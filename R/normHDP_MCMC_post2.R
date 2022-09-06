
#' Extension of normHDP_mcmc_post
#' 
#' To obtain more samples then the ones obtained from normHDP_mcmc_post.
#' 
#' @param normHDP_post_output Output from function normHDP_mcmc_post.
#' @param number_iter Total number of iterations.
#' @return A list of posterior samples of unique parameters; acceptance probabilities and variance-covariance
#' for adaptive Metropolis Hastings; and all the fixed parameters required for the simulation of unique parameters.
#' @export
normHDP_mcmc_post2 <- function(normHDP_post_output,
                               number_iter){
  
  num.cores <- normHDP_post_output$num.cores
  registerDoParallel(cl <- makeCluster(num.cores))
  
  
  ##-- Unique mean expressions and dispersion
  mu_star_1_J_output <- normHDP_post_output$mu_star_1_J_output
  phi_star_1_J_output <- normHDP_post_output$phi_star_1_J_output
  mu_star_1_J_new <- normHDP_post_output$mu_star_1_J_new
  phi_star_1_J_new <- normHDP_post_output$phi_star_1_J_new
  tilde_s_unique_new <- normHDP_post_output$tilde_s_unique_new
  mean_X_unique_new <- normHDP_post_output$mean_X_unique_new
  covariance_unique_new <- normHDP_post_output$covariance_unique_new
  unique_count <- normHDP_post_output$unique_count
  acceptance_prob_vec <- normHDP_post_output$acceptance_prob_vec
  
  ##-- Other parameters
  output_index <- normHDP_post_output$output_index
  final_iter <- normHDP_post_output$final_iter
  b_posterior_mean <- normHDP_post_output$b_posterior_mean
  Beta_posterior_mean <- normHDP_post_output$Beta_posterior_mean
  alpha_mu_2 <- normHDP_post_output$alpha_mu_2
  quadratic <- normHDP_post_output$quadratic
  adaptive_prop <- normHDP_post_output$adaptive_prop
  auto.save <- normHDP_post_output$auto.save
  partial.save.name <- normHDP_post_output$partial.save.name
  iter_update <- normHDP_post_output$iter_update
  alpha_phi_2 <- normHDP_post_output$alpha_phi_2
  
  Y <- normHDP_post_output$Y
  G <- normHDP_post_output$G
  J <- normHDP_post_output$J
  D <- normHDP_post_output$D
  C <- normHDP_post_output$D
  Z <- normHDP_post_output$Z
  
  ##-- Update acceptance probability 
  acceptance_prob_vec <- c(acceptance_prob_vec[1:(final_iter-1)],
                           rep(0, number_iter-final_iter))
  
  
  ##-- MCMC continue from iterations - final iteration + 1
  for(iter in (final_iter+1):number_iter){
    
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
    acceptance_prob_list$unique_accept[iter-1] <- unique_count/((iter-1)*J*G)
    
    
    if(update == TRUE){
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      phi_star_1_J_output[[output_index]] <- phi_star_1_J_new
    }else{}
    
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
                      'Y' = Y,
                      'Z' = Z,
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
                  'Y' = Y,
                  'Z' = Z,
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
                  'alpha_phi_2' = alpha_phi_2,
                  'num.cores' = num.cores)
  
  return(my_list)
  
  
 
}