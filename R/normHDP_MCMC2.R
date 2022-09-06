
#' Continuation of the normHDP_mcmc function
#' 
#' Function to be used if want to obtain more posterior samples after running normHDP_mcmc
#' 
#' @param normHDP_output Output from the normHDP_mcmc function.
#' @param number_iter Total number of iterations.
#' @return A list of posterior samples of parameters of interests; acceptance probabilities and variance-covariance
#' structure for parameters which need adaptive Metropolis Hastings; and all the fixed parameters.
#' @export
normHDP_mcmc_2 <- function(normHDP_output, 
                           number_iter){
  
  ##-- Mean and variance hyper-parameters
  b_output <- normHDP_output$b_output
  alpha_phi_2_output <- normHDP_output$alpha_phi_2_output
  
  ##-- Acceptance probability
  acceptance_prob_list <- normHDP_output$acceptance_prob_list
  
  ##-- Allocation
  Z_new <- normHDP_output$Z_new
  Z_output <- normHDP_output$Z_output
  allocation_prob_output <- normHDP_output$allocation_prob_output
  
  ##-- Dataset specific component probability
  P_J_D_output <- normHDP_output$P_J_D_output
  P_J_D_new <- normHDP_output$P_J_D_new
  
  ##-- Component probability
  P_new <- normHDP_output$P_new
  tilde_s_component_new <- normHDP_output$tilde_s_component_new
  mean_X_component_new <- normHDP_output$mean_X_component_new
  covariance_component_new <- normHDP_output$covariance_component_new
  P_count <- normHDP_output$P_count
  P_output <- normHDP_output$P_output
  
  ##-- Alpha
  alpha_new <- normHDP_output$alpha_new
  mean_X_alpha_new <- normHDP_output$mean_X_alpha_new
  M_2_alpha_new <- normHDP_output$M_2_alpha_new
  variance_alpha_new <- normHDP_output$variance_alpha_new
  alpha_count <- normHDP_output$alpha_count
  alpha_output <- normHDP_output$alpha_output
  
  ##-- Alpha_zero
  alpha_zero_new <- normHDP_output$alpha_zero_new
  mean_X_alpha_zero_new <- normHDP_output$mean_X_alpha_zero_new
  M_2_alpha_zero_new <- normHDP_output$M_2_alpha_zero_new
  alpha_zero_output <- normHDP_output$alpha_zero_output
  variance_alpha_zero_new <- normHDP_output$variance_alpha_zero_new
  alpha_zero_count <- normHDP_output$alpha_zero_count
  
  ##-- Unique parameters
  mu_star_1_J_new <- normHDP_output$mu_star_1_J_new
  phi_star_1_J_new <- normHDP_output$phi_star_1_J_new
  tilde_s_unique_new <- normHDP_output$tilde_s_unique_new
  mean_X_unique_new <- normHDP_output$mean_X_unique_new
  covariance_unique_new <- normHDP_output$covariance_unique_new
  unique_count <- normHDP_output$unique_count
  mu_star_1_J_output <- normHDP_output$mu_star_1_J_output
  phi_star_1_J_output <- normHDP_output$phi_star_1_J_output
  
  ##-- Capture efficiencies
  Beta_new <- normHDP_output$Beta_new
  mean_X_capture_new <- normHDP_output$mean_X_capture_new
  M_2_capture_new <- normHDP_output$M_2_capture_new
  variance_capture_new <- normHDP_output$variance_capture_new
  Beta_count <- normHDP_output$Beta_count
  Beta_output <- normHDP_output$Beta_output
  
  ##-- Fixed parameters
  alpha_mu_2 <- normHDP_output$alpha_mu_2
  a_d_beta <- normHDP_output$a_d_beta
  b_d_beta <- normHDP_output$b_d_beta
  v_1 <- normHDP_output$v_1
  v_2 <- normHDP_output$v_2
  m_b <- normHDP_output$m_b
  adaptive_prop <- normHDP_output$adaptive_prop
  thinning <- normHDP_output$thinning
  burn_in <- normHDP_output$burn_in
  final_iter <- normHDP_output$final_iter
  iter_update <- normHDP_output$iter_update
  beta.mean <- beta.mean
  J <- normHDP_output$J
  D <- normHDP_output$D
  C <- normHDP_output$C
  G <- normHDP_output$G
  
  ##-- Other parameters
  quadratic <- normHDP_output$quadratic
  output_index <- normHDP_output$output_index
  empirical <- normHDP_output$empirical
  partial.save.name <- normHDP_output$partial.save.name
  auto.save <- normHDP_output$auto.save
  BB_SIZE <- normHDP_output$BB_SIZE
  print_Z <- normHDP_output$print_Z
  
  ##-- Update acceptance probability table
  ##-- Here we need to consider 2 cases:
  ##-- Case 1: previous MCMC did not fully run (i.e.: did not run until the specified final iteration)
  ##-- Case 2: previous MCMC did fully run
  acceptance_prob_list <- rbind(acceptance_prob_list[1:(final_iter-1),],
                                
                                data.frame(P_accept = rep(0, number_iter-final_iter), 
                                           alpha_accept = rep(0, number_iter-final_iter),
                                           alpha_zero_accept = rep(0, number_iter-final_iter),
                                           unique_accept = rep(0, number_iter-final_iter),
                                           Beta_accept = rep(0, number_iter-final_iter)))
  
  ##-- Register cores
  num.cores <- normHDP_output$num.cores
  registerDoParallel(cl <- makeCluster(num.cores))
  
  
  ######################################### Running MCMC ####################################
  
  for(iter in (final_iter+1):number_iter){
    
    # 0) Starting value of the output index = 1
    # If the current iteration is greater than the burn in and divisible by the thinning index
    if(iter >= burn_in & iter%%thinning == 0){
      output_index <- output_index + 1
      update <- TRUE
    }else{
      update <- FALSE
    }
    
    if(iter %% iter_update == 0){
      print(paste("iter:", iter))
    }
    
    ##---------------------------------- Regression parameters ------------------------------------
    
    mean_dispersion_output <- mean_dispersion_mcmc(mu_star_1_J = mu_star_1_J_new, 
                                                   phi_star_1_J = phi_star_1_J_new,
                                                   v_1 = v_1,
                                                   v_2 = v_2,
                                                   m_b = m_b,
                                                   quadratic = quadratic)
    
    alpha_phi_2_new <- mean_dispersion_output$alpha_phi_2
    b_new <- mean_dispersion_output$b
    
    ##----------------------------------- Allocation variables ------------------------------------
    
    
    allocation_output <- allocation_variables_mcmc(P_J_D = P_J_D_new, 
                                                   mu_star_1_J = mu_star_1_J_new,
                                                   phi_star_1_J = phi_star_1_J_new, 
                                                   Y = Y, 
                                                   Beta = Beta_new,
                                                   iter_num = iter)
    
    Z_new <- allocation_output$Z
    allocation.prob <- allocation_output$allocation.prob
    
    ##-- If print_Z is TRUE, we print out a summary table of allocations for every iteration
    if(isTRUE(print_Z)){
      for(d in 1:D) print(table(Z_new[[d]]))
    }
    
    ##------------------------ Dataset specific component probabilities ------------------------------
    
    dataset_specfic_output <- dataset_specific_mcmc(Z = Z_new, 
                                                    P = P_new, 
                                                    alpha = alpha_new)
    
    P_J_D_new <- dataset_specfic_output
    
    ##------------------------ Component probabilities -----------------------------------------------
    
    component_output <- component_probabilities_mcmc(P = P_new, 
                                                     P_J_D = P_J_D_new, 
                                                     alpha_zero = alpha_zero_new,
                                                     alpha = alpha_new, 
                                                     covariance = covariance_component_new,
                                                     mean_x = mean_X_component_new, 
                                                     tilde_s = tilde_s_component_new,
                                                     iter_num = iter,
                                                     adaptive_prop = adaptive_prop)
    
    P_new <- component_output$P_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    P_count <- P_count + component_output$accept
    acceptance_prob_list$P_accept[iter-1] <- P_count/(iter-1)
    
    ##-------------------------- Update alpha ----------------------------------
    
    alpha_output_sim <- alpha_mcmc(P_J_D = P_J_D_new, 
                                   P = P_new, 
                                   alpha = alpha_new,
                                   X_mean = mean_X_alpha_new, 
                                   M_2 = M_2_alpha_new, 
                                   variance = variance_alpha_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)
    
    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept
    acceptance_prob_list$alpha_accept[iter-1] <- alpha_count/(iter-1)
    
    ##----------------------- Update alpha_zero --------------------------------
    
    alpha_zero_output_sim <- alpha_zero_mcmc(P = P_new, 
                                             alpha_zero = alpha_zero_new,
                                             X_mean = mean_X_alpha_zero_new,
                                             M_2 = M_2_alpha_zero_new, 
                                             variance = variance_alpha_zero_new, 
                                             iter_num = iter,
                                             adaptive_prop = adaptive_prop)
    
    alpha_zero_new <- alpha_zero_output_sim$alpha_zero_new
    mean_X_alpha_zero_new <- alpha_zero_output_sim$X_mean_new
    M_2_alpha_zero_new <- alpha_zero_output_sim$M_2_new
    variance_alpha_zero_new <- alpha_zero_output_sim$variance_new
    alpha_zero_count <- alpha_zero_count + alpha_zero_output_sim$accept
    acceptance_prob_list$alpha_zero_accept[iter-1] <- alpha_zero_count/(iter-1)
    
    ##---------------------- Unique parameters ----------------------------------
    
    unique_output_sim <- unique_parameters_mcmc(mu_star_1_J = mu_star_1_J_new, 
                                                phi_star_1_J = phi_star_1_J_new,
                                                mean_X_mu_phi = mean_X_unique_new,
                                                tilde_s_mu_phi = tilde_s_unique_new,
                                                Z = Z_new,
                                                b = b_new, 
                                                alpha_phi_2 = alpha_phi_2_new,
                                                Beta = Beta_new,
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
    
    
    ##-------------------- Capture efficiencies ----------------------------------
    
    capture_output_sim <- capture_efficiencies_mcmc(Beta = Beta_new,
                                                    Y = Y,
                                                    Z = Z_new,
                                                    mu_star_1_J = mu_star_1_J_new,
                                                    phi_star_1_J = phi_star_1_J_new, 
                                                    a_d_beta = a_d_beta,
                                                    b_d_beta = b_d_beta, 
                                                    iter_num = iter, 
                                                    M_2 = M_2_capture_new,
                                                    mean_X = mean_X_capture_new,
                                                    variance = variance_capture_new,
                                                    adaptive_prop = adaptive_prop)
    
    Beta_new <- capture_output_sim$Beta_new
    mean_X_capture_new <- capture_output_sim$mean_X_new
    M_2_capture_new <- capture_output_sim$M_2_new
    variance_capture_new <- capture_output_sim$variance_new
    Beta_count <- Beta_count + capture_output_sim$accept_count
    acceptance_prob_list$Beta_accept[iter-1] <- Beta_count/((iter-1)*(sum(C)))
    
    #-------------------------- Update simulated values ------------------------
    
    if(update == TRUE){
      alpha_phi_2_output[output_index] <- alpha_phi_2_new
      b_output[[output_index]] <- as.vector(b_new)
      Z_output[[output_index]] <- Z_new
      P_J_D_output[[output_index]] <- P_J_D_new
      P_output[[output_index]] <- P_new
      alpha_output[output_index] <- alpha_new
      alpha_zero_output[output_index] <- alpha_zero_new
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      phi_star_1_J_output[[output_index]] <- phi_star_1_J_new
      Beta_output[[output_index]] <- Beta_new
      allocation_prob_output[[output_index]] <- allocation.prob
    }
    
    ##----------- Save out result every 100th iteration if auto.save == TRUE ------
    
    if(iter %% 100 == 0 && auto.save == TRUE){
      my_list <- list('b_output' = b_output,
                      'alpha_phi2_output' = alpha_phi_2_output,
                      'Z_output' = Z_output, 
                      'Z_new' = Z_new,
                      'allocation_prob_output' = allocation_prob_output,
                      'P_J_D_output' = P_J_D_output,
                      'P_output' = P_output, 
                      'alpha_output' = alpha_output,
                      'alpha_zero_output' = alpha_zero_output, 
                      'mu_star_1_J_output' = mu_star_1_J_output,
                      'phi_star_1_J_output' = phi_star_1_J_output, 
                      'Beta_output' = Beta_output,
                      'acceptance_prob_list' = acceptance_prob_list,
                      'P_new' = P_new, 
                      'tilde_s_component_new' = tilde_s_component_new,
                      'mean_X_component_new' = mean_X_component_new,
                      'covariance_component_new' = covariance_component_new,
                      'P_count' = P_count,
                      'alpha_new' = alpha_new,
                      'mean_X_alpha_new' = mean_X_alpha_new,
                      'M_2_alpha_new' = M_2_alpha_new,
                      'variance_alpha_new' = variance_alpha_new,
                      'alpha_count' = alpha_count, 
                      'alpha_zero_new' = alpha_zero_new,
                      'mean_X_alpha_zero_new' = mean_X_alpha_zero_new, 
                      'M_2_alpha_zero_new' = M_2_alpha_zero_new,
                      'variance_alpha_zero_new' = variance_alpha_zero_new,
                      'alpha_zero_count' = alpha_zero_count,
                      'mu_star_1_J_new' = mu_star_1_J_new,
                      'phi_star_1_J_new' = phi_star_1_J_new,
                      'tilde_s_unique_new' = tilde_s_unique_new,
                      'mean_X_unique_new' = mean_X_unique_new, 
                      'covariance_unique_new' = covariance_unique_new,
                      'unique_count' = unique_count,
                      'Beta_new' = Beta_new,
                      'mean_X_capture_new' = mean_X_capture_new,
                      'M_2_capture_new' = M_2_capture_new,
                      'variance_capture_new' = variance_capture_new, 
                      'Beta_count' = Beta_count,
                      'output_index' = output_index, 
                      'P_J_D_new' = P_J_D_new,
                      'alpha_mu_2' = alpha_mu_2,
                      'a_d_beta' = a_d_beta,
                      'b_d_beta' = b_d_beta,
                      'v_1' = v_1,
                      'v_2' = v_2,
                      'm_b' = m_b,
                      'quadratic' = quadratic,
                      'final_iter' = iter,
                      'adaptive_prop' = adaptive_prop,
                      'J' = J,
                      'D' = D,
                      'C' = C,
                      'thinning' = thinning,
                      'empirical' = empirical,
                      'burn_in' = burn_in,
                      'iter_update' = iter_update,
                      'beta.mean' = beta.mean,
                      'partial.save.name' = partial.save.name,
                      'auto.save' = auto.save,
                      'BB_SIZE' = BB_SIZE,
                      'print_Z' = print_Z,
                      'num.cores' = num.cores)
      
      save(my_list, file=partial.save.name)
    }
  }
  
  ##---------------------- Return final list -------------------------------
  
  my_list <- list('b_output' = b_output, 
                  'alpha_phi2_output' = alpha_phi_2_output, 
                  'Z_output' = Z_output,
                  'Z_new' = Z_new,
                  'allocation_prob_output' = allocation_prob_output,
                  'P_J_D_output' = P_J_D_output, 
                  'P_output' = P_output,
                  'alpha_output' = alpha_output,
                  'alpha_zero_output' = alpha_zero_output,
                  'mu_star_1_J_output' = mu_star_1_J_output,
                  'phi_star_1_J_output' = phi_star_1_J_output,
                  'Beta_output' = Beta_output,
                  'acceptance_prob_list' = acceptance_prob_list, 
                  'P_new' = P_new,
                  'tilde_s_component_new' = tilde_s_component_new,
                  'mean_X_component_new' = mean_X_component_new,
                  'covariance_component_new' = covariance_component_new,
                  'P_count' = P_count,
                  'alpha_new' = alpha_new, 
                  'mean_X_alpha_new' = mean_X_alpha_new,
                  'M_2_alpha_new' = M_2_alpha_new,
                  'variance_alpha_new' = variance_alpha_new, 
                  'alpha_count' = alpha_count, 
                  'alpha_zero_new' = alpha_zero_new,
                  'mean_X_alpha_zero_new' = mean_X_alpha_zero_new,
                  'M_2_alpha_zero_new' = M_2_alpha_zero_new, 
                  'variance_alpha_zero_new' = variance_alpha_zero_new,
                  'alpha_zero_count' = alpha_zero_count, 
                  'mu_star_1_J_new' = mu_star_1_J_new, 
                  'phi_star_1_J_new' = phi_star_1_J_new,
                  'tilde_s_unique_new' = tilde_s_unique_new,
                  'mean_X_unique_new' = mean_X_unique_new, 
                  'covariance_unique_new' = covariance_unique_new,
                  'unique_count' = unique_count, 
                  'Beta_new' = Beta_new, 
                  'mean_X_capture_new' = mean_X_capture_new,
                  'M_2_capture_new' = M_2_capture_new, 
                  'variance_capture_new' = variance_capture_new, 
                  'Beta_count' = Beta_count,
                  'output_index' = output_index,
                  'P_J_D_new' = P_J_D_new,
                  'alpha_mu_2' = alpha_mu_2,
                  'a_d_beta' = a_d_beta,
                  'b_d_beta' = b_d_beta,
                  'quadratic' = quadratic,
                  'final_iter' = iter,
                  'adaptive_prop' = adaptive_prop,
                  'J' = J,
                  'D' = D,
                  'C' = C,
                  'thinning' = thinning,
                  'empirical' = empirical,
                  'burn_in' = burn_in,
                  'iter_update' = iter_update,
                  'partial.save.name' = partial.save.name,
                  'auto.save' = auto.save,
                  'BB_SIZE' = BB_SIZE,
                  'print_Z' = print_Z,
                  'num.cores' = num.cores)
  
  return(my_list)
  
}