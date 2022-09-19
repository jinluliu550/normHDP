#' MCMC Simulation
#'
#' This function generates posterior estimates of unknown parameters in the normHDP model.
#'
#' @param Y Input dataset.
#' @param J Total number of components.
#' @param number_iter Total number of iterations.
#' @param thinning Default is to set thinning to 5; saves posterior estimates for every fifth iteration.
#' @param burn_in Default is to set to 3000; the first 3000 iterations are ignored.
#' @param quadratic Default is set to FALSE, i.e.: assumption that the relationship between the unique parameters are linear.
#' @param iter_update Controls for text, i.e.: if iter_update is set to 100, then the algorithm will output a message every 100 iterations.
#' @param beta.mean Mean of capture efficiency in the prior.
#' @param alpha_mu_2 Variation of log mean expressions.
#' @param partial.save.name Default is to set to NULL, an option to give a name to the output data if set auto.save = TRUE.
#' @param auto.save An option to save the output every 100 iterations. Default is to set to FALSE.
#' @param adpative_prop Extra variation of the random walk.
#' @param BB_SIZE An option to give to the bayNorm function to determine the method used to estimate dispersion in the prior.
#' @param print_Z An option to output a summary of the allocations after each iteration.
#' @param num.cores Number of cores to be used by the function.
#' @return A list of posterior samples of parameters of interests; acceptance probabilities and variance-covariance
#' structure for parameters which need adaptive Metropolis Hastings; and all the fixed parameters.
#' @export
normHDP_mcmc <- function(Y,
                         J,
                         number_iter,
                         thinning = 5,
                         empirical = TRUE,
                         burn_in = 3000,
                         quadratic=FALSE,
                         iter_update=100,
                         beta.mean = 0.06,
                         alpha_mu_2 = NULL,
                         partial.save.name = NULL,
                         auto.save = FALSE,
                         adaptive_prop = 0.1,
                         BB_SIZE = TRUE,
                         print_Z = FALSE,
                         num.cores = 4){

  ##--------------- dimensions --------------------

  D <- length(Y)
  G <- nrow(Y[[1]])
  C <- unlist(lapply(1:D, function(d) ncol(Y[[d]])))

  ##--------------- bayNorm estimates -------------

  # Beta estimate
  baynorm.beta <- lapply(1:D, function(d){

    bayNorm::BetaFun(Data = Y[[d]],
                     MeanBETA = beta.mean)
  })

  # Apply bayNorm to individual datasets
  baynorm_ind <- lapply(1:D, function(d){

    bayNorm(Data = Y[[d]],
            BETA_vec = baynorm.beta[[d]]$BETA,
            mode_version = TRUE,
            mean_version = FALSE,
            BB_SIZE = BB_SIZE)
  })

  # Apply bayNorm to the combined dataset
  baynorm_tot <- bayNorm(Data = do.call(cbind, Y),
                         BETA_vec = unlist(baynorm.beta),
                         mode_version = TRUE,
                         mean_version = FALSE,
                         BB_SIZE = BB_SIZE)

  # Estimate of mu
  baynorm.mu.estimate <- baynorm_tot$PRIORS$MME_prior[,1]

  # Estimate of phi
  if(isFALSE(BB_SIZE)){
    baynorm.phi.estimate <- baynorm_tot$PRIORS$MME_prior[,2]
  }else{
    baynorm.phi.estimate <- baynorm_tot$PRIORS$MME_SIZE_adjust
  }

  #------------------------ Prepare for outputs -----------------

  b_output <- NULL
  alpha_phi_2_output <- 0
  Z_output <- NULL
  P_J_D_output <- NULL
  P_output <- NULL
  alpha_output <- 0
  alpha_zero_output <- 0
  mu_star_1_J_output <- NULL
  phi_star_1_J_output <- NULL
  Beta_output <- NULL
  allocation_prob_output <- NULL

  #----------------------- Initial values -----------

  ## Set initial values

  if(quadratic){
    b_initial <- c(0,1,0)
  }else{
    b_initial <- c(0,1)
  }

  ## Random allocation
  Z_initial <- lapply(1:D, function(d) sample(1:J, size = C[d], replace = TRUE))

  P_J_D_initial <- matrix(1/J, nrow = J, ncol = D)
  P_initial <- rep(1/J, J)
  alpha_initial <- 1
  alpha_zero_initial <- 1

  #---------------------------- Hyper-parameters -----------------------

  ##----------- alpha_mu_2
  if(is.null(alpha_mu_2)){
    log.mu.finite.index <- which(is.finite(log(baynorm.mu.estimate)))
    alpha_mu_2 <- 2*log(mean(baynorm.mu.estimate[log.mu.finite.index]))
  }

  ##----------- a_d_beta, b_d_beta
  loop.result <- lapply(1:D, function(d){

    ##-- Initial mean and variance
    baynorm_mean_capeff <- beta.mean
    baynorm_var_capeff <- 0.5

    ##-- Initial a_beta and b_beta
    a_beta <- ((1-baynorm_mean_capeff)/baynorm_var_capeff - 1/baynorm_mean_capeff)*baynorm_mean_capeff^2
    b_beta <- a_beta*(1/baynorm_mean_capeff - 1)

    while((baynorm_var_capeff >= baynorm_mean_capeff*(1-baynorm_mean_capeff)) | (a_beta < 1) | (b_beta < 1)){

      ##-- While the above condition is true, half the variance
      baynorm_var_capeff <- baynorm_var_capeff/2

      ##-- Compute the new a_beta and b_beta
      a_beta <- ((1-baynorm_mean_capeff)/baynorm_var_capeff - 1/baynorm_mean_capeff)*baynorm_mean_capeff^2
      b_beta <- a_beta*(1/baynorm_mean_capeff - 1)

    }

    return(c('a_beta' = a_beta,
             'b_beta' = b_beta))
  })

  a_d_beta <- unlist(lapply(loop.result, function(i) i['a_beta']))
  b_d_beta <- unlist(lapply(loop.result, function(i) i['b_beta']))

  ##------------ v_1, v_2 and m_b

  # Fit a linear model - regression coefficients
  baynorm.phi.estimate[is.infinite(log(baynorm.mu.estimate))] <- NA
  baynorm.mu.estimate[is.infinite(log(baynorm.mu.estimate))] <- NA

  x <- log(baynorm.mu.estimate)
  y <- log(baynorm.phi.estimate)

  if(isFALSE(quadratic)){
    lm.1 <- lm(y ~ x)
  }else{
    lm.1 <- lm(y ~ x+I(x^2))
  }

  # Summary statistics
  rse.lm.1.squared <- deviance(lm.1)/df.residual(lm.1)
  variance <- 5
  v_1 <- rse.lm.1.squared^2/variance + 2
  v_2 <- (v_1 - 1)*rse.lm.1.squared

  if(isFALSE(empirical)){
    v_1 <- 2
    v_2 <- 1
  }

  ##-- alpha_phi2 initial value
  alpha_phi_2_initial <- rse.lm.1.squared

  ##-- m_b
  if(isTRUE(empirical)){
    m_b = as.numeric(coef(lm.1))
  }else{

    if(isTRUE(quadratic)){
      m_b <- c(-1,2,0)
    }else{
      m_b <- c(-1,2)
    }
  }


  ##-- Initialize mean and dispersion
  baynorm.mu.estimate <- ifelse(baynorm.mu.estimate==0, 0.01, baynorm.mu.estimate)
  baynorm.phi.estimate <- ifelse(baynorm.phi.estimate==0, 0.01, baynorm.phi.estimate)

  mu_star_1_J_initial <- t(matrix(baynorm.mu.estimate, nrow=G, ncol=J))
  phi_star_1_J_initial <- t(matrix(baynorm.phi.estimate, nrow=G, ncol=J))


  #--------------------------------------------------------------------------------

  # Acceptance probability
  acceptance_prob_list <- data.frame(P_accept = rep(0,number_iter-1),
                                     alpha_accept = rep(0,number_iter-1),
                                     alpha_zero_accept = rep(0,number_iter-1),
                                     unique_accept = rep(0,number_iter-1),
                                     Beta_accept = rep(0,number_iter-1))


  #----------------------- Step 2: Set the initial_values as new values ----------------------------
  b_new <- b_initial
  alpha_phi_2_new <- alpha_phi_2_initial
  Z_new <- Z_initial
  P_J_D_new <- P_J_D_initial
  P_new <- P_initial
  alpha_new <- alpha_initial
  alpha_zero_new <- alpha_zero_initial
  mu_star_1_J_new <- mu_star_1_J_initial
  phi_star_1_J_new <- phi_star_1_J_initial
  Beta_new <- lapply(1:D, function(d) baynorm.beta[[d]]$BETA)

  # Total count of acceptance
  P_count <- 0
  alpha_count <- 0
  alpha_zero_count <- 0
  unique_count <- 0
  Beta_count <- 0

  #----------------------- Step 3: Prepare for the covariance update -----------------------------

  # 1) For Component probabilities
  mean_X_component_new <- matrix(log(P_initial[1:J-1]/P_initial[J]), nrow = 1)
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  covariance_component_new <- matrix(0, nrow = J-1, ncol = J-1)

  # 2) For alpha
  mean_X_alpha_new <- log(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0

  # 3) For alpha_zero
  mean_X_alpha_zero_new <- log(alpha_zero_new)
  M_2_alpha_zero_new <- 0
  variance_alpha_zero_new <- 0

  # 4) Unique parameters
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


  # 5) Capture efficiency
  variance_capture_new <- lapply(1:D, function(d) rep(0,C[d]))
  mean_X_capture_new <- lapply(1:D, function(d) log(Beta_new[[d]]/(1-Beta_new[[d]])))
  M_2_capture_new <- lapply(1:D, function(d) rep(0,C[d]))

  ##-- Register Cores
  registerDoParallel(cl <- makeCluster(num.cores))

  output_index <- 0

  #----------------------- Step 4: Updates -----------------------------

  # Iteration starts with number 2
  for(iter in 2:number_iter){

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

    #-------------------------- Step 5: Update simulated values ------------------------
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
    }else{}

    ##-- Save out result every 100th iteration if auto.save == TRUE
    if(iter %% 100 == 0 && auto.save == TRUE){
      my_list <- list('b_output' = b_output,
                      'alpha_phi2_output' = alpha_phi_2_output,
                      'Z_output' = Z_output,
                      'Z_new' = Z_new,
                      'allocation_prob_output' = allocation_prob_output,
                      'allocation.prob' = allocation.prob,
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
                      'G' = G,
                      'thinning' = thinning,
                      'empirical' = empirical,
                      'burn_in' = burn_in,
                      'iter_update' = iter_update,
                      'beta.mean' = beta.mean,
                      'partial.save.name' = partial.save.name,
                      'auto.save' = auto.save,
                      'BB_SIZE' = BB_SIZE,
                      'print_Z' = print_Z,
                      'mu_star_1_J_initial' = mu_star_1_J_initial,
                      'phi_star_1_J_initial' = phi_star_1_J_initial,
                      'num.cores' = num.cores)

      save(my_list, file=partial.save.name)
    }
  }

  ## Return the list
  my_list <- list('b_output' = b_output,
                  'alpha_phi2_output' = alpha_phi_2_output,
                  'Z_output' = Z_output,
                  'Z_new' = Z_new,
                  'P_J_D_output' = P_J_D_output,
                  'allocation_prob_output' = allocation_prob_output,
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
                  'G' = G,
                  'thinning' = thinning,
                  'empirical' = empirical,
                  'burn_in' = burn_in,
                  'iter_update' = iter_update,
                  'partial.save.name' = partial.save.name,
                  'auto.save' = auto.save,
                  'BB_SIZE' = BB_SIZE,
                  'print_Z' = print_Z,
                  'mu_star_1_J_initial' = mu_star_1_J_initial,
                  'phi_star_1_J_initial' = phi_star_1_J_initial,
                  'num.cores' = num.cores)

  return(my_list)

}

