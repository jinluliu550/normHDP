

normHDP_mcmc_fixed_z <- function(Y,
                                 Z,
                                 number_iter,
                                 thinning = 5,
                                 empirical = TRUE,
                                 burn_in = 3000,
                                 quadratic = FALSE,
                                 iter_update = 100,
                                 beta.mean = 0.06,
                                 alpha_mu_2 = NULL,
                                 adaptive_prop = 0.1,
                                 BB_SIZE = TRUE,
                                 num.cores = 4,
                                 baynorm.beta = NULL,
                                 auto.save = FALSE,
                                 auto.save.every.n = 0,
                                 auto.save.name = NULL,
                                 run.on.pc = TRUE){


  ##--------------- dimensions --------------------

  D <- length(Y)
  G <- nrow(Y[[1]])
  C <- unlist(lapply(1:D, function(d) ncol(Y[[d]])))
  Y <- sapply(1:D,
              function(d) as.matrix(Y[[d]]))

  # Set Z to the total number of clusters
  J <- max(unlist(Z))

  ##--------------- bayNorm estimates -------------

  # Beta estimate
  if(is.null(baynorm.beta)){

    baynorm.beta <- lapply(1:D, function(d){

      bayNorm::BetaFun(Data = Y[[d]],
                       MeanBETA = beta.mean)$BETA
    })
  }


  # Apply bayNorm to individual datasets
  baynorm_ind <- lapply(1:D, function(d){

    bayNorm(Data = Y[[d]],
            BETA_vec = baynorm.beta[[d]],
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
  P_J_D_output <- NULL
  P_output <- NULL
  alpha_output <- 0
  alpha_zero_output <- 0
  mu_star_1_J_output <- NULL
  phi_star_1_J_output <- NULL
  Beta_output <- NULL

  #----------------------- Initial values -----------

  ## Set initial values

  if(quadratic){
    b_initial <- c(0,1,0)
  }else{
    b_initial <- c(0,1)
  }


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
  P_J_D_new <- P_J_D_initial
  P_new <- P_initial
  alpha_new <- alpha_initial
  alpha_zero_new <- alpha_zero_initial
  mu_star_1_J_new <- mu_star_1_J_initial
  phi_star_1_J_new <- phi_star_1_J_initial
  Beta_new <- baynorm.beta

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

    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      regression_detail <- list('alpha_phi_2_new' = alpha_phi_2_new,
                                'b_new' = b_new)
    }

    ##------------------------ Dataset specific component probabilities ------------------------------

    dataset_specfic_output <- dataset_specific_mcmc(Z = Z,
                                                    P = P_new,
                                                    alpha = alpha_new)

    P_J_D_new <- dataset_specfic_output

    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      P_J_D_details <- list('P_J_D_new' = P_J_D_new)
    }

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


    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      component_details <- list('P_new' = P_new,
                                'tilde_s_component_new' = tilde_s_component_new,
                                'mean_X_component_new' = mean_X_component_new,
                                'covariance_component_new' = covariance_component_new,
                                'P_count' = P_count)
    }

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

    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      alpha_details <- list('alpha_new' = alpha_new,
                            'mean_X_alpha_new' = mean_X_alpha_new,
                            'M_2_alpha_new' = M_2_alpha_new,
                            'variance_alpha_new' = variance_alpha_new,
                            'alpha_count' = alpha_count)

    }

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


    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      alpha_zero_details <- list('alpha_zero_new' = alpha_zero_new,
                                 'mean_X_alpha_zero_new' = mean_X_alpha_zero_new,
                                 'M_2_alpha_zero_new' = M_2_alpha_zero_new,
                                 'variance_alpha_zero_new' = variance_alpha_zero_new,
                                 'alpha_zero_count' = alpha_zero_count)


    }


    ##---------------------- Unique parameters ----------------------------------

    unique_output_sim <- unique_parameters_mcmc(mu_star_1_J = mu_star_1_J_new,
                                                phi_star_1_J = phi_star_1_J_new,
                                                mean_X_mu_phi = mean_X_unique_new,
                                                tilde_s_mu_phi = tilde_s_unique_new,
                                                Z = Z,
                                                b = b_new,
                                                alpha_phi_2 = alpha_phi_2_new,
                                                Beta = Beta_new,
                                                alpha_mu_2 = alpha_mu_2,
                                                covariance = covariance_unique_new,
                                                Y = Y,
                                                iter_num = iter,
                                                quadratic = quadratic,
                                                adaptive_prop = adaptive_prop,
                                                num.cores = num.cores,
                                                run.on.pc = run.on.pc)

    mu_star_1_J_new <- unique_output_sim$mu_star_1_J_new
    phi_star_1_J_new <- unique_output_sim$phi_star_1_J_new
    tilde_s_unique_new <- unique_output_sim$tilde_s_mu_phi_new
    mean_X_unique_new <- unique_output_sim$mean_X_mu_phi_new
    covariance_unique_new <- unique_output_sim$covariance_new
    unique_count <- unique_count + unique_output_sim$accept_count
    acceptance_prob_list$unique_accept[iter-1] <- unique_count/((iter-1)*J*G)



    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      mu_phi_details <- list('mu_star_1_J_new' = mu_star_1_J_new,
                             'phi_star_1_J_new' = phi_star_1_J_new,
                             'tilde_s_unique_new' = tilde_s_unique_new,
                             'mean_X_unique_new' = mean_X_unique_new,
                             'covariance_unique_new' = covariance_unique_new,
                             'unique_count' = unique_count)


    }

    ##-------------------- Capture efficiency ----------------------------------

    capture_output_sim <- capture_efficiencies_mcmc(Beta = Beta_new,
                                                    Y = Y,
                                                    Z = Z,
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

    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){

      beta_details <- list('Beta_new' = Beta_new,
                           'mean_X_capture_new' = mean_X_capture_new,
                           'M_2_capture_new' = M_2_capture_new,
                           'variance_capture_new' = variance_capture_new,
                           'Beta_count' = Beta_count)



    }

    #-------------------------- Step 5: Update simulated values ------------------------

    if(update == TRUE){
      alpha_phi_2_output[output_index] <- alpha_phi_2_new
      b_output[[output_index]] <- as.vector(b_new)
      P_J_D_output[[output_index]] <- P_J_D_new
      P_output[[output_index]] <- P_new
      alpha_output[output_index] <- alpha_new
      alpha_zero_output[output_index] <- alpha_zero_new
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      phi_star_1_J_output[[output_index]] <- phi_star_1_J_new
      Beta_output[[output_index]] <- Beta_new
    }


    if(isTRUE(auto.save) & (iter %% auto.save.every.n == 0)){


      my_list <- list('b_output' = b_output,
                      'alpha_phi2_output' = alpha_phi_2_output,
                      'P_J_D_output' = P_J_D_output,
                      'P_output' = P_output,
                      'alpha_output' = alpha_output,
                      'alpha_zero_output' = alpha_zero_output,
                      'mu_star_1_J_output' = mu_star_1_J_output,
                      'phi_star_1_J_output' = phi_star_1_J_output,
                      'Beta_output' = Beta_output,
                      'acceptance_prob_list' = acceptance_prob_list,

                      'J' = J,
                      'D' = D,
                      'C' = C,
                      'G' = G,
                      'Z' = Z,

                      'regression_detail' = regression_detail,
                      'P_J_D_details' = P_J_D_details,
                      'component_details' = component_details,
                      'alpha_details' = alpha_details,
                      'alpha_zero_details' = alpha_zero_details,
                      'mu_phi_details' = mu_phi_details,
                      'beta_details' = beta_details)

      save(my_list, file = paste0(auto.save.name, '.RData'))


    }


  }

  ## Return the list
  my_list <- list('b_output' = b_output,
                  'alpha_phi2_output' = alpha_phi_2_output,
                  'P_J_D_output' = P_J_D_output,
                  'P_output' = P_output,
                  'alpha_output' = alpha_output,
                  'alpha_zero_output' = alpha_zero_output,
                  'mu_star_1_J_output' = mu_star_1_J_output,
                  'phi_star_1_J_output' = phi_star_1_J_output,
                  'Beta_output' = Beta_output,
                  'acceptance_prob_list' = acceptance_prob_list,

                  'J' = J,
                  'D' = D,
                  'C' = C,
                  'G' = G,
                  'Z' = Z)


  return(my_list)

}

