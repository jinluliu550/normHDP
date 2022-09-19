
#' Simulation of regression hyper-parameters
#'
#' @param mu_star_1_J mean expression matrix.
#' @param phi_star_1_J dispersion matrix.
#' @param m_b Fixed parameter in prior to model slope between log mean expression and dispersion.
#' @param v_1 Fixed parameter in prior, shape parameter of the prior for alpha_phi_2.
#' @param v_2 Fixed parameter in prior, scale parameter of the prior for alpha_phi_2.
#' @param quadratic Relationship between mean expressions and dispersion on the log-scale. Default setting is to assume quadratic is FALSE.
#' @return A list of 2 items: estimated variance (alpha_phi_2) and the regression parameter (b).
#' @export
mean_dispersion_mcmc <- function(mu_star_1_J,
                                 phi_star_1_J,
                                 m_b,
                                 v_1,
                                 v_2,
                                 quadratic=FALSE){

  J <- nrow(mu_star_1_J)
  G <- ncol(mu_star_1_J)

  # Setting initial values
  parameter0 <- diag(x=1, nrow=length(m_b), ncol=length(m_b))
  parameter1 <- matrix(m_b,ncol=1)
  parameter2 <- 0

  # For each j, update the parameters
  if(isFALSE(quadratic)){

    ## non-quadratic
    for(j in 1:J) {
      mu_tilde_j <- matrix(c(rep(1,G),log(mu_star_1_J[j,])), ncol=2, nrow=G)
      parameter0 <- parameter0 + t(mu_tilde_j)%*%(mu_tilde_j)
      parameter1 <- parameter1 + t(mu_tilde_j)%*%(matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
      parameter2 <- parameter2 + (matrix(log(phi_star_1_J[j,]), nrow=1, ncol=G))%*%
        (matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
    }
  }else{

    ## quadratic
    for(j in 1:J){
      mu_tilde_j <- matrix(c(rep(1,G),log(mu_star_1_J[j,]),(log(mu_star_1_J[j,]))^2), ncol=3, nrow=G)
      parameter0 <- parameter0 + t(mu_tilde_j)%*%(mu_tilde_j)
      parameter1 <- parameter1 + t(mu_tilde_j)%*%(matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
      parameter2 <- parameter2 + (matrix(log(phi_star_1_J[j,]), nrow=1, ncol=G))%*%
        (matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
    }
  }

  # Setting parameters
  V_tilde_b <- solve(parameter0)
  m_tilde_b <- (V_tilde_b)%*%(parameter1)
  v_tilde_1 <- v_1 + J*G/2
  v_tilde_2 <- v_2 + 1/2*(parameter2-t(m_tilde_b)%*%solve(V_tilde_b)%*%(m_tilde_b) + sum(m_b^2))

  # Simulate alpha_phi_2
  alpha_phi_2 <- rinvgamma(n = 1,
                           alpha = v_tilde_1,
                           beta = v_tilde_2) # rinvgamma comes from the extraDistr package

  # Simulate b
  b <- rmvnorm(n = 1,
               mean = m_tilde_b,
               sigma = alpha_phi_2*(V_tilde_b))

  # Returning alpha_phi_2 and b
  return(list('alpha_phi_2' = alpha_phi_2,
              'b' = b))
}


#' Simulation of allocations
#'
#' @param P_J_D Dataset-specific component probability.
#' @param mu_star_1_J Mean expression matrix.
#' @param phi_star_1_J Dispersion matrix.
#' @param Y Input dataset.
#' @param Beta Capture efficiencies.
#' @param iter_num The iteration number of the current run.
#' @return A list of 2 items: The allocation of all cells (Z, a list of length D, each with a vector of length C[d]) and the allocation probability of
#' each cell to the assigned cluster (allocation.prob, same structure as Z).
#' @export
allocation_variables_mcmc <- function(P_J_D,
                                      mu_star_1_J,
                                      phi_star_1_J,
                                      Y,
                                      Beta,
                                      iter_num){

  # Create an empty list
  D <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(P_J_D)
  G <- ncol(mu_star_1_J)

  ## Initialize Z
  Z <- NULL
  allocation.prob <- NULL

  for(d in 1:D){

    ##-- Put Y into an array
    Y_array <- array(0, c(G, C[d], J))
    for(j in 1:J) Y_array[, , j] <- as.matrix(Y[[d]])

    ##-- Put phi into an array
    phi_array <- array(0, c(G, C[d], J))
    for(c in 1:C[d]) phi_array[, c, ] <- as.matrix(t(phi_star_1_J))

    ##-- Put mu into an array
    mu_array <- array(0, c(G, C[d], J))
    for(c in 1:C[d]) mu_array[, c, ] <- as.matrix(t(mu_star_1_J))

    ##-- Put Beta into an array
    Beta_array <- array(0, c(G, C[d], J))
    for(c in 1:C[d])Beta_array[, c, ] <- matrix(Beta[[d]][c], nrow = G, ncol = J)

    ##-- Compute mu x Beta
    muBeta_array <- array(as.vector(mu_array)*as.vector(Beta_array),
                          c(G, C[d], J))

    rm(mu_array, Beta_array)

    ##-- Log probability for each gene, component and cell
    output.log <- array(dnbinom(as.vector(Y_array),
                                mu = as.vector(muBeta_array),
                                size = as.vector(phi_array),
                                log = TRUE),
                        c(G, C[d], J))

    ##-- Log probability sum over all genes + prior
    LP.matrix <- apply(output.log, c(2,3), sum)
    LP.matrix <- t(apply(LP.matrix, 1, function(x) x + as.vector(log(P_J_D[,d]))))

    ##-- Find the maximum of LP.matrix for each gene (row)
    nc <- -apply(LP.matrix, 1, max)

    ##-- Lp + nc
    LP_plus_nc <- LP.matrix + matrix(nc, nrow = C[d], ncol = J)

    ##-- Allocation probability for cell c in dataset d
    P <- t(apply(LP_plus_nc, 1, function(x) exp(x)/sum(exp(x))))

    ##-- For each row, obtain a sample from posterior
    Z[[d]] <- apply(P, 1, function(x) rcat(n = 1, prob = x))

    ##-- Allocation probability
    allocation.prob[[d]] <- sapply(1:C[d], function(c) P[c,Z[[d]][c]])
  }

  ##-- Return Z and allocation.prob
  return(list('Z' = Z,
              'allocation.prob' = allocation.prob))
}


#' Simulation of dataset-specific component probabilities
#'
#' @param Z Allocation of all cells.
#' @param P Component probabilities.
#' @param alpha Hyper-parameter to indicate variation of dataset-specific component probabilities between datasets.
#' @return Dataset specific component probability; a matrix with J rows and D columns.
#' @export
dataset_specific_mcmc <- function(Z,
                                  P,
                                  alpha){

  J <- length(P)
  D <- length(Z)

  loop.result <- lapply(1:D, function(d) {

    parameter <- as.vector(table(factor(Z[[d]], levels = c(1:J)))) + alpha*P
    P_J_d <- rdirichlet(n = 1, alpha = parameter)

    ## If any of the dataset-specific component probability equal to zero, then give it a small value
    P_J_d <- ifelse(P_J_d == 0, 0.001, P_J_d)
    P_J_d <- P_J_d/sum(P_J_d)

    return(P_J_d)
  })

  P_J_D <- matrix(unlist(loop.result), nrow=J, ncol=D)

  # Returning P_J_D
  return(P_J_D)
}


#' Posterior density of the component probabilities
#'
#' @param P Component probabilities.
#' @param P_J_D Dataset specific component probabilities.
#' @param alpha_zero Hyper-parameter for the component probabilities.
#' @param alpha Hyper-parameters for the dataset specific component probabilities.
#' @return The log probability of posterior density of the component probabilities.
#' @export
component_log_prob <- function(P,
                               P_J_D,
                               alpha_zero,
                               alpha){

  J <- nrow(P_J_D)
  D <- ncol(P_J_D)

  lprod <- sum((alpha_zero/J-1)*log(P)) +

    sum(unlist(lapply(1:D, function(d){

      sum(alpha*P*log(P_J_D[,d])-lgamma(alpha*P))
      })))

  # outputs
  return(lprod)
}


#' Simulation of component probabilities
#'
#' @param P Component probabilities.
#' @param P_J_D Dataset specific component probabilities.
#' @param alpha_zero Hyper-parameter for the component probabilities.
#' @param alpha Hyper-parameters for the dataset specific component probabilities.
#' @param covariance Covariance structure of the component probabilities.
#' @param mean_x Mean of component probabilities.
#' @param tilde_s Sum of cross products of the component probabilities.
#' @param iter_num The iteration number of the current run.
#' @param adaptive_prop Additional variance when simulating new values.
#' @return A list of 5 items: the new component probabilities; updated tilde_s by adding the new observation; updated mean_x; updated
#' covariance structure; and an indicator to indicate whether the value simulated in the current iteration is accepted.
#' @export
component_probabilities_mcmc <- function(P,
                                         P_J_D,
                                         alpha_zero,
                                         alpha,
                                         covariance,
                                         mean_x,
                                         tilde_s,
                                         iter_num,
                                         adaptive_prop = 0.01){

  J <- nrow(P_J_D)
  D <- ncol(P_J_D)

  # Define the inputs
  P_old <- P
  covariance_old <- covariance
  mean_X_d_old <- mean_x
  tilde_s_old <- tilde_s

  X_d_old <- log(P_old[1:J-1]/P_old[J]) # Length = J-1

  # Specify the iteration number
  n <- iter_num

  # Adaptive step - g0 = 100
  if(n <= 100){
    X_new <- rmvnorm(n=1, mean = X_d_old, sigma = diag(x=1, nrow=J-1, ncol=J-1))

  }else{
    X_new <- rmvnorm(n=1, mean = X_d_old, sigma = 2.4^2/(J-1)*(covariance_old + adaptive_prop*diag(1, nrow = J-1, ncol = J-1)))

  }

  # Compute P_new (Length = J) from X_new
  P_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))

  # Compute acceptance probability
  log_acceptance <- component_log_prob(P_new, P_J_D, alpha_zero, alpha) -
    component_log_prob(P_old, P_J_D, alpha_zero, alpha) +
    sum(log(P_new)-log(P_old))


  outcome <- rbinom(n = 1, size = 1, prob=min(1, exp(log_acceptance)))

  if(is.na(outcome) == TRUE | outcome == 0){
    X_new <- X_d_old
    P_new <- P_old
    accept <- 0

  }else{
    accept <- 1
  }


  # Update covariance, mean_x and tilde_s
  tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
  mean_x_new <- mean_X_d_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
  covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new

  # output
  return(list('P_new' = P_new,
              'tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'accept' = accept))
}


#' posterior density of alpha
#'
#' @param P Component probabilities.
#' @param P_J_D Dataset specific component probabilities.
#' @param alpha Hyper-parameters for the dataset specific component probabilities.
#' @return posterior density of alpha on the log-scale.
#' @export
alpha_log_prob <- function(P_J_D,
                           P,
                           alpha){

  D <- ncol(P_J_D)

  lprod <- -alpha + D*lgamma(alpha) +

    sum(unlist(lapply(1:D, function(d) {
      sum(alpha*P*log(P_J_D[,d])-lgamma(alpha*P))
    })))

  # output
  return(lprod)
}


#' Simulation of alpha
#'
#' @param P Component probabilities.
#' @param P_J_D Dataset specific component probabilities.
#' @param alpha Hyper-parameters for the dataset specific component probabilities.
#' @param X_mean Mean of alpha.
#' @param M_2 Sum of square differences between observation and mean.
#' @param variance Variance of alpha.
#' @param iter_num The iteration number of the current run.
#' @param adaptive_prop Additional variance when simulating new values.
#' @return A list of 5 items: the new alpha; updated X_mean; updated M_2; updated variance and an indicator to indicate
#' whether the new simulation is accepted.
#' @export
alpha_mcmc <- function(P_J_D,
                       P,
                       alpha,
                       X_mean,
                       M_2,
                       variance,
                       iter_num,
                       adaptive_prop = 0.01){

  D <- ncol(P_J_D)

  # Defining the inputs
  alpha_old <- alpha
  X_d_old <- log(alpha_old)
  X_mean_old <- X_mean
  M_2_old <- M_2
  variance_old <- variance

  # Defining the dimensions
  n <- iter_num

  # Apply AMH based on the iterative number of the current iteration
  # to simulated new value of X
  if(n <= 100){
    X_new <- rnorm(n = 1, mean = X_d_old, sd = 1)
  }else{
    X_new <- rnorm(n = 1, mean = X_d_old, sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }

  # Transform the new value of X back to new value of alpha, namely alpha_new
  alpha_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- alpha_log_prob(P_J_D, P, alpha = alpha_new) -
    alpha_log_prob(P_J_D, P, alpha = alpha_old) +
    log(alpha_new) - log(alpha_old)

  # Update X_alpha
  outcome <- rbinom(n = 1, size = 1, prob = min(1, exp(log_acceptance)))

  if(is.na(outcome) == TRUE | outcome == 0){
    X_new <- X_d_old
    alpha_new <- alpha_old
    accept <- 0
  }else{
    accept <- 1
  }

  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  variance_new <- 1/(n-1)*M_2_new

  # output
  return(list('alpha_new' = alpha_new,
              'X_mean_new' = X_mean_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new,
              'accept' = accept))
}


#' Posterior density of alpha_zero
#'
#' @param P Component probabilities.
#' @param alpha_zero Hyper-parameter for the component probabilities.
#' @return posterior density of alpha_zero on the log-scale.
#' @export
alpha_zero_log_prob <- function(P,
                                alpha_zero){

  # dimension
  J <- length(P)

  # log probability
  lprob <- -alpha_zero + lgamma(alpha_zero) - J*lgamma(alpha_zero/J) + sum(alpha_zero/J*log(P))

  # output
  return(lprob)
}


#' Simulation of alpha_zero
#'
#' @param P Component probabilities.
#' @param alpha_zero Hyper-parameter for the component probabilities.
#' @param X_mean Mean of alpha_zero.
#' @param M_2 Sum of square differences between observation and mean.
#' @param variance Variance of alpha_zero.
#' @param iter_num The iteration number of the current run.
#' @param adaptive_prop Extra variation.
#' @return A list of 5 items: the new alpha_zero; updated mean; updated M_2; updated variance and an indicator
#' to indicate whether the new simulation is accepted.
#' @export
alpha_zero_mcmc <- function(P,
                            alpha_zero,
                            X_mean,
                            M_2,
                            variance,
                            iter_num,
                            adaptive_prop = 0.01){

  # dimension
  J <- length(P)
  n <- iter_num

  # previous values
  alpha_zero_old <- alpha_zero
  X_d_old <- log(alpha_zero_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean

  # adaptive MH
  if(n <= 100){
    X_new <- rnorm(n = 1, mean = X_d_old, sd = 1)
  }else{
    X_new <- rnorm(n = 1,mean = X_d_old, sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }

  # Obtain the new simulated value for alpha_zero
  alpha_zero_new <- exp(X_new)

  # Compute acceptance probability
  log_acceptance <- alpha_zero_log_prob(P,alpha_zero = alpha_zero_new) -
    alpha_zero_log_prob(P, alpha_zero = alpha_zero_old) +
    log(alpha_zero_new) - log(alpha_zero_old)

  outcome <- rbinom(n = 1, size = 1, prob = min(1, exp(log_acceptance)))

  if(is.na(outcome) == TRUE | outcome == 0){
    X_new <- X_d_old
    alpha_zero_new <- alpha_zero_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update adaptive MH parameters
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  variance_new <- 1/(n-1)*M_2_new

  # output
  return(list('alpha_zero_new' = alpha_zero_new,
              'X_mean_new' = X_mean_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new,
              'accept' = accept))
}

#' Posterior density of mean expressions and dispersion
#'
#' @param mu_star Mean expression for component j and gene g.
#' @param phi_star Dispersion for component h and gene g.
#' @param Z Allocations.
#' @param b Relationship between mean expressions and dispersions.
#' @param alpha_phi_2 Variance for dispersion given mean expressions.
#' @param Y Dataset.
#' @param Beta Capture efficiency.
#' @param j Index of the component.
#' @param g Index of the gene.
#' @param alpha_mu_2 Fixed parameter to control variation of mean expressions.
#' @param quadratic Relationship between mean expressions and dispersion on the log-scale. Default setting is to assume quadratic is FALSE.
#' @return probability density of mean expressions and dispersions for component j and gene g on the log-scale.
#' @export
unique_parameters_log_prob <- function(mu_star,
                                       phi_star,
                                       Z,
                                       b,
                                       alpha_phi_2,
                                       Y,
                                       Beta,
                                       j,
                                       g,
                                       alpha_mu_2,
                                       quadratic=FALSE){

  D <- length(Y)

  # Define B_c_d_j to be the vector of Beta_c_d, such that Z_c_d = j
  # The algorithm select the cells separately for each dataset
  B_c_d_j <- unlist(lapply(1:D, function(d){
    Beta[[d]][which(Z[[d]]==j)]
  }))

  # Define Y_c_g_d_j to be the vector of Y_c_g_d, such that Z_c_d = j and g given
  # The algorithm select the cells separately for the 1st and 2nd dataset
  Y_c_g_d_j <- unlist(lapply(1:D, function(d){
    Y[[d]][g,which(Z[[d]]==j)]
  }))

  if(quadratic==FALSE){
    lprod1 <- -log(mu_star*phi_star) - 1/(2*alpha_mu_2)*(log(mu_star))^2 -
      (log(phi_star)-(b[1]+b[2]*log(mu_star)))^2/(2*alpha_phi_2)
  }else{
    lprod1 <- -log(mu_star*phi_star) - 1/(2*alpha_mu_2)*(log(mu_star))^2 -
      (log(phi_star)-(b[1]+b[2]*log(mu_star)+b[3]*(log(mu_star))^2))^2/(2*alpha_phi_2)
  }



  if(length(Y_c_g_d_j) != 0){

    # if cluster j is non-empty
    lprod2 <- sum(phi_star*log(phi_star/(mu_star*B_c_d_j+phi_star)) + lgamma(Y_c_g_d_j + phi_star) - lgamma(phi_star) - lgamma(Y_c_g_d_j + 1) +
                    Y_c_g_d_j*log(mu_star/(mu_star*B_c_d_j+phi_star)))
    lprod <- lprod1 + lprod2

  }else{

    # if cluster j is empty
    lprod <- lprod1
  }

  # output
  return(lprod)
}


#' Simulation of mean expressions and dispersions
#'
#' @param mu_star_1_J Mean expression matrix.
#' @param phi_star_1_J Dispersion matrix.
#' @param Z Allocations.
#' @param b Relationship between mean expressions and dispersions.
#' @param alpha_phi_2 Variance for dispersion given mean expressions.
#' @param Y Dataset.
#' @param Beta Capture efficiency.
#' @param alpha_mu_2 Fixed parameter to control variation of mean expressions.
#' @param mean_X_mu_phi Mean of mean expressions and dispersion.
#' @param tilde_s_mu_phi sum of cross products of mean expressions and dispersions.
#' @param covariance Variance-covariance for mean expressions and dispersions.
#' @param quadratic Relationship between mean expressions and dispersion on the log-scale. Default setting is to assume quadratic is FALSE.
#' @param iter_num Current iteration.
#' @param adaptive_prop Extra variation.
#' @return A list of 5 items: new mean expression matrix; new dispersion matrix; total number of accepted values; update tilde_s, mean and
#' covariance structure of mean expressions and dispersion.
#' @export
unique_parameters_mcmc <- function(mu_star_1_J,
                                   phi_star_1_J,
                                   mean_X_mu_phi,
                                   tilde_s_mu_phi,
                                   Z,
                                   b,
                                   alpha_phi_2,
                                   Beta,
                                   alpha_mu_2,
                                   covariance,
                                   iter_num,
                                   quadratic=FALSE,
                                   Y,
                                   adaptive_prop = 0.01){

  # dimensions
  J <- nrow(mu_star_1_J)
  G <- ncol(mu_star_1_J)
  n <- iter_num

  # previous values
  mu_star_1_J_new <- matrix(0, nrow = J, ncol = G)
  phi_star_1_J_new <- matrix(0, nrow = J, ncol = G)

  mean_X_mu_phi_old <- mean_X_mu_phi
  covariance_old <- covariance
  tilde_s_mu_phi_old <- tilde_s_mu_phi

  # Prepare for outputs
  covariance_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  tilde_s_mu_phi_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  mean_X_mu_phi_new <- rep(list(rep(list(matrix(0,nrow=1,ncol=2)),G)),J)

  accept_count_tot <- 0

  for(j in 1:J){

    loop.result <- foreach(g = 1:G,
                           .packages = c('mvtnorm'),
                           .export = c('unique_parameters_log_prob')) %dopar% {

      mu_star_old <- mu_star_1_J[j,g]
      phi_star_old <- phi_star_1_J[j,g]
      X_mu_phi_star_old <- log(c(mu_star_old, phi_star_old))

      # Simulate new value of X, based on the previous value and covariance structure
      if(any(unlist(Z)==j)){

        # Cluster j is occupied
        if(n <= 100){
          X_mu_phi_star_new <- rmvnorm(n = 1,
                                       mean = X_mu_phi_star_old,
                                       sigma = diag(1,nrow = 2, ncol = 2))
        }else{
          X_mu_phi_star_new <- rmvnorm(n = 1,
                                       mean = X_mu_phi_star_old,
                                       sigma = (2.4^2)/2*(covariance_old[[j]][[g]] + adaptive_prop*diag(1,nrow = 2, ncol = 2)))
        }
      }else{

        # Cluster j is non-occupied
        mu_star_new <- rlnorm(n = 1,
                              meanlog = 0,
                              sdlog = sqrt(alpha_mu_2))

        phi_star_new <- rlnorm(n = 1,
                               meanlog = b[1] + b[2]*log(mu_star_new),
                               sdlog = sqrt(alpha_phi_2))

        X_mu_phi_star_new <- log(c(mu_star_new, phi_star_new))
      }


      # Convert back to mu_new and phi_new and based on these values
      # to compute the acceptance probability
      mu_star_new <- exp(X_mu_phi_star_new[1])
      phi_star_new <- exp(X_mu_phi_star_new[2])

      if(any(unlist(Z)==j)){

        # Compute log acceptance probability for occupied cluster
        acceptance_prob_log <- unique_parameters_log_prob(mu_star = mu_star_new,
                                                          phi_star = phi_star_new,
                                                          Z = Z,
                                                          b = b,
                                                          alpha_phi_2 = alpha_phi_2,
                                                          Y = Y,
                                                          Beta = Beta,
                                                          j = j,
                                                          g = g,
                                                          alpha_mu_2 = alpha_mu_2,
                                                          quadratic = quadratic) -

          unique_parameters_log_prob(mu_star = mu_star_old,
                                     phi_star = phi_star_old,
                                     Z = Z,
                                     b = b,
                                     alpha_phi_2 = alpha_phi_2,
                                     Y = Y,
                                     Beta = Beta,
                                     j = j,
                                     g = g,
                                     alpha_mu_2 = alpha_mu_2,
                                     quadratic = quadratic) -

          log(mu_star_old) -log(phi_star_old) + log(mu_star_new) + log(phi_star_new)

        outcome <- rbinom(n = 1,
                          size = 1,
                          prob = min(1,exp(acceptance_prob_log)))

      }else{

        # For non occupied cluster, we simulate from the prior, hence always accept
        outcome <- 1
      }


      if(is.na(outcome) == TRUE | outcome == 0){
        X_mu_phi_star_new <- X_mu_phi_star_old
        mu_star_new <- mu_star_old
        phi_star_new <- phi_star_old
        accept_count <- 0
      }else{
        accept_count <- 1
      }

      ## Return outputs
      list('mu_star_new' = mu_star_new,
           'phi_star_new' = phi_star_new,
           'X_mu_phi_star_new' = X_mu_phi_star_new,
           'accept_count' = accept_count)
    }

    ##-- Redistribute outcome from loop.result

    for(g in 1:G){

      mu_star_1_J_new[j,g] <- loop.result[[g]]$mu_star_new
      phi_star_1_J_new[j,g] <- loop.result[[g]]$phi_star_new

      X_mu_phi_star_new <- loop.result[[g]]$X_mu_phi_star_new

      tilde_s_mu_phi_new_j_g <- tilde_s_mu_phi_old[[j]][[g]] + matrix(X_mu_phi_star_new,ncol=1) %*%
        matrix(X_mu_phi_star_new,nrow=1)
      mean_X_mu_phi_new_j_g <- mean_X_mu_phi_old[[j]][[g]]*(1-1/n) + 1/n*matrix(X_mu_phi_star_new,nrow=1)
      covariance_new_j_g <- 1/(n-1)*tilde_s_mu_phi_new_j_g - n/(n-1)*t(mean_X_mu_phi_new_j_g)%*%mean_X_mu_phi_new_j_g

      covariance_new[[j]][[g]] <- covariance_new_j_g
      tilde_s_mu_phi_new[[j]][[g]] <- tilde_s_mu_phi_new_j_g
      mean_X_mu_phi_new[[j]][[g]] <- mean_X_mu_phi_new_j_g

      ## Count the total number of accepted simulation
      accept_count_tot <- accept_count_tot + loop.result[[g]]$accept_count
    }
  }

  # Outputs
  return(list('mu_star_1_J_new' = mu_star_1_J_new,
              'phi_star_1_J_new' = phi_star_1_J_new,
              'accept_count' = accept_count_tot,
              'tilde_s_mu_phi_new' = tilde_s_mu_phi_new,
              'mean_X_mu_phi_new' = mean_X_mu_phi_new,
              'covariance_new' = covariance_new))
}


#' Posterior density of capture efficiency
#'
#' @param Beta_d Capture efficiency of cells in dataset d.
#' @param Y Dataset.
#' @param Z Allocations.
#' @param mu_star_1_J Mean expression matrix.
#' @param phi_star_1_J Dispersion matrix.
#' @param d The index of the dataset.
#' @param a_d_beta Fixed shape parameter 1 of the prior for capture efficiency.
#' @param b_d_beta Fixed shape parameter 1 of the prior for capture efficiency.
#' @return Posterior density of capture efficiency for cells in dataset d on the log-scale.
#' @export
capture_efficiencies_log_prob <- function(Beta_d,
                                          Y,
                                          Z,
                                          mu_star_1_J,
                                          phi_star_1_J,
                                          d,
                                          a_d_beta,
                                          b_d_beta){

  ## Allocation of each cell
  j <- Z[[d]]

  # Construct the log-probability
  lprod1 <- (a_d_beta[d]-1)*log(Beta_d) + (b_d_beta[d]-1)*log(1-Beta_d)

  # mu x beta
  mu_star_beta <- apply(mu_star_1_J[j,], 2, function(x) x*Beta_d)

  # Y x beta
  Y_beta <- apply(t(Y[[d]]), 2, function(x) x*Beta_d)

  lprod2 <- rowSums((phi_star_1_J[j,]+t(Y[[d]]))*log(phi_star_1_J[j,] + mu_star_beta) - Y_beta)
  lprod <- lprod1 - lprod2

  # output
  return(lprod)
}


#' Simulation of capture efficiency
#'
#' @param Beta Capture efficiency.
#' @param Y Dataset.
#' @param Z Allocations.
#' @param mu_star_1_J Mean expression matrix.
#' @param phi_star_1_J Dispersion matrix.
#' @param d The index of the dataset.
#' @param a_d_beta Fixed shape parameter 1 of the prior for capture efficiency.
#' @param b_d_beta Fixed shape parameter 1 of the prior for capture efficiency.
#' @param iter_num The iteration of the current run.
#' @param M_2 Sum of square differences between observation and mean for each cell.
#' @param mean_X Mean of capture efficiency of all cells.
#' @param variance Variance of capture efficiency of all cells.
#' @param adpative_prop Extra variation.
#' @return A list of 5 items: The new simulated capture effiency; total number of accpeted simulations; updated mean, M_2 and variance.
#' @export
capture_efficiencies_mcmc <- function(Beta,
                                      Y,
                                      Z,
                                      mu_star_1_J,
                                      phi_star_1_J,
                                      a_d_beta,
                                      b_d_beta,
                                      iter_num,
                                      M_2,
                                      mean_X,
                                      variance,
                                      adaptive_prop = 0.01){

  # dimensions
  D <- length(Z)
  C <- unlist(lapply(Z, length))
  n <- iter_num

  # old variance structure
  M_2_old <- M_2
  mean_X_old <- mean_X
  variance_old <- variance

  # output name
  M_2_new <- NULL
  Beta_new <- NULL
  X_new <- NULL
  mean_X_new <- NULL
  variance_new <- NULL

  accept_count_tot <- 0

  for(d in 1:D){

    Beta_d_old <- Beta[[d]]
    X_d_old <- log(Beta_d_old/(1-Beta_d_old))

    # Simulate X_d_new
    if(n <= 100){
      X_d_new <- rnorm(n = length(X_d_old),
                     mean = X_d_old,
                     sd = 1)
    }else{
      X_d_new <- rnorm(n = length(X_d_old),
                     mean = X_d_old,
                     sd = sqrt(2.4^2*(variance_old[[d]] + adaptive_prop)))
    }

    Beta_d_new <- exp(X_d_new)/(1+exp(X_d_new))
    Beta_d_new <- ifelse(Beta_d_new == 1, 0.99, Beta_d_new)

    # log acceptance probability
    acceptance_prob_log <- capture_efficiencies_log_prob(Beta_d = Beta_d_new,
                                                         Y = Y,
                                                         Z = Z,
                                                         mu_star_1_J = mu_star_1_J,
                                                         phi_star_1_J = phi_star_1_J,
                                                         d = d,
                                                         a_d_beta = a_d_beta,
                                                         b_d_beta = b_d_beta) -

      capture_efficiencies_log_prob(Beta_d = Beta_d_old,
                                    Y = Y,
                                    Z = Z,
                                    mu_star_1_J = mu_star_1_J,
                                    phi_star_1_J = phi_star_1_J,
                                    d = d,
                                    a_d_beta = a_d_beta,
                                    b_d_beta = b_d_beta) +

      log(Beta_d_new) + log(1-Beta_d_new) - log(Beta_d_old) - log(1-Beta_d_old)

    # random bernoulli
    outcome <- rbinom(n = C[d],
                      size = 1,
                      prob = min(1, exp(acceptance_prob_log)))

    accept.num <- 0

    for(c in 1:C[d]){

      if(is.na(outcome[c]) == TRUE | outcome[c] == 0){

        X_d_new[c] <- X_d_old[c]
        Beta_d_new[c] <- Beta_d_old[c]
      }else{

        accept.num <- accept.num + 1
      }
    }

    Beta_new[[d]] <- Beta_d_new
    X_new[[d]] <- X_d_new

    mean_X_new[[d]] <- (1-1/n)*mean_X_old[[d]] + (1/n)*X_new[[d]]
    M_2_new[[d]] <- M_2_old[[d]] + (X_new[[d]]-mean_X_old[[d]])*(X_new[[d]]-mean_X_new[[d]])
    variance_new[[d]] <- 1/(n-1)*M_2_new[[d]]

    accept_count_tot <- accept_count_tot + accept.num
  }

  ## return
  return(list('Beta_new' = Beta_new,
              'accept_count' = accept_count_tot,
              'mean_X_new' = mean_X_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new))
}

