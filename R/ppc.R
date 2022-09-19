#' Posterior predictive check with multiple datasets
#'
#' @param normHDP_output Output from function normHDP_mcmc.
#' @param Y Input dataset.
#' @param number_rep Number of replicated datasets to produce.
#' @return A list of 2 items: statistics of replicated datasets (rep_Y_statistics); and statistics of the
#' observed dataset (Y_statistics).
#' @export
ppc <- function(normHDP_output,
                Y,
                number_rep){

  ##-- Input data
  D <- normHDP_output$D
  G <- normHDP_output$G
  C <- normHDP_output$C
  J <- normHDP_output$J

  ##-- Trace from MCMC
  mu_star_1_J_output <- normHDP_output$mu_star_1_J_output
  alpha_phi_2_output <- normHDP_output$alpha_phi2_output
  Beta_output <- normHDP_output$Beta_output
  Z_output <- normHDP_output$Z_output
  b_output <- normHDP_output$b_output
  phi_star_1_J_output <- normHDP_output$phi_star_1_J_output
  P_J_D_output <- normHDP_output$P_J_D_output

  ##-- Total number of samples
  sample.num <- length(mu_star_1_J_output)

  ##-- Generate some random index of length number_rep
  index <- sample(1:sample.num, number_rep, replace = FALSE)

  ##-- Store replicated information
  rep_Y_statistics <- NULL

  ##-- Allocation probability
  allo.prob <- normHDP_output$allocation_prob_output


  ##-- Find cells with high mean allocation probability
  high.allo.index <- NULL

  for(d in 1:D){

    list.result <- lapply(1:number_rep, function(i) allo.prob[[i]][[d]])
    df.result <- do.call(rbind, list.result)

    ##-- Mean allocation probability of each cell
    allo.prob.mean <- colSums(df.result)/number_rep

    ##-- Set a threshold
    threshold <- 0.8

    ##-- Find cells with mean allocation probability greater than the threshold
    high.allo.index[[d]] <- which(allo.prob.mean > threshold)
  }

  ##-- The new cell number
  C_new <- sapply(1:D, function(d) length(high.allo.index[[d]]))

  ##----------------------- Compute discrepancy measures ----------------------------------

  num.index <- 1

  for(i in index){

    ##-- Set theta
    theta <- list('mu' = mu_star_1_J_output[[i]],
                  'b' = b_output[[i]],
                  'alpha_phi_2' = alpha_phi_2_output[[i]],
                  'Z' = Z_output[[i]],
                  'Beta' = Beta_output[[i]])

    ##-- Truncate Z and Beta
    for(d in 1:D){
      theta$Z[[d]] <- theta$Z[[d]][high.allo.index[[d]]]
      theta$Beta[[d]] <- theta$Beta[[d]][high.allo.index[[d]]]
    }


    ##-------------------- Step 2: For replicated data -----------------

    ##-- Simulate phi from mu, b and alpha_phi_2
    if(length(theta$b) == 3){

      phi <- matrix(rlnorm(n = J*G,
                           meanlog = theta$b[1]+theta$b[2]*log(theta$mu)+theta$b[3]*log(theta$mu)^2,
                           sdlog = theta$alpha_phi_2),
                    nrow = J,
                    ncol = G)
    }else{

      phi <- matrix(rlnorm(n = J*G,
                           meanlog = theta$b[1]+theta$b[2]*log(theta$mu),
                           sdlog = theta$alpha_phi_2),
                    nrow = J,
                    ncol = G)
    }

    # Replicated Y
    Y_rep <- lapply(1:D, function(d){

      # Expectation
      Expectations <- t(apply(theta$mu[theta$Z[[d]],1:G],
                              2,
                              function(x) x*theta$Beta[[d]])
      )

      # Size
      size <- t(phi[theta$Z[[d]],1:G])

      matrix(rnbinom(n = G*C_new[d],
                     mu = as.vector(Expectations),
                     size = as.vector(size)),

             nrow = G)
    })

    # Statistics for replicated Y
    Y_rep_statistics_t <- lapply(1:D, function(d){

      data.frame(dataset = d,
                 gene = 1:G,
                 mean.log.shifted.counts = apply(Y_rep[[d]], 1, function(x) mean(log(x+1))),
                 sd.log.shifted.counts = apply(Y_rep[[d]], 1, function(x) sd(log(x+1))),
                 log.mean.counts = apply(Y_rep[[d]], 1, function(x) log(mean(x))),
                 dropout.probability = apply(Y_rep[[d]], 1, function(x) length(which(x == 0))/C_new[d]))
    })

    rep_Y_statistics[[num.index]] <- do.call(rbind, Y_rep_statistics_t)
    rep_Y_statistics[[num.index]]$t <- num.index

    ##-- Update index
    num.index <- num.index + 1
  }

  rep_Y_statistics <- do.call(rbind, rep_Y_statistics)

  ##-- For observed Y
  Y_statistics <- lapply(1:D, function(d){

    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,
                           mean.log.shifted.counts = apply(Y[[d]][, high.allo.index[[d]]],
                                                           1,
                                                           function(x) mean(log(x+1))),

                           sd.log.shifted.counts = apply(Y[[d]][, high.allo.index[[d]]],
                                                         1,
                                                         function(x) sd(log(x+1))),

                           log.mean.counts = apply(Y[[d]][, high.allo.index[[d]]],
                                                   1,
                                                   function(x) log(mean(x))),

                           dropout.probability = apply(Y[[d]][, high.allo.index[[d]]],
                                                       1,
                                                       function(x) length(which(x == 0))/C_new[d]))
  })

  Y_statistics <- do.call(rbind, Y_statistics)
  Y_statistics$t <- 0

  ##-- Return final output
  return(list('rep_Y_statistics' = rep_Y_statistics,
              'Y_statistics' = Y_statistics))
}
