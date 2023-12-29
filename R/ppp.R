Discrepancy <- function(Y,
                        theta){

  # Dimensions
  D <- length(Y)
  C <- unlist(lapply(Y,
                     function(x) ncol(x)))
  G <- nrow(Y[[1]])

  summary.df <- lapply(1:D,
                       function(d){

                         # D1
                         summary.mean <- log(abs(Y[[d]] - t(theta$mu[theta$Z[[d]],])*matrix(rep(theta$beta[[d]],
                                                                                                each = G),

                                                                                            nrow = G,
                                                                                            ncol = C[d])))


                         summary.mean <- matrixStats::rowMedians(summary.mean)

                         # D2
                         summary.sd <- log(abs(Y[[d]] - t(theta$mu[theta$Z[[d]],])*matrix(rep(theta$beta[[d]],
                                                                                              each = G),

                                                                                          nrow = G,
                                                                                          ncol = C[d])))

                         summary.sd <- matrixStats::rowIQRs(summary.sd)


                         # D3
                         summary.dropout <- log(abs(ifelse(Y[[d]]==0, 1, 0) - dnbinom(x = 0,
                                                                                      mu = t(theta$mu[theta$Z[[d]],])*matrix(rep(theta$beta[[d]],
                                                                                                                                 each = G),

                                                                                                                             nrow = G,
                                                                                                                             ncol = C[d]),
                                                                                      size = t(theta$phi[theta$Z[[d]],])
                         )
                         )
                         )

                         summary.dropout <- rowMeans(summary.dropout)

                         # Output
                         data.frame(dataset = d,
                                    gene = 1:G,
                                    D1 = summary.mean,
                                    D2 = summary.sd,
                                    D3 = summary.dropout)
                       })

  return(do.call(rbind,
                 summary.df))
}

ppp_calculation <- function(normHDP_post_output,
                            Y,
                            num.cores = 1,
                            rep.number,
                            type = c('standard', 'mixed'),
                            run.on.pc = TRUE){



  # Index of the run
  random.index <- sample(1:length(normHDP_post_output$b_output),
                             size = rep.number,
                             replace = TRUE)

  # Dimension
  J <- normHDP_post_output$J
  G <- normHDP_post_output$G
  D <- normHDP_post_output$D
  C <- normHDP_post_output$C

  # Allocation
  Z <- normHDP_post_output$Z

  # Posterior sample
  if(type == 'standard'){

    beta.sample <- lapply(random.index,
                          function(t) normHDP_post_output$Beta_output[[t]])

    mu.sample <- lapply(random.index,
                        function(t) normHDP_post_output$mu_star_1_J_output[[t]])

    phi.sample <- lapply(random.index,
                         function(t) normHDP_post_output$phi_star_1_J_output[[t]])

    # Register cores
    if(run.on.pc == FALSE){

      cl <- makeCluster(num.cores,
                        type = "FORK",
                        .packages = 'stats')

    }else{

      cl <- makeCluster(num.cores,
                        .packages = 'stats')
    }

    # Replicated Y
    y.rep <- pblapply(1:rep.number,
                      cl = cl,
                      FUN = function(i){

                        lapply(1:D,
                               function(d){

                                 matrix(rnbinom(n = G*C[d],
                                                mu = as.vector(t(mu.sample[[i]][Z[[d]],]))*rep(beta.sample[[i]][[d]], each = G),
                                                size = as.vector(t(phi.sample[[i]][Z[[d]],]))),

                                        nrow = G,
                                        ncol = C[d])
                               })

                      })

    stopCluster(cl)


  }else if(type == 'mixed'){

    # Samples
    beta.sample <- lapply(random.index,
                          function(t) normHDP_post_output$Beta_output[[t]])

    mu.sample <- lapply(random.index,
                        function(t) normHDP_post_output$mu_star_1_J_output[[t]])

    alpha.phi2.sample <- sapply(random.index,
                                function(t) normHDP_post_output$alpha_phi2_output[t])

    b.sample <- lapply(random.index,
                       function(t) normHDP_post_output$b_output[[t]])

    # sample of phi

    if(length(b.sample[[1]]) == 3){

      phi.sample <- lapply(1:rep.number,
                           function(i){

                             matrix(rlnorm(n = J*G,
                                           meanlog = b.sample[[i]][1]+b.sample[[i]][2]*log(mu.sample[[i]])+b.sample[[i]][3]*log(mu.sample[[i]])^2,
                                           sdlog = sqrt(alpha.phi2.sample[i])),
                                    nrow = J,
                                    ncol = G)
                           })

    }else{

      phi.sample <- lapply(1:rep.number,
                           function(i){

                             matrix(rlnorm(n = J*G,
                                           meanlog = b.sample[[i]][1]+b.sample[[i]][2]*log(mu.sample[[i]]),
                                           sdlog = sqrt(alpha.phi2.sample[i])),
                                    nrow = J,
                                    ncol = G)

                           })
    }


    # Register cores
    if(run.on.pc == FALSE){

      cl <- makeCluster(num.cores,
                        type = "FORK",
                        .packages = 'stats')

    }else{

      cl <- makeCluster(num.cores,
                        .packages = 'stats')
    }

    # Replicated Y
    y.rep <- pblapply(1:rep.number,
                      cl = cl,
                      FUN = function(i){

                        lapply(1:D,
                               function(d){

                                 matrix(rnbinom(n = G*C[d],
                                                mu = as.vector(t(mu.sample[[i]][Z[[d]],]))*rep(beta.sample[[i]][[d]], each = G),
                                                size = as.vector(t(phi.sample[[i]][Z[[d]],]))),

                                        nrow = G,
                                        ncol = C[d])
                               })

                      })

    stopCluster(cl)

  }

  # Register cores
  if(run.on.pc == FALSE){

    cl <- makeCluster(num.cores,
                      type = "FORK",
                      .export = 'Discrepancy',
                      .packages = 'matrixStats')

  }else{

    cl <- makeCluster(num.cores,
                      .export = 'Discrepancy',
                      .packages = 'matrixStats')

  }


  # Discrepancy for the replicated data
  y.rep.Dis <- pblapply(1:rep.number,
                       cl = cl,
                       FUN = function(i){

                         Discrepancy(Y = y.rep[[i]],

                                     theta = list('Z' = Z,
                                                  'mu' = mu.sample[[i]],
                                                  'phi' = phi.sample[[i]],
                                                  'beta' = beta.sample[[i]])
                         )

                       })

  stopCluster(cl)
  y.rep.Dis <- do.call(rbind, y.rep.Dis)



  # Register cores
  if(run.on.pc == FALSE){

    cl <- makeCluster(num.cores,
                      type = "FORK",
                      .export = 'Discrepancy',
                      .packages = 'matrixStats')

  }else{

    cl <- makeCluster(num.cores,
                      .export = 'Discrepancy',
                      .packages = 'matrixStats')

  }

  # Discrepancy for the observed data
  y.Dis <- pblapply(1:rep.number,
                    cl = cl,
                    FUN = function(i){

                      Discrepancy(Y = Y,

                                  theta = list('Z' = Z,
                                               'mu' = mu.sample[[i]],
                                               'phi' = phi.sample[[i]],
                                               'beta' = beta.sample[[i]])
                      )

                    })

  stopCluster(cl)
  y.Dis <- do.call(rbind, y.Dis)


  #--------------------------------- Compare statistics -----------------------------

  p_value <- NULL

  for(d in 1:D){

    if(run.on.pc == FALSE){

      cl <- makeCluster(num.cores,
                        type = "FORK",
                        .packages = c('tidyr', 'tidyverse'))
    }else{

      cl <- makeCluster(num.cores,
                        .packages = c('tidyr', 'tidyverse'))
    }


    p_value[[d]] <- pblapply(1:G,
                             cl = cl,
                             function(g){


                             # For y rep
                             p_value_rep_dg <- y.rep.Dis %>%

                               filter(dataset == d,
                                      gene == g)

                             # For y obs
                             p_value_obs_dg <- y.Dis %>%

                               filter(dataset == d,
                                      gene == g)

                             # Output a data frame
                             data.frame(dataset = d,
                                        gene = g,
                                        Dis.1 = mean(p_value_rep_dg$D1 <= p_value_obs_dg$D1),
                                        Dis.2 = mean(p_value_rep_dg$D2 <= p_value_obs_dg$D2),
                                        Dis.3 = mean(p_value_rep_dg$D3 <= p_value_obs_dg$D3))

                           })

    stopCluster(cl)

    p_value[[d]] <- do.call(rbind,
                            p_value[[d]])




  }

  return(do.call(rbind,
                 p_value))

}



p_value_plot <- function(ppc_output){

  par(mfrow = c(1,3))

  # Discrepancy measure 1
  hist(ppc_output$Dis.1,
       main = "median",
       xlab = "p-values",
       freq = FALSE,
       xlim = c(0,1))

  # Discrepancy measure 2
  hist(ppc_output$Dis.2,
       main = "IQR",
       xlab = "p-values",
       freq = FALSE,
       xlim = c(0,1))

  # Discrepancy measure 2
  hist(ppc_output$Dis.3,
       main = "dropout probability",
       xlab = "p-values",
       freq = FALSE,
       xlim = c(0,1))


}

