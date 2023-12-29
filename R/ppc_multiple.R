
# Function 1: Compute statistics of the observed and replicated datasets
ppc_multiple <- function(normHDP_post_output,
                         Y,
                         number_rep,
                         num.cores,
                         run.on.pc = TRUE){


  # Input data
  D <- normHDP_post_output$D
  G <- normHDP_post_output$G
  C <- normHDP_post_output$C
  J <- normHDP_post_output$J
  Z <- normHDP_post_output$Z

  # Trace from MCMC
  mu_star_1_J_output <- normHDP_post_output$mu_star_1_J_output
  alpha_phi_2_output <- normHDP_post_output$alpha_phi2_output
  Beta_output <- normHDP_post_output$Beta_output
  b_output <- normHDP_post_output$b_output


  # Generate some random index of length number_rep
  index <- sample(1:length(mu_star_1_J_output),
                  number_rep,
                  replace = FALSE)

  ##----------------------- Compute discrepancy measures ----------------------------------

  if(run.on.pc == FALSE){

    cl <- makeCluster(num.cores,
                      type = "FORK",
                      .packages = 'stats')

  }else{

    cl <- makeCluster(num.cores,
                      .packages = 'stats')
  }

  rep_Y_statistics <- pblapply(1:number_rep,
                               cl = cl,
                               FUN = function(i){


                                 t <- index[i]

                                 # Set theta
                                 theta <- list('mu' = mu_star_1_J_output[[t]],
                                               'b' = b_output[[t]],
                                               'alpha_phi_2' = alpha_phi_2_output,
                                               'Beta' = Beta_output[[t]])

                                 # Simulate phi
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
                                 Y_rep <- lapply(1:D,
                                                 function(d){

                                                   matrix(rnbinom(n = G*C[d],
                                                                  mu = as.vector(t(theta$mu[Z[[d]],]))*rep(theta$Beta[[d]], each = G),
                                                                  size = as.vector(t(phi[Z[[d]],]))),

                                                          nrow = G,
                                                          ncol = C[d])
                                                 })

                                 # Statistics for replicated Y
                                 Y_rep_statistics_t <- lapply(1:D, function(d){

                                   data.frame(dataset = d,
                                              gene = 1:G,
                                              mean.log.shifted.counts = apply(Y_rep[[d]], 1, function(x) mean(log(x+1))),
                                              sd.log.shifted.counts = apply(Y_rep[[d]], 1, function(x) sd(log(x+1))),
                                              log.mean.counts = apply(Y_rep[[d]], 1, function(x) log(mean(x))),
                                              dropout.probability = apply(Y_rep[[d]], 1, function(x) length(which(x == 0))/C[d]))
                                 })

                                 # Combine all data
                                 Y_rep_statistics_t <- do.call(rbind,
                                                               Y_rep_statistics_t)
                                 Y_rep_statistics_t$t <- i

                                 return(Y_rep_statistics_t)
                               })


  stopCluster(cl)

  # Combine all replicated data
  rep_Y_statistics <- do.call(rbind,
                              rep_Y_statistics)

  ##-- For observed Y
  Y_statistics <- lapply(1:D, function(d){

    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,
                           mean.log.shifted.counts = apply(Y[[d]],
                                                           1,
                                                           function(x) mean(log(x+1))),

                           sd.log.shifted.counts = apply(Y[[d]],
                                                         1,
                                                         function(x) sd(log(x+1))),

                           log.mean.counts = apply(Y[[d]],
                                                   1,
                                                   function(x) log(mean(x))),

                           dropout.probability = apply(Y[[d]],
                                                       1,
                                                       function(x) length(which(x == 0))/C[d]))
  })

  Y_statistics <- do.call(rbind,
                          Y_statistics)
  Y_statistics$t <- 0

  ##-- Return final output
  return(list('rep_Y_statistics' = rep_Y_statistics,
              'Y_statistics' = Y_statistics,
              'number_rep' = number_rep))
}


# Function 2: Compare statistics of the observed data against replicated data
ppc_multiple_plot <- function(ppc_multiple_df,
                              title){

  # Statistics
  Y_statistics <- ppc_multiple_df$Y_statistics
  rep_Y_statistics <- ppc_multiple_df$rep_Y_statistics


  Y_statistics$t <- 'observed data'
  rep_Y_statistics$t <- paste('replicated data', rep_Y_statistics$t)


  # Dimension
  D <- max(Y_statistics$dataset)

  # Kernel density plot
  df <- rbind(Y_statistics,
              rep_Y_statistics)

  colnames(df)[7] <- 'source'

  for(d in 1:D){

    df_d <- df %>%
      filter(dataset == d)

    plot1 <- ggplot()+
      geom_density(mapping = aes(x = mean.log.shifted.counts,
                                 colour = source),

                   data = df_d,
                   size = 1.2)+
      theme_bw()+
      xlab('mean of log shifted counts')+
      scale_color_manual(values=c(rep('grey',ppc_multiple_df$number_rep), 'red'))+
      theme(legend.position = "none")+
      ggtitle(title[d])


    plot2 <- ggplot()+
      geom_density(mapping = aes(x = sd.log.shifted.counts,
                                 colour = source),

                   data = df_d,
                   size = 1.2)+
      theme_bw()+
      xlab('standard deviaion of log shifted counts')+
      scale_color_manual(values=c(rep('grey',ppc_multiple_df$number_rep), 'red'))+
      theme(legend.position = "none")+
      ggtitle(title[d])


    plot3 <- ggplot()+
      geom_density(mapping = aes(x = dropout.probability,
                                 colour = source),

                   data = df_d,
                   size = 1.2)+
      theme_bw()+
      xlab('dropout probabilities')+
      scale_color_manual(values=c(rep('grey',ppc_multiple_df$number_rep), 'red'))+
      theme(legend.position = "none")+
      ggtitle(title[d])

    ggarrange(plot1, plot2, plot3, nrow = 1)

  }



}
