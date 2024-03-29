

ppc_single_plot <- function(normHDP_post_output,
                            Y,
                            data.name){


  ##-- Dimensions
  G <- normHDP_post_output$G
  D <- normHDP_post_output$D
  C <- normHDP_post_output$C
  J <- normHDP_post_output$J


  ##-- Imports
  mu_star_1_J_output <- normHDP_post_output$mu_star_1_J_output
  b_output <- normHDP_post_output$b_output
  alpha_phi2_output <- normHDP_post_output$alpha_phi2_output
  Beta_output <- normHDP_post_output$Beta_output
  Z <- normHDP_post_output$Z

  mcmc.length <- length(mu_star_1_J_output)


  ##------------------------------------------------------------------------------------------------


  ##-- Data frame for the real data
  rel.df <- data.frame(dataset = numeric(),
                       gene = numeric(),
                       mean.log.shifted.counts = numeric(),
                       sd.log.shifted.counts = numeric(),
                       log.mean.counts = numeric(),
                       dropout.probability = numeric())

  ##-- Relationship (1) between mean and standard deviation of log shifted counts
  ##-- Relationship (2) between log of mean counts and dropout probabilities

  for(d in 1:D){

    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,

                           ##-- Relationship (1)
                           mean.log.shifted.counts = apply(Y[[d]],
                                                           1,
                                                           function(x) mean(log(x+1))),
                           sd.log.shifted.counts = apply(Y[[d]],
                                                         1,
                                                         function(x) sd(log(x+1))),

                           ##-- Relationship (2)
                           log.mean.counts = apply(Y[[d]],
                                                   1,
                                                   function(x) log(mean(x))),

                           dropout.probability = apply(Y[[d]],
                                                       1,
                                                       function(x) length(which(x == 0))/C[d]))


    ##-- Append for each dataset
    rel.df <- rbind(rel.df, rel.df.d)
  }

  ##-- Add a source
  rel.df$source <- 'Observed Y'

  ##------------------------------------------------------------------------------------

  # Parameters
  t <- sample(1:mcmc.length, 1)

  mu.sample <- mu_star_1_J_output[[t]]
  b.sample <- b_output[[t]]
  alpha_phi2.sample <- alpha_phi2_output[t]
  Beta.sample <- Beta_output[[t]]



  ##----------------------------------------------------------------------------

  # Simulate phi
  if(length(b.sample) == 3){

    ##-- If quadratic
    phi <- matrix(rlnorm(n = J*G,
                         meanlog = b.sample[1]+b.sample[2]*log(mu.sample)+b.sample[3]*log(mu.sample)^2,
                         sdlog = sqrt(alpha_phi2.sample)),
                  nrow = J,
                  ncol = G)
  }else{

    ##-- If linear
    phi <- matrix(rlnorm(n = J*G,
                         meanlog = as.vector(b.sample[1]+b.sample[2]*log(mu.sample)),
                         sdlog = sqrt(alpha_phi2.sample)),
                  nrow = J,
                  ncol = G)
  }


  ##-- Simulate Y_rep
  Y_rep <- lapply(1:D,
                  function(d){

                    matrix(rnbinom(n = G*C[d],
                                   mu = as.vector(t(mu.sample[Z[[d]],]))*rep(Beta.sample[[d]], each=G),
                                   size = as.vector(t(phi[Z[[d]],]))
                                   ),

                           nrow = G,
                           ncol = C[d]
                           )
                  })


  ##-- Find the corresponding statistics for Y_rep
  rel.rep.df <- data.frame(dataset = numeric(),
                           gene = numeric(),
                           mean.log.shifted.counts = numeric(),
                           sd.log.shifted.counts = numeric(),
                           log.mean.counts = numeric(),
                           dropout.probability = numeric())

  for(d in 1:D){


    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,

                           mean.log.shifted.counts = apply(Y_rep[[d]],
                                                           1,
                                                           function(x) mean(log(x+1))),

                           sd.log.shifted.counts = apply(Y_rep[[d]],
                                                         1,
                                                         function(x) sd(log(x+1))),

                           log.mean.counts = apply(Y_rep[[d]],
                                                   1,
                                                   function(x) log(mean(x))),

                           dropout.probability = apply(Y_rep[[d]],
                                                       1,
                                                       function(x) length(which(x == 0))/C[d]))

    ##-- Append for each dataset
    rel.rep.df <- rbind(rel.rep.df, rel.df.d)
  }

  rel.rep.df$source <- 'Replicated Y'

  ##-- Combine rel.df with rel.rep.df
  df <- rbind(rel.df, rel.rep.df)

  ##----------------------------- ggplot -----------------------------------

  ##-- Plot 2 relationships for each dataset

  ##-- Set dimension
  par(mfrow=c(1,2*length(data.name)))

  ##-- Relationship 1
  for(d in 1:D){

    rel.df.d <- filter(rel.df,
                       dataset == d)

    rel.rep.df.d <- filter(rel.rep.df,
                           dataset == d)

    # Replicated data in black
    plot(x = rel.rep.df.d$mean.log.shifted.counts,
         y = rel.rep.df.d$sd.log.shifted.counts,
         xlab = 'Mean of log-shifted counts',
         ylab = 'Standard deviation of log-shifted counts',
         pch = 19,
         cex = 0.2,
         main = data.name[d],
         col = 'black',
         cex.lab = 1.5)

    # Observed data in red
    points(x = rel.df.d$mean.log.shifted.counts,
           y = rel.df.d$sd.log.shifted.counts,
           pch = 19,
           cex = 0.2,
           col = 'red')
  }

  ##-- Relationship 2
  for(d in 1:D){

    rel.df.d <- filter(rel.df,
                       dataset == d)

    rel.rep.df.d <- filter(rel.rep.df,
                           dataset == d)

    # Replicated data in black
    plot(x = rel.rep.df.d$log.mean.counts,
         y = rel.rep.df.d$dropout.probability,
         xlab = 'log of mean gene counts',
         ylab = 'Dropout probabilities',
         pch = 19,
         cex = 0.2,
         main = data.name[d],
         col = 'black',
         cex.lab = 1.5)

    # Observed data in red
    points(x = rel.df.d$log.mean.counts,
           y = rel.df.d$dropout.probability,
           pch = 19,
           cex = 0.2,
           col = 'red')
  }

  ##-- Gene-wise comparison

  ##-- Gene-wise comparison of statistics
  df2 <- data.frame(mean.log.shifted.counts.observed = rel.df$mean.log.shifted.counts,
                    sd.log.shifted.counts.observed = rel.df$sd.log.shifted.counts,
                    dropout.probability.observed = rel.df$dropout.probability,

                    mean.log.shifted.counts.fitted = rel.rep.df$mean.log.shifted.counts,
                    sd.log.shifted.counts.fitted = rel.rep.df$sd.log.shifted.counts,
                    dropout.probability.fitted = rel.rep.df$dropout.probability)

  ##-- Plots on the 2nd row
  par(mfrow=c(1,3))

  plot(x = df2$mean.log.shifted.counts.observed,
       y = df2$mean.log.shifted.counts.fitted,
       xlab = 'Observed dataset',
       ylab = 'Replicated dataset',
       main = 'Mean of log-shifted counts',
       pch = 19,
       cex = 0.2,
       cex.lab = 1.5,
       col = 'grey')

  abline(a = 0, b = 1, lty = 2)


  plot(x = df2$sd.log.shifted.counts.observed,
       y = df2$sd.log.shifted.counts.fitted,
       xlab = 'Observed dataset',
       ylab = 'Replicated dataset',
       main = 'Standard deviation of log-shifted counts',
       pch = 19,
       cex = 0.2,
       cex.lab = 1.5,
       col = 'grey')

  abline(a = 0, b = 1, lty = 2)


  plot(x = df2$dropout.probability.observed,
       y = df2$dropout.probability.fitted,
       xlab = 'Observed dataset',
       ylab = 'Replicated dataset',
       main = 'Dropout probabilities',
       pch = 19,
       cex = 0.2,
       cex.lab = 1.5,
       col = 'grey')

  abline(a = 0, b = 1, lty = 2)


}
