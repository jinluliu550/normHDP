
#' Plots for posterior predictive checks with single replicate.
#' 
#' @param normHDP_output Output from normHDP_mcmc.
#' @param Y Input dataset.
#' @param data.name Name of each dataset.
#' @return For each dataset, compare relationship between mean of log shifted counts and standard deviation of log shifted counts
#' for the observed dataset against the replicated dataset. In addition, compare relationship between log of mean counts and dropout
#' probabilities between observed and replicated dataset. Also compare the statistics gene-wise.
#' @export
ppc_single <- function(normHDP_output,
                       Y,
                       data.name){
  
  ##-- Dimensions
  G <- normHDP_output$G
  D <- normHDP_output$D
  C <- normHDP_output$C
  J <- normHDP_output$J
  
  ##-- Allocation probabilities of each cell in each sample
  allocation_prob_output <- normHDP_output$allocation_prob_output
  
  ##-- Imports
  mu_star_1_J_output <- normHDP_output$mu_star_1_J_output
  b_output <- normHDP_output$b_output
  alpha_phi2_output <- normHDP_output$alpha_phi2_output
  Beta_output <- normHDP_output$Beta_output
  Z_output <- normHDP_output$Z_output
  
  mcmc.length <- length(Z_output)
  
  
  ##------------------------ Find cells with low allocation probabilities ----------------------
  
  ##-- Find cells with high mean allocation probability
  high.allo.index <- NULL
  
  for(d in 1:D){
    
    ##-- For each dataset
    list.result <- lapply(1:mcmc.length, 
                          function(i) allocation_prob_output[[i]][[d]])
    df.result <- do.call(rbind, list.result)
    
    ##-- Mean allocation probability of each cell
    allo.prob.mean <- colSums(df.result)/nrow(df.result)
    
    ##-- Set a threshold
    threshold <- 0.8
    
    ##-- Find cells with mean allocation probability greater than the threshold
    high.allo.index[[d]] <- which(allo.prob.mean > threshold)
  }
  
  ##-- Number of cells with high allocation probabilities
  C_new <- sapply(1:D, 
                  function(d) length(high.allo.index[[d]]))
  
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
    
    ##-- Cells with high allocation probabilities
    Y_d <- Y[[d]][,high.allo.index[[d]]]
    
    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,
                           
                           ##-- Relationship (1)
                           mean.log.shifted.counts = apply(Y_d, 
                                                           1,
                                                           function(x) mean(log(x+1))),
                           sd.log.shifted.counts = apply(Y_d,
                                                         1,
                                                         function(x) sd(log(x+1))),
                           
                           ##-- Relationship (2)
                           log.mean.counts = apply(Y_d,
                                                   1,
                                                   function(x) log(mean(x))),
                           dropout.probability = apply(Y_d, 
                                                       1,
                                                       function(x) length(which(x == 0))/C[d]))
    
    
    ##-- Append for each dataset
    rel.df <- rbind(rel.df, rel.df.d)
  }
  
  ##-- Add a source
  rel.df$source <- 'Observed Y'
  
  ##------------------------------------------------------------------------------------
  
  ##-- Replicated data
  index <- sample(1:mcmc.length, 1)
  mu.sample <- mu_star_1_J_output[[index]]
  b.sample <- b_output[[index]]
  alpha_phi2.sample <- alpha_phi2_output[index]
  
  Z.sample <- lapply(1:D, 
                     function(d) Z_output[[index]][[d]][high.allo.index[[d]]])
  
  Beta.sample <- lapply(1:D, 
                        function(d) Beta_output[[index]][[d]][high.allo.index[[d]]])
  
  
  ##----------------------------------------------------------------------------
  
  ##-- Find phi
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
  Y_rep <- NULL
  
  for(d in 1:D){
    
    loop.result <- lapply(1:C_new[d], function(c){
      
      as.vector(rnbinom(n=G, 
                        mu = mu.sample[Z.sample[[d]][c],]*Beta.sample[[d]][c], 
                        size = phi[Z.sample[[d]][c],]))
    })
    
    ##-- Combine all items from the list
    loop.result <- do.call(rbind, loop.result)
    Y_rep[[d]] <- t(loop.result)
  }
  
  ##-- Find the corresponding statistics for Y_rep
  rel.rep.df <- data.frame(dataset = numeric(),
                           gene = numeric(),
                           mean.log.shifted.counts = numeric(),
                           sd.log.shifted.counts = numeric(),
                           log.mean.counts = numeric(),
                           dropout.probability = numeric())
  
  for(d in 1:D){
    
    Y_rep_d <- Y_rep[[d]]
    
    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,
                           
                           mean.log.shifted.counts = apply(Y_rep_d,
                                                           1,
                                                           function(x) mean(log(x+1))),
                           sd.log.shifted.counts = apply(Y_rep_d, 
                                                         1,
                                                         function(x) sd(log(x+1))),
                           
                           log.mean.counts = apply(Y_rep_d,
                                                   1,
                                                   function(x) log(mean(x))),
                           dropout.probability = apply(Y_rep_d,
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
    
    ##-- Plot
    plot(x = rel.df.d$mean.log.shifted.counts,
         y = rel.df.d$sd.log.shifted.counts,
         xlab = 'Mean of log-shifted counts',
         ylab = 'Standard deviation of log-shifted counts',
         pch = 19,
         cex = 0.2,
         main = data.name[d],
         col = 'black',
         cex.lab = 1.5)
    
    points(x = rel.rep.df.d$mean.log.shifted.counts,
           y = rel.rep.df.d$sd.log.shifted.counts,
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
    
    plot(x = rel.df.d$log.mean.counts,
         y = rel.df.d$dropout.probability,
         xlab = 'log of mean gene counts',
         ylab = 'Dropout probabilities',
         pch = 19,
         cex = 0.2,
         main = data.name[d],
         col = 'black',
         cex.lab = 1.5)
    
    points(x = rel.rep.df.d$log.mean.counts,
           y = rel.rep.df.d$dropout.probability,
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
