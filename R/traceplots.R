traceplots <- function(normHDP_mcmc_output,
                       normHDP_post_mcmc_output,
                       parameter,
                       j = NULL,
                       g = NULL,
                       c = NULL,
                       d = NULL){

  # Overall chain length
  chain.length <- length(normHDP_mcmc_output$b_output)
  chain.length.post <- length(normHDP_post_mcmc_output$mu_star_1_J_output)

  # Trace plot of alpha
  if(parameter == 'alpha'){

    plot(1:chain.length,
         normHDP_mcmc_output$alpha_output,
         xlab = 'Index',
         ylab = 'alpha')
  }

  # Trace plot of alpha_zero
  if(parameter == 'alpha_zero'){

    plot(1:chain.length,
         normHDP_mcmc_output$alpha_zero_output,
         xlab = 'Index',
         ylab = 'alpha_zero')
  }

  # Trace plot of beta
  if(parameter == 'beta'){

    plot(1:chain.length,
         sapply(1:chain.length, function(i) normHDP_mcmc_output$Beta_output[[i]][[d]][c]),
         xlab = 'Index',
         ylab = paste('beta', d, c))
  }

  # Trace plot of b
  if(parameter == 'b'){

    length.b <- length(normHDP_mcmc_output$b_output[[1]])
    for(q in 1:length.b){

      plot(1:chain.length,
           sapply(1:chain.length, function(i) normHDP_mcmc_output$b_output[[i]][q]),
           xlab = 'Index',
           ylab = paste('b', q-1))
    }
  }

  # Trace plot of alpha_phi_2
  if(parameter == 'alpha_phi_2'){

    plot(1:chain.length,
         sapply(1:chain.length, function(i) normHDP_mcmc_output$alpha_phi2_output),
         xlab = 'Index',
         ylab = 'alpha_phi2')
  }

  # Trace plot of mu
  if(parameter == 'mu'){

    plot(1:chain.length.post,
         sapply(1:chain.length.post,
                function(i) normHDP_post_mcmc_output$mu_star_1_J_output[[i]][j,g]),
         xlab = 'Index',
         ylab = paste('mu', j, g))
  }

  # Trace plot of phi
  if(parameter == 'phi'){

    plot(1:chain.length.post,
         sapply(1:chain.length.post,
                function(i) normHDP_post_mcmc_output$phi_star_1_J_output[[i]][j,g]),
         xlab = 'Index',
         ylab = paste('phi', j, g))
  }

}
