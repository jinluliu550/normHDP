
# Function 1: Classify clusters into three types: Stable, over-represented and under-represented clusters.
beta_analysis <- function(normHDP_post_output,
                          cluster.labels){


  # Dimensions
  G <- normHDP_post_output$G
  D <- normHDP_post_output$D
  C <- normHDP_post_output$C
  Z <- normHDP_post_output$Z
  J <- normHDP_post_output$J


  # Posterior mean of beta
  Beta_mean <- lapply(1:D,
                      function(d){

                        temp_vec <- lapply(1:length(normHDP_post_output$Beta_output),
                                           function(iter){

                                             normHDP_post_output$Beta_output[[iter]][[d]]
                                           })

                        colMeans(do.call(rbind,
                                         temp_vec))

                      })

  # Cluster labels
  df <- data.frame(beta = unlist(Beta_mean),
                   cluster = factor(unlist(Z),
                                    levels = 1:J),
                   label = cluster.labels$cluster.type[unlist(Z)])


  # Box plot in ggplot
  df %>%
    ggplot(mapping = aes(x = cluster,
                         y = beta,
                         color = label,
                         fill = label))+
    geom_boxplot()+
    ylab('Estimated capture efficiency')+
    theme_bw()

}



# Function 2: Compare baynorm estimates with posterior values
beta_compare <- function(Y,
                         beta.mean = 0.06,
                         normHDP_post_output,
                         data.names){


  D <- normHDP_post_output$D
  C <- normHDP_post_output$C

  # Prior estimate beta
  baynorm.beta <- lapply(1:D, function(d){

    bayNorm::BetaFun(Data = Y[[d]],
                     MeanBETA = beta.mean)$BETA
  })

  # Posterior mean
  beta_CD <- lapply(1:D,
                    function(d){

                      beta_D <- lapply(1:length(normHDP_post_output$Beta_output),
                                       function(t) normHDP_post_output$Beta_output[[t]][[d]])

                      colMeans(do.call(rbind,
                                       beta_D))
                    })

  # df
  df <- lapply(1:D,
               function(d){

                 data.frame(baynorm.estimate = baynorm.beta[[d]],
                            posterior.estimate = beta_CD[[d]],
                            data.name = data.names[d])
               })

  df <- do.call(rbind, df)

  # Plot
  df %>%
    ggplot(mapping = aes(x = baynorm.estimate,
                         y = posterior.estimate))+
    geom_point(color = 'grey', size = 0.1)+
    facet_wrap(~data.name)+
    theme_bw()+
    geom_abline(intercept = 0, slope = 1, col = 'blue', linetype = 'dashed')+
    xlab('baynorm estimates')+
    ylab('posterior estimates')

  # Histogram
  df %>%
    pivot_longer(!data.name,
                 names_to = 'source',
                 values_to = 'beta') %>%
    ggplot()+
    geom_histogram(mapping = aes(x = beta, fill = source, color = source), position = 'identity', alpha = 0.3)+
    theme_bw()+
    facet_wrap(~data.name)+
    xlab('capture efficiencies')



}

