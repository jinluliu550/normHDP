
# Compare difference in data-specific component probabilities
# Input: Trace of data-specific component probabilities

difference_in_p_jd <- function(p_jd_trace,
                               alpha = 0.95){

  trace.length <- length(p_jd_trace)

  J <- nrow(p_jd_trace[[1]])

  # Number of datasets
  D <- ncol(p_jd_trace[[1]])

  # Unique set of cluster
  unique.j <- 1:J
  unique.j.length <- length(unique.j)

  # All possible combinations of d
  combo <- expand.grid(1:D, 1:D)
  combo <- combo[!duplicated(t(apply(combo, 1, sort))),]
  combo <- combo[!t(apply(combo, 1, function(x) x[1] == x[2])),]

  # data to save
  df0 <- lapply(unique.j,

                function(j){

                  df00 <- lapply(1:nrow(combo),

                                 function(i){

                                   data1.index <- combo[i,1]
                                   data2.index <- combo[i,2]

                                   data.j.trace <- lapply(1:trace.length,
                                                          function(s) p_jd_trace[[s]][j,c(data1.index,
                                                                                          data2.index)])

                                   data.j.trace <- do.call(rbind,
                                                           data.j.trace)

                                   # Mean absolute difference
                                   prob <- mean(data.j.trace[,1] > data.j.trace[,2])
                                   prob <- max(prob, 1-prob)
                                   mean.abs.difference <- mean(abs(data.j.trace[,1] - data.j.trace[,2]))


                                   data.frame(cluster = paste('cluster', j),
                                              data1 = data1.index,
                                              data2 = data2.index,
                                              mean.abs.difference = mean.abs.difference,
                                              prob = prob,
                                              different = ifelse(prob > 0.95, 'Yes', 'No'))
                                 })

                  df00 <- do.call(rbind,
                                  df00)
                })

  df0_all <- do.call(rbind,
                     df0)



  par(mfrow=c(1,1))

  # Plot
  plot <- df0_all %>%
    ggplot(mapping = aes(x = mean.abs.difference,
                         y = prob))+
    geom_point(mapping = aes(color = different), size = 3)+
    geom_label_repel(aes(label = cluster),
                     size = 3)+
    geom_hline(yintercept = 0.95)+
    theme_bw()+
    xlab('mean absolute difference')+
    ylab('Probability')



  # Select the cluster and data with different p_jd
  df0_significant <- df0_all %>%
    filter(different == 'Yes')

  return(list('all_obs' = df0_all,
              'significant_obs' = df0_significant,
              'plot' = plot))

}


# Function should only be used when the total number of datasets equal to 2
difference_in_p_jd2 <- function(p_jd_trace,
                                alpha = 0.95){

  trace.length <- length(p_jd_trace)

  J <- nrow(p_jd_trace[[1]])

  # Number of datasets
  D <- ncol(p_jd_trace[[1]])

  # Unique set of cluster
  unique.j <- 1:J
  unique.j.length <- length(unique.j)

  # All possible combinations of d
  combo <- expand.grid(1:D, 1:D)
  combo <- combo[!duplicated(t(apply(combo, 1, sort))),]
  combo <- combo[!t(apply(combo, 1, function(x) x[1] == x[2])),]

  # data to save
  df0 <- lapply(unique.j,

                function(j){

                  df00 <- lapply(1:nrow(combo),

                                 function(i){

                                   data1.index <- combo[i,1]
                                   data2.index <- combo[i,2]

                                   data.j.trace <- lapply(1:trace.length,
                                                          function(s) p_jd_trace[[s]][j,c(data1.index,
                                                                                          data2.index)])

                                   data.j.trace <- do.call(rbind,
                                                           data.j.trace)

                                   # Mean absolute difference
                                   prob <- mean(data.j.trace[,1] > data.j.trace[,2])
                                   mean.abs.difference <- mean(abs(data.j.trace[,1] - data.j.trace[,2]))

                                   if(prob < 0.05){

                                     cluster.type = 'over-represented'

                                   }else if(prob > 0.95){

                                     cluster.type = 'under-represented'
                                   }else{

                                     cluster.type = 'stable'
                                   }

                                   data.frame(cluster = paste('cluster', j),
                                              data1 = data1.index,
                                              data2 = data2.index,
                                              mean.abs.difference = mean.abs.difference,
                                              prob = prob,
                                              cluster.type = cluster.type)
                                 })

                  df00 <- do.call(rbind,
                                  df00)
                })

  df0_all <- do.call(rbind,
                     df0)



  par(mfrow=c(1,1))

  # Plot
  plot <- df0_all %>%
    ggplot(mapping = aes(x = mean.abs.difference,
                         y = prob))+
    geom_point(size = 3,
               mapping = aes(color = cluster.type))+
    geom_label_repel(aes(label = cluster),
                     size = 3)+
    geom_hline(yintercept = c(0.05, 0.95),
               linetype = 'dashed',
               color = 'red')+
    theme_bw()+
    xlab('mean absolute difference')+
    ylab('Probability')

  df0_significant <- df0_all %>%
    filter(prob < 0.05 | prob > 0.95)



  return(list('all_obs' = df0_all,
              'plot' = plot,
              'significant_obs' = df0_significant))

}

