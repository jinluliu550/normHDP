consensus_clustering <- function(Z_output,
                                 num.cores,
                                 run.on.pc = TRUE){

  # Width
  Width <- length(Z_output)

  # Length
  Length <- length(Z_output[[1]])

  # Z_trace
  Z_trace <- lapply(1:Width,
                    function(w) Z_output[[w]][[Length]])

  # psm
  psm <- similarity_matrix(Z_trace = Z_trace,
                           num.cores = num.cores,
                           run.on.pc = run.on.pc)

  # Clustering estimate
  cluster.estimate <- opt.clust(Z_trace = Z_trace,
                                psm_output =  psm)

  return(list('cluster.estimate' = cluster.estimate,
              'psm' = psm))
}

consensus_clustering_plot <- function(Z_output,
                                      W.breaks,
                                      L.breaks,
                                      num.cores,
                                      run.on.pc = TRUE){

  W.breaks <- W.breaks
  L.breaks <- L.breaks-1


  # df0
  df0 <- data.frame(width = numeric(),
                    length = numeric(),
                    mean_difference = numeric())

  for(w in W.breaks){

    for(l in L.breaks){

      Z.trace.wl <- lapply(1:w,
                           function(w0) Z_output[[w0]][[l]]
                           )

      Z.trace.wl_minus_1 <- lapply(1:w,
                                   function(w0) Z_output[[w0]][[l-1]])

      psm.wl <- similarity_matrix(Z_trace = Z.trace.wl,
                                  num.cores = num.cores,
                                  run.on.pc = run.on.pc)

      psm.wl_minus_1 <- similarity_matrix(Z_trace = Z.trace.wl_minus_1,
                                          num.cores = num.cores,
                                          run.on.pc = run.on.pc)

      # Update df0
      df0 <- rbind(df0,

                   data.frame(width = w,
                              length = l+1,
                              mean_difference = mean(abs(psm.wl$psm.combined - psm.wl_minus_1$psm.combined))))
    }
  }

  # Plot
  df0$width <- as.factor(df0$width)

  df0 %>%
    ggplot(mapping = aes(x = length, y = mean_difference, group = width))+
    geom_line(aes(color = width))+
    geom_point(aes(color = width))+
    theme_bw()+
    xlab('Chain depth')+
    ylab('Mean aboslute difference of similarity')

}



