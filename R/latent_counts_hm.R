
latent_counts <- function(Y,
                          normHDP_post_output,
                          gmg_output = NULL,
                          num.cores = 1,
                          run.on.pc = TRUE){

  ##-- Dimensions
  D <- length(Y)
  C <- sapply(1:D, function(d){ncol(Y[[d]])})
  G <- nrow(Y[[1]])
  J <- normHDP_post_output$J


  #-- Posterior sample
  mu_sample <- normHDP_post_output$mu_star_1_J_output
  phi_sample <- normHDP_post_output$phi_star_1_J_output
  beta_sample <- normHDP_post_output$Beta_output
  Z <- normHDP_post_output$Z

  sample.length <- length(mu_sample)

  #-- Latent Y
  Y_latent <- NULL

  for(d in 1:D){

    if(run.on.pc == FALSE){

      cl <- makeCluster(num.cores,
                        type = "FORK")

    }else{

      cl <- makeCluster(num.cores)
    }


    loop.result <- pblapply(1:sample.length,
                            cl = cl,
                            FUN = function(t){

                              matrix(as.vector(Y[[d]])*(as.vector(t(mu_sample[[t]][Z[[d]],]))+as.vector(t(phi_sample[[t]][Z[[d]],])))/
                                       (as.vector(t(mu_sample[[t]][Z[[d]],]))*rep(beta_sample[[t]][[d]],each=G)+as.vector(t(phi_sample[[t]][Z[[d]],]))) +

                                       as.vector(t(mu_sample[[t]][Z[[d]],]))*(as.vector(t(phi_sample[[t]][Z[[d]],]))*(1-rep(beta_sample[[t]][[d]],each=G)))/
                                       (as.vector(t(mu_sample[[t]][Z[[d]],]))*rep(beta_sample[[t]][[d]],each=G)+as.vector(t(phi_sample[[t]][Z[[d]],]))),

                                     nrow = G,
                                     ncol = C[d]
                              )

                            })

    stopCluster(cl)

    Y_latent[[d]] <- Reduce('+',loop.result)/sample.length
  }






  ##-- Create an empty matrix
  matrix.output <- matrix(0,
                          ncol = 0,
                          nrow = G)

  unique.Z <- 1:J

  index <- NULL

  ##------------------------- Reorder by cells ------------------------------

  for(clust in unique.Z){
    for(d in 1:D){

      matrix.output <- cbind(matrix.output, Y_latent[[d]][,which(Z[[d]] == clust)])

      ##-- Update index
      index <- c(index, length(which(Z[[d]] == clust)))
    }
  }


  ##-- Vertical line solid
  index <- cumsum(index)
  index.solid <- index[seq(from = D,
                                 to = D*J,
                                 D)]
  index.solid <- index.solid[-length(index.solid)]


  ##-- Vertical line dashed
  index.dashed <- index[which(index %notin% index.solid)]

  ##-- If only reorder by cells
  if(is.null(gmg_output)){

    return(list('matrix.output' = matrix.output,
                'index.solid' = index.solid,
                'index.dashed' = index.dashed,
                'Y_latent' = Y_latent,
                'Z' = Z,
                'G' = G))
  }else{

    ##-------------- If also reorder by gene ----------------

    ##-- Reorder genes by trail probability
    marker_DE <- arrange(gmg_output$marker_DE,
                         tail.prob.mu)

    marker_DD <- arrange(gmg_output$marker_DD,
                         tail.prob.phi)

    ##-- Reorder by tail probabilities
    Y_latent_DE <- matrix.output[marker_DE$gene,]
    Y_latent_DD <- matrix.output[marker_DD$gene,]

    ##-- Number of DE and DD genes
    DE_number <- length(which(marker_DE$tail.prob.mu >= gmg_output$alpha_M))
    DD_number <- length(which(marker_DD$tail.prob.phi >= gmg_output$alpha_D))

    ##-- Return
    return(list('Y_latent_DE' = Y_latent_DE,
                'Y_latent_DD' = Y_latent_DD,
                'DE_number' = DE_number,
                'DD_number' = DD_number,
                'index.solid' = index.solid,
                'index.dashed' = index.dashed,
                'Y_latent' = Y_latent,
                'Z' = Z,
                'G' = G))
  }

}


latent_counts_plot <- function(latent_counts_output,
                               width.plot = 1000,
                               height.plot = 400,
                               file.name.1,
                               file.name.2 = NULL){

  Z <- latent_counts_output$Z
  G <- latent_counts_output$G

  if(!is.null(latent_counts_output$DE_number)){

    # Plot for DE
    png(filename = file.name.1,
        width = width.plot,
        height = height.plot)

    image.plot(1:length(unlist(Z)),
               1:G,
               t(log(latent_counts_output$Y_latent_DE+1)),
               xlab = 'Cell',
               ylab = 'Genes',
               main = 'Latent counts of HET and HOM')

    abline(v = latent_counts_output$index.solid, lty = 1, lwd = 1.5, col = 'yellow')
    abline(v = latent_counts_output$index.dashed, lty = 2, lwd = 1.5, col = 'yellow')
    abline(h = G - latent_counts_output$DE_number, lty = 2, lwd = 2, col = 'red')
    dev.off()

    # Plot for DD
    png(filename = file.name.2,
        width = width.plot,
        height = height.plot)

    image.plot(1:length(unlist(Z)),
               1:G,
               t(log(latent_counts_output$Y_latent_DD+1)),
               xlab = 'Cell',
               ylab = 'Genes',
               main = 'Latent counts of HET and HOM')
    abline(v = latent_counts_output$index.solid, lty = 1, lwd = 1.5, col = 'yellow')
    abline(v = latent_counts_output$index.dashed, lty = 2, lwd = 1.5, col = 'yellow')
    abline(h = G - latent_counts_output$DD_number, lty = 2, lwd = 2, col = 'red')
    dev.off()

  }else{


    png(filename = file.name.1,
        width = width.plot,
        height = height.plot)

    image.plot(1:length(unlist(Z)),
               1:G,
               t(log(latent_counts_output$matrix_output+1)),
               xlab = 'Cell',
               ylab = 'Genes',
               main = 'Latent counts of HET and HOM')

    abline(v = latent_counts_output$index.solid, lty = 1, lwd = 1.5, col = 'yellow')
    abline(v = latent_counts_output$index.dashed, lty = 2, lwd = 1.5, col = 'yellow')
    dev.off()
  }

}

