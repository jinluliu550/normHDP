
observed_counts <- function(normHDP_post_output,
                            Y,
                            gmg_output = NULL){


  ##-- Dimensions
  D <- length(Y)
  C <- sapply(Y, ncol)
  G <- nrow(Y[[1]])
  J <- normHDP_post_output$J

  Z <- normHDP_post_output$Z

  ##-- Create an empty matrix
  matrix.output <- matrix(0,
                          ncol = 0,
                          nrow = G)

  unique.Z <- sort(unique(unlist(Z)))

  index <- NULL


  ##------------------------- Reorder by cells ------------------------------

  for(clust in unique.Z){
    for(d in 1:D){

      matrix.output <- cbind(matrix.output, Y[[d]][,which(Z[[d]] == clust)])

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
                'Z' = Z))

  }else{

    ##--------------------- If also reorder by genes ---------------------

    ##-- Reorder genes by trail probability
    marker_DE <- arrange(gmg_output$marker_DE,
                         tail.prob.mu)

    marker_DD <- arrange(gmg_output$marker_DD,
                         tail.prob.phi)

    ##-- Reorder by tail probabilities
    Y_DE <- matrix.output[marker_DE$gene,]
    Y_DD <- matrix.output[marker_DD$gene,]

    ##-- Number of DE and DD genes
    DE_number <- length(which(marker_DE$tail.prob.mu >= gmg_output$alpha_M))
    DD_number <- length(which(marker_DD$tail.prob.phi >= gmg_output$alpha_D))

    ##-- Return
    return(list('Y_DE' = Y_DE,
                'Y_DD' = Y_DD,
                'DE_number' = DE_number,
                'DD_number' = DD_number,
                'index.solid' = index.solid,
                'index.dashed' = index.dashed,
                'Z' = Z))
  }
}


observed_counts_plot <- function(observed_counts_output,
                                 width.plot = 1000,
                                 height.plot = 400,
                                 file.name.1,
                                 file.name.2 = NULL){

  Z <- observed_counts_output$Z



  if(!is.null(observed_counts_output$DE_number)){

    # Plot for DE
    png(filename = file.name.1,
        width = width.plot,
        height = height.plot)

    image.plot(1:length(unlist(Z)),
               1:nrow(observed_counts_output$Y_DE),
               as.matrix(t(log(observed_counts_output$Y_DE+1))),
               xlab = 'Cell',
               ylab = 'Genes',
               main = 'Observed counts of HET and HOM')

    abline(v = observed_counts_output$index.solid, lty = 1, lwd = 1.5, col = 'yellow')
    abline(v = observed_counts_output$index.dashed, lty = 2, lwd = 1.5, col = 'yellow')
    abline(h = nrow(observed_counts_output$Y_DE) - observed_counts_output$DE_number, lty = 2, lwd = 2, col = 'red')
    dev.off()



    # Plot for DD
    png(filename = file.name.2,
        width = width.plot,
        height = height.plot)

    image.plot(1:length(unlist(Z)),
               1:nrow(observed_counts_output$Y_DD),
               as.matrix(t(log(observed_counts_output$Y_DD+1))),
               xlab = 'Cell',
               ylab = 'Genes',
               main = 'Observed counts of HET and HOM')
    abline(v = observed_counts_output$index.solid, lty = 1, lwd = 1.5, col = 'yellow')
    abline(v = observed_counts_output$index.dashed, lty = 2, lwd = 1.5, col = 'yellow')
    abline(h = nrow(observed_counts_output$Y_DD) - observed_counts_output$DD_number, lty = 2, lwd = 2, col = 'red')

    dev.off()

  }else{

    png(filename = file.name.1,
        width = width.plot,
        height = height.plot)

    image.plot(1:length(unlist(Z)),
               1:nrow(observed_counts_output$matrix.output),
               as.matrix(t(log(observed_counts_output$matrix.output+1))),
               xlab = 'Cell',
               ylab = 'Genes',
               main = 'Observed counts of HET and HOM')

    abline(v = observed_counts_output$index.solid, lty = 1, lwd = 1.5, col = 'yellow')
    abline(v = observed_counts_output$index.dashed, lty = 2, lwd = 1.5, col = 'yellow')
    dev.off()

  }

}

