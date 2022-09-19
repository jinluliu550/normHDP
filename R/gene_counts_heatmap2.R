
#' Observed gene counts
#'
#' @param Z_CD Cell allocations.
#' @param Y Input dataset.
#' @param gmg_output Output from global_marker_genes. If this input is provided, the function will reorder datasets by both
#' the tail probabilities and clustering estimate. If this input is not provided, the function will only reorder cells
#' based on the clustering estimate.
#' @return If gmg_output is not provided, the function return a list of 3 items. The observed gene counts reordered by Z, the index
#' to separate cells from different clusters (index.solid), and the index to separate cells from different datasets (index.dashed).
#' If gmg_output is provided, the function return a list of 6 items. The first 2 items are the observed counts
#' with genes reordered by the tail probabilities corresponding to mean expression and dispersion, respectively. Item 3
#' shows the total number of DE genes. Item 4 shows the total number of DD genes. Item 5 and 6 are index.solid and index.dashed,
#' as described above.
#' @export
observed_counts <- function(Z_CD,
                            Y,
                            gmg_output = NULL){


  ##-- Dimensions
  D <- length(Y)
  C <- sapply(Y, ncol)
  G <- nrow(Y[[1]])

  ##-- Create an empty matrix
  matrix.output <- matrix(0,
                          ncol = 0,
                          nrow = G)

  unique.Z <- sort(unique(unlist(Z_CD)))

  index <- NULL


  ##------------------------- Reorder by cells ------------------------------

  for(clust in unique.Z){
    for(d in 1:D){

      matrix.output <- cbind(matrix.output, Y[[d]][,which(Z_CD[[d]] == clust)])

      ##-- Update index
      index <- c(index, length(which(Z_CD[[d]] == clust)))
    }
  }


  ##-- Vertical line solid
  index.solid <- index[seq(1, length(index), D)]
  index.solid <- index.solid[-length(index.solid)]

  ##-- Vertical line dashed
  index.dashed <- index[which(index %notin% index.solid)]

  ##-- If only reorder by cells
  if(is.null(gmg_output)){

    return(list('matrix.output' = matrix.output,
                'index.solid' = index.solid,
                'index.dashed' = index.dashed))

  }else{

    ##--------------------- If also reorder by genes ---------------------

    ##-- Reorder genes by trail probability
    marker_DE <- arrange(gmg_output$marker_DE,
                         tail.prob.mu)

    marker_DD <- arrange(gmg_output$marker_DD,
                         tail.prob.phi)

    ##-- Reorder by tail probabilities
    Y_DE <- lapply(1:D,
                   function(d) matrix.output[marker_DE$gene,])

    Y_DD <- lapply(1:D,
                   function(d) matrix.output[marker_DE$gene,])

    ##-- Number of DE and DD genes
    DE_number <- length(which(marker_DE$tail.prob.mu >= gmg_output$alpha_M))
    DD_number <- length(which(marker_DD$tail.prob.phi >= gmg_output$alpha_D))

    ##-- Return
    return(list('Y_DE' = Y_DE,
                'Y_DD' = Y_DD,
                'DE_number' = DE_number,
                'DD_number' = DD_number,
                'index.solid' = index.solid,
                'index.dashed' = index.dashed))
  }
}



