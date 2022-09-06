#' Latent gene counts
#' 
#' @param mu_JG Mean expression matrix.
#' @param phi_JG Mean dispersion matrix.
#' @param Z_CD Cell allocations.
#' @param beta_CD Capture efficiencies.
#' @param Y Input dataset.
#' @param gmg_output Output from global_marker_genes. If this input is provided, the latent counts estimate will also
#' reorder genes by ascending tail probabilities. If this input is not provided, the latent counts estimate will only
#' be reordered by the allocations.
#' @return If gmg_output is not provided, the function return a list of 3 items. The estimated latent counts, the index
#' to separate cells from different clusters (index.solid), and the index to separate cells from different datasets (index.dashed).
#' If gmg_output is provided, the function return a list of 6 items. The first 2 items are the estimated latent counts
#' with genes reordered by the tail probabilities corresponding to mean expression and dispersion, respectively. Item 3
#' shows the total number of DE genes. Item 4 shows the total number of DD genes. Item 5 and 6 are index.solid and index.dashed,
#' as described above.
#' @export
latent_counts <- function(mu_JG,
                          phi_JG,
                          Z_CD,
                          beta_CD,
                          Y,
                          gmg_output = NULL){
  
  ##-- Dimensions
  D <- length(Y)
  C <- sapply(Y, ncol)
  G <- nrow(Y[[1]])
  
  ##-- Latent Y
  Y_latent <- NULL
  
  for(d in 1:D){
    
    loop.result <- lapply(1:C[d], function(c){
      
      ##-- Allocation
      Z.c <- Z_CD[[d]][c]
      
      matrix(as.vector(Y[[d]][,c])*(as.vector(mu_JG[Z.c,])+as.vector(phi_JG[Z.c,]))/
               (as.vector(mu_JG[Z.c,])*beta_CD[[d]][c]+as.vector(phi_JG[Z.c,])) + 
               
               as.vector(mu_JG[Z.c,])*(as.vector(phi_JG[Z.c,])*(1-beta_CD[[d]][c]))/
               (as.vector(mu_JG[Z.c,])*beta_CD[[d]][c]+as.vector(phi_JG[Z.c,])),
             
             ncol = 1)
      
    })
    
    loop.result <- do.call(cbind, loop.result)
    Y_latent[[d]] <- loop.result
  }
  
  ##-- Create an empty matrix
  matrix.output <- matrix(0,
                          ncol = 0,
                          nrow = G)
  
  unique.Z <- sort(unique(unlist(Z_CD)))
  
  index <- NULL
  
  ##------------------------- Reorder by cells ------------------------------
  
  for(clust in unique.Z){
    for(d in 1:D){
      
      matrix.output <- cbind(matrix.output, Y_latent[[d]][,which(Z_CD[[d]] == clust)])
      
      ##-- Update index
      index <- c(index, length(which(Z[[d]] == clust)))
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
  }
  
  
  
  ##---------------------- Reorder by genes ----------------------------------
  
  ##-- If reorder by both genes and cells
  if(is.null(gmg_output)){
    
    ##-- Reorder genes by trail probability
    marker_DE <- arrange(gmg_output$marker_DE, 
                         tail.prob.mu)
    
    marker_DD <- arrange(gmg_output$marker_DD, 
                         tail.prob.phi)
    
    ##-- Reorder by tail probabilities
    Y_latent_DE <- lapply(1:D, 
                          function(d) matrix.output[[d]][marker_DE$gene,])
    
    Y_latent_DD <- lapply(1:D,
                          function(d) matrix.output[[d]][marker_DE$gene,])
    
    ##-- Number of DE and DD genes
    DE_number <- length(which(marker_DE$tail.prob.mu >= gmg_output$alpha_M))
    DD_number <- length(which(marker_DD$tail.prob.phi >= gmg_output$alpha_D))
    
    ##-- Return
    return(list('Y_latent_DE' = Y_latent_DE,
                'Y_latent_DD' = Y_latent_DD,
                'DE_number' = DE_number,
                'DD_number' = DD_number,
                'index.solid' = index.solid,
                'index.dashed' = index.dashed))
  }
  
}

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
      index <- c(index, length(which(Z[[d]] == clust)))
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
  }
  
  
  
  ##---------------------- Reorder by genes ----------------------------------
  
  ##-- If reorder by both genes and cells
  if(is.null(gmg_output)){
    
    ##-- Reorder genes by trail probability
    marker_DE <- arrange(gmg_output$marker_DE, 
                         tail.prob.mu)
    
    marker_DD <- arrange(gmg_output$marker_DD, 
                         tail.prob.phi)
    
    ##-- Reorder by tail probabilities
    Y_DE <- lapply(1:D, 
                    function(d) matrix.output[[d]][marker_DE$gene,])
    
    Y_DD <- lapply(1:D,
                   function(d) matrix.output[[d]][marker_DE$gene,])
    
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



