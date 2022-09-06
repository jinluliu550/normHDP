#' Plot to show posterior similarity between cells.
#' 
#' @param psm.ind Posterior similarity matrix between 2 datasets or within the same dataset. This is equivalent to the
#' first item of the output from function similarity_matrix.
#' @param psm.tot Posterior similarity matrix between cells in all datasets.
#' @param method Method used to reorder cells, default is to set method to 'complete'.
#' @param ... All other controls for the plot.
#' @return The output includes 2 plots. Plot 1: Similarity between cells within distinguishing between datasets.
#' Plot 2: similarity between cells with cells from different datasets separated by horizontal and vertical lines.
#' @export
plotpsm <- function(psm.ind, 
                    psm.tot,
                    method="complete",...){

  if(any(psm.tot !=t(psm.tot)) | any(psm.tot >1) | any(psm.tot < 0) | sum(diag(psm.tot)) != nrow(psm.tot) ){
    stop("psm.tot must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}

  ##-------- Setting up the graphics
  par(mfrow=c(1,2))

  ##------------------------------ without distinguishing between datasets ------------------------------
  hc=hclust(as.dist(1-psm.tot), method = method, members = NULL)
  psm_hc=psm.tot
  n=nrow(psm.tot)
  psm_hc[1:n,]=psm_hc[hc$order,]
  psm_hc[,1:n]=psm_hc[,hc$order]

  image(1:n,
        1:n,
        1-psm_hc,
        col=heat.colors(20), ...)

  # Dimension of datasets
  D <- length(psm.ind)
  n <- rep(0,D)
  for(d in 1:D){
    n[d] <- nrow(psm.ind[[d]][[d]])
  }
  n <- c(0,cumsum(n))

  ##----------------------------------- Distinguishing between datasets ----------------------------------
  psm_matrix_output <- matrix(0, nrow=max(n), ncol=max(n))
  hc_diag <- NULL

  ##-- For diagonal
  for(d in 1:D){
    hc_diag[[d]] <- hclust(as.dist(1-psm.ind[[d]][[d]]), method=method, members=NULL)
    psm_hc_diag <- psm.ind[[d]][[d]]
    psm_hc_diag[1:nrow(psm_hc_diag),] <- psm_hc_diag[rev(hc_diag[[d]]$order),]
    psm_hc_diag[,1:nrow(psm_hc_diag)] <- psm_hc_diag[,hc_diag[[d]]$order]
    psm_matrix_output[(n[d]+1):n[d+1], (max(n)-n[d+1]+1):(max(n)-n[d])] <- psm_hc_diag
  }

  ##-- For off-diagonal
  for(d1 in 2:D){
    for(d2 in 1:(d1-1)){
      psm_d1_d2 <- t(psm.ind[[d1]][[d2]])
      psm_d1_d2[1:nrow(psm_d1_d2),] <- psm_d1_d2[rev(hc_diag[[d2]]$order),]
      psm_d1_d2[,1:ncol(psm_d1_d2)] <- psm_d1_d2[,hc_diag[[d1]]$order]
      psm_matrix_output[(n[d2]+1):n[d2+1], (max(n)-n[d1+1]+1):(max(n)-n[d1])] <- psm_d1_d2
    }
  }

  for(d1 in 1:(D-1)){
    for(d2 in (d1+1):D){
      psm_d1_d2 <- psm.ind[[d2]][[d1]]
      psm_d1_d2[1:nrow(psm_d1_d2),] <- psm_d1_d2[rev(hc_diag[[d2]]$order),]
      psm_d1_d2[,1:ncol(psm_d1_d2)] <- psm_d1_d2[,hc_diag[[d1]]$order]
      psm_matrix_output[(n[d2]+1):n[d2+1], (max(n)-n[d1+1]+1):(max(n)-n[d1])] <- psm_d1_d2
    }
  }

  n.max <- max(n)
  image(1:n.max, 1:n.max, 1-psm_matrix_output, col=heat.colors(100), ...)
  abline(v=n[-c(1,D+1)], lwd=3)
  abline(h=(max(n)-n)[-c(1,D+1)], lwd=3)
}

