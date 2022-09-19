#' Optimal clustering
#'
#' @param normHDP_output Output from function normHDP_mcmc.
#' @param psm_output Output from function similarity_matrix.
#' @return A point estimate of allocations of cells; a list of length D, with each item having a length C[d].
#' @export
opt.clust <- function(normHDP_output,
                      psm_output){

  ##-- Dimension
  C <- normHDP_output$C
  D <- normHDP_output$D

  ##-- Cumulative sum
  C_cusum <- c(0,cumsum(C))

  ##-- Sample length
  sample.length <- length(normHDP_output$b_output)

  ##-- Combine all allocations
  Z_CD_sample_tot <- lapply(1:sample.length, function(i) unlist(normHDP_output$Z_output[[i]]))

  ##-- Convert list to a matrix
  Z_CD_sample_tot <- do.call(rbind, Z_CD_sample_tot)

  ##-- Optimal clustering
  Z_point_estimate <- minVI(psm = psm_output$psm.combined,
                            cls.draw = Z_CD_sample_tot,
                            method = 'all')
  Z_point_estimate <- Z_point_estimate$cl[1,]

  Z_point_estimate_D <- NULL
  for(d in 1:D) Z_point_estimate_D[[d]] <- Z_point_estimate[(C_cusum[d]+1):C_cusum[d+1]]

  return(Z_point_estimate_D)
}



