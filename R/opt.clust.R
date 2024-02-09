# Function 1: Obtain point estimate of clustering based on trace of z and posterior similarity matrix
opt.clust <- function(Z_trace,
                      psm_output){

  ##-- Dimension
  D <- length(Z_trace[[1]])

  C <- sapply(1:D,
              function(d) length(Z_trace[[1]][[d]]))


  ##-- Cumulative sum
  C_cusum <- c(0,cumsum(C))

  ##-- Sample length
  sample.length <- length(Z_trace)

  ##-- Combine all allocations
  Z_CD_sample_tot <- lapply(1:sample.length,
                            function(i) unlist(Z_trace[[i]]))

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
