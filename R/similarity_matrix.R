
#' Simlarity between cells
#'
#' @param Z_trace Allocation samples.
#' @return The output contains 2 items: the first item is the posterior similarity matrix between 2 datasets or within the same dataset.
#' For example: psm.within[[i]][[j]] is the similarity between cells in dataset i and j. Note that to save memory space, i is never greater
#' than j to avoid repeats. The second item is the similarity matrix to compare cells in all datasets.
#' @export
similarity_matrix <- function(Z_trace,
                              num.cores,
                              run.on.pc = TRUE){


  ##-- Dimension
  D <- length(Z_trace[[1]])

  C <- sapply(1:D,
              function(d) length(Z_trace[[1]][[d]]))

  ##-- Cumulative sums of cells
  C_cum <- c(0, cumsum(C))

  ##-------------------- Similarity between cells from the same dataset and across datasets --------------

  psm.within <- NULL
  for(d in 1:D) psm.within[[d]] <- NULL

  ##-- For all possible combinations of datasets
  for(d1 in rev(1:D)){
    for(d2 in 1:d1){

      # Register Cores
      if(run.on.pc == FALSE){

        cl <- makeCluster(num.cores,
                          type = "FORK")
      }else{

        cl <- makeCluster(num.cores)
      }

      loop.result <- pblapply(1:length(Z_trace),
                              cl = cl,
                              FUN = function(iter){

                                Z_trace_d1_iter <- Z_trace[[iter]][[d1]]
                                Z_trace_d2_iter <- Z_trace[[iter]][[d2]]

                                psm.empty <- matrix(0, nrow = C[d1], ncol = C[d2])
                                for(i in 1:C[d1]){
                                  psm.empty[i,] <- psm.empty[i,] + ifelse(Z_trace_d2_iter == Z_trace_d1_iter[i],
                                                                          1,
                                                                          0)
                                }

                                psm.empty

                              })

      # Stop parallel computing
      stopCluster(cl)


      # Sum
      loop.result <- Reduce('+', loop.result)
      psm.within[[d1]][[d2]] <- loop.result/length(Z_trace)

    }
  }

  ##-------------------------------- Computing a combined matrix -------------------------

  combined_matrix <- matrix(0, nrow = max(C_cum), ncol = max(C_cum))

  for(d1 in 1:D){
    for(d2 in 1:d1){

      combined_matrix[(C_cum[d1]+1):C_cum[d1+1], (C_cum[d2]+1):C_cum[d2+1]] <- psm.within[[d1]][[d2]]
    }

    if(d1 != D){
      for(d2 in (d1+1):D){

        combined_matrix[(C_cum[d1]+1):C_cum[d1+1], (C_cum[d2]+1):C_cum[d2+1]] <- t(psm.within[[d2]][[d1]])
      }
    }
  }

  ##-------------------------------- Return both -----------------------------------------

  return(list('psm.within' = psm.within,
              'psm.combined' = combined_matrix))

}
