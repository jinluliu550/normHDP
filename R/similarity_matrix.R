
#' Simlarity between cells
#' 
#' @param normHDP_output Output from the function normHDP_mcmc
#' @return The output contains 2 items: the first item is the posterior similarity matrix between 2 datasets or within the same dataset.
#' For example: psm.within[[i]][[j]] is the similarity between cells in dataset i and j. Note that to save memory space, i is never greater
#' than j to avoid repeats. The second item is the similarity matrix to compare cells in all datasets.
#' @export
similarity_matrix <- function(normHDP_output){
  
  ##-- Input from MCMC
  D <- normHDP_output$D
  C <- normHDP_output$C
  allocation_list <- normHDP_output$Z_output
  
  ##-- Cumulative sums of cells
  C_cum <- c(0, cumsum(C))
  
  ##-------------------- Similarity between cells from the same dataset and across datasets --------------
  
  psm.within <- NULL
  for(d in 1:D) psm.within[[d]] <- NULL
  
  ##-- For all possible combinations of datasets
  for(d1 in rev(1:D)){
    for(d2 in 1:d1){
      
      ##-- For each iteration
      loop.result <- lapply(1:length(allocation_list), function(iter){
        
        ##-- Allocations of cells in dataset d1 and iteration = iter
        allocation_list_d1_iter <- allocation_list[[iter]][[d1]]
        allocation_list_d2_iter <- allocation_list[[iter]][[d2]]
        
        psm.empty <- matrix(0, nrow = C[d1], ncol = C[d2])
        for(i in 1:C[d1]){
          psm.empty[i,] <- psm.empty[i,] + ifelse(allocation_list_d2_iter == allocation_list_d1_iter[i],
                                                  1,
                                                  0)
        }
        
        psm.empty
      })
      
      ##-- Sum up the matrices
      loop.result <- Reduce('+', loop.result)
      psm.within[[d1]][[d2]] <- loop.result/length(allocation_list)
      
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


