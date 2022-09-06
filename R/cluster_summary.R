#' Optimal clustering
#' 
#' @param normHDP_output Output from function normHDP_mcmc.
#' @return A point estimate of allocations of cells; a list of length D, with each item having a length C[d].
#' @export
opt.clust <- function(normHDP_output){
  
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
  Z_point_estimate <- minVI(psm = psm$psm.combined,
                            cls.draw = Z_CD_sample_tot,
                            method = 'all')
  Z_point_estimate <- Z_point_estimate$cl[1,]
  
  Z_point_estimate_D <- NULL
  for(d in 1:D) Z_point_estimate_D[[d]] <- Z_point_estimate[(C_cusum[d]+1):C_cusum[d+1]]
  
  return(Z_point_estimate_D)
}


#' Cluster summary
#' 
#' A summary of cell allocations, note the function below should only be used if there are 2 datasets in total.
#' 
#' @param Z_point A point estimate of cell allocations.
#' @return The output include 2 plots: Plot 1 is a bar chart summarizing the total number of cells allocated
#' to each identified cluster. Plot 2 shows the proportion of cells from the first and second dataset for
#' each cluster, and compare this proportion with the overall true proportion of cells.
#' @export
cluster_summary <- function(Z_point){
  
  ##-- Unique clusters
  Z.unique <- unique(unlist(Z_point))
  
  ##-- Number of cells in each cluster
  loop.result <- unlist(lapply(Z.unique, function(j) length(which(unlist(Z_point) == j))))
  
  number_of_cells <- data.frame(cluster = as.factor(Z.unique),
                                counts = loop.result)
  
  ##-- Plot
  plot1 <- ggplot(number_of_cells, aes(x=cluster, y=counts))+
    geom_bar(stat = 'identity')+
    theme_bw()+
    ggtitle('Number of cells in each cluster')+
    ylab('frequency')
  
  ##-- Plot the proportion of data 1 and data 2 cells in each cluster
  loop.result <- lapply(Z.unique, function(j){
    
    data.frame(cluster = j,
               size_data1 = length(which(Z_point[[1]] == j)),
               size_data2 = length(which(Z_point[[2]] == j)))
  })
  
  number_of_cells_per_dataset <- do.call(rbind, loop.result)
  number_of_cells_per_dataset$cluster <- as.factor(number_of_cells_per_dataset$cluster)
  
  ##-- Pivot wider
  summary_cluster_2 <- data.frame(cluster = rep(number_of_cells_per_dataset$cluster,2),
                                  
                                  dataset = c(rep('data 1',nrow(number_of_cells_per_dataset)), 
                                              rep('data 2',nrow(number_of_cells_per_dataset))),
                                  
                                  count = c(number_of_cells_per_dataset$size_data1, 
                                            number_of_cells_per_dataset$size_data2))
  
  ##-- Convert counts into percentages
  summary_cluster_2 <- summary_cluster_2 %>%
    group_by(cluster) %>%
    mutate(percentage_per_cluster = count/sum(count))
  
  ##-- Overall percentage of HOM cells
  data2.prop.overall <- length(Z_point[[2]])/(length(unlist(Z_point)))
  data1.prop.overall <- 1-data2.prop.overall
  
  ##-- Plot to show composition of each cluster
  plot2 <- ggplot(summary_cluster_2, aes(x=cluster, y=percentage_per_cluster, fill=dataset))+
    geom_bar(stat = 'identity')+
    ylab('Percentage')+
    ggtitle('Percentage of data 1/ data 2 cells in each cluster')+
    theme_bw()+
    geom_hline(yintercept = data2.prop.overall)
  
  ggarrange(plot1, plot2, nrow = 1)
  
}

