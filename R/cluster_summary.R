

# Cluster summary plot
cluster_summary <- function(Z,
                            data.name){

  ##-- Unique clusters
  Z.unique <- unique(unlist(Z))

  ##-- Number of cells in each cluster
  loop.result <- unlist(lapply(Z.unique, function(j) length(which(unlist(Z) == j))))

  number_of_cells <- data.frame(cluster = as.factor(Z.unique),
                                counts = loop.result)

  ##-- Plot
  plot1 <- ggplot(number_of_cells, aes(x=cluster, y=counts))+
    geom_bar(stat = 'identity')+
    theme_bw()+
    ggtitle('Number of cells in each cluster')+
    ylab('frequency')

  # If length is 2
  if(length(Z) == 2){

    ##-- Plot the proportion of data 1 and data 2 cells in each cluster
    loop.result <- lapply(Z.unique, function(j){

      data.frame(cluster = j,
                 size_data1 = length(which(Z[[1]] == j)),
                 size_data2 = length(which(Z[[2]] == j)))
    })

    number_of_cells_per_dataset <- do.call(rbind, loop.result)
    number_of_cells_per_dataset$cluster <- as.factor(number_of_cells_per_dataset$cluster)

    ##-- Pivot wider
    summary_cluster_2 <- data.frame(cluster = rep(number_of_cells_per_dataset$cluster,2),

                                    dataset = c(rep(data.name[1],nrow(number_of_cells_per_dataset)),
                                                rep(data.name[2],nrow(number_of_cells_per_dataset))),

                                    count = c(number_of_cells_per_dataset$size_data1,
                                              number_of_cells_per_dataset$size_data2))

    ##-- Convert counts into percentages
    summary_cluster_2 <- summary_cluster_2 %>%
      group_by(cluster) %>%
      mutate(percentage_per_cluster = count/sum(count))

    ##-- Overall percentage of HOM cells
    data2.prop.overall <- length(Z[[2]])/(length(unlist(Z)))
    data1.prop.overall <- 1-data2.prop.overall

    ##-- Plot to show composition of each cluster
    plot2 <- ggplot(summary_cluster_2, aes(x=cluster, y=percentage_per_cluster, fill=dataset))+
      geom_bar(stat = 'identity')+
      ylab('Percentage')+
      ggtitle(paste0('Percentage of ', data.name[1], '/ ', data.name[2], ' cells in each cluster'))+
      theme_bw()+
      geom_hline(yintercept = data2.prop.overall)

    ggarrange(plot1, plot2, nrow = 1)

  }else{

    plot1
  }



}

