
heatmap_important_genes <- function(important_genes,
                                    normHDP_post_output,
                                    gene_list){


  # Important genes index
  important_genes_index <- sapply(1:length(important_genes),
                                  function(g) which(gene_list == important_genes[g]))



  # Posterior samples
  sample.length <- length(normHDP_post_output$mu_star_1_J_output)

  mu_JG <- Reduce('+',
                  normHDP_post_output$mu_star_1_J_output)/length(normHDP_post_output$mu_star_1_J_output)

  phi_JG <- Reduce('+',
                   normHDP_post_output$phi_star_1_J_output)/length(normHDP_post_output$phi_star_1_J_output)


  df0 <- data.frame(log.mu = as.vector(log(mu_JG[,important_genes_index])),
                    log.phi = as.vector(log(phi_JG[,important_genes_index])),
                    cluster = rep(1:nrow(mu_JG), length(important_genes_index)),
                    gene.name = rep(important_genes, each = J))

  df0$cluster <- factor(df0$cluster,
                        levels = 1:nrow(mu_JG))

  plot1 <- ggplot(df0, aes(x = cluster,
                           y = gene.name,
                           fill = log.mu))+
    geom_tile(color = 'black')+
    scale_fill_gradient(low = "white", high = "red")

  plot2 <- ggplot(df0, aes(x = cluster,
                           y = gene.name,
                           fill = log.phi))+
    geom_tile(color = 'black')+
    scale_fill_gradient(low = "white", high = "red")

  ggarrange(plot1, plot2)

}

