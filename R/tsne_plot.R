
# UMAP plot
tsne_plot <- function(normHDP_post_output,
                      Y_latent,
                      Y){


  # Dimensions
  Z <- normHDP_post_output$Z


  # Combines latent Y
  Y_latent_cbind <- do.call(cbind, Y_latent)

  # umap for latent counts
  umap_Y_latent <- Rtsne::Rtsne(as.matrix(t(Y_latent_cbind)))

  # df for plotting
  umap_Y_latent_df <- data.frame(tsne1 = umap_Y_latent$Y[,1],
                                 tsne2 = umap_Y_latent$Y[,2],
                                 cluster = as.factor(unlist(Z)))

  # umap for the latent counts
  plot1 <- umap_Y_latent_df %>%
    ggplot()+
    geom_point(mapping = aes(x = tsne1,
                             y = tsne2,
                             colour = cluster),
               size = 0.1)+
    theme_bw()+
    ggtitle('Estimated latent counts')



  # umap for observed counts
  umap_Y <- Rtsne::Rtsne(as.matrix(t(do.call(cbind, Y))))

  umap_Y_df <- data.frame(tsne1 = umap_Y$Y[,1],
                          tsne2 = umap_Y$Y[,2],
                          cluster = as.factor(unlist(Z)))


  plot2 <- umap_Y_df %>%
    ggplot()+
    geom_point(mapping = aes(x = tsne1,
                             y = tsne2,
                             colour = cluster),
               size = 0.1)+
    theme_bw()+
    ggtitle('Observed counts')

  ggarrange(plot1, plot2)

}

tsne_global_marker_genes <- function(gmg_output,
                                     Y,
                                     normHDP_post_output,
                                     Y_latent){

  # Select the index of global marker genes
  global.DE.index <- gmg_output$marker_DE %>%
    filter(class.mu == 'DE') %>%
    select(gene) %>%
    pull()

  global.DD.index <- gmg_output$marker_DD %>%
    filter(class.phi == 'DD') %>%
    select(gene) %>%
    pull()

  global.marker.index <- intersect(global.DE.index,
                                   global.DD.index)

  Z <- normHDP_post_output$Z
  D <- length(Y)

  Y_latent_marker <- lapply(1:D,
                            function(d) Y_latent[[d]][global.marker.index,])


  # umap for latent counts
  umap_Y_latent <- Rtsne::Rtsne(as.matrix(t(do.call(cbind, Y_latent_marker))))

  umap_Y_latent_df <- data.frame(tsne1 = umap_Y_latent$Y[,1],
                                 tsne2 = umap_Y_latent$Y[,2],
                                 cluster = as.factor(unlist(Z)))

  plot1 <- umap_Y_latent_df %>%
    ggplot()+
    geom_point(mapping = aes(x = tsne1,
                             y = tsne2,
                             colour = cluster),
               size = 0.1)+
    theme_bw()+
    ggtitle('estimated latent counts')

  # umap for the observed counts
  Y_marker <- lapply(1:D,
                     function(d){

                       Y[[d]][global.marker.index,]
                     })

  umap_Y <- Rtsne::Rtsne(as.matrix(t(do.call(cbind, Y_marker))))


  umap_Y_df <- data.frame(tsne1 = umap_Y$Y[,1],
                          tsne2 = umap_Y$Y[,2],
                          cluster = as.factor(unlist(Z)))

  plot2 <- umap_Y_df %>%
    ggplot()+
    geom_point(mapping = aes(x = tsne1,
                             y = tsne2,
                             colour = cluster),
               size = 0.1)+
    theme_bw()+
    ggtitle('observed counts')

  ggarrange(plot1, plot2)

}

