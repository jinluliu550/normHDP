
# Function 1: Finding the global marker genes
global_marker_genes <- function(normHDP_post_output,
                                threshold,
                                alpha_M = NULL,
                                alpha_D = NULL,
                                num.cores,
                                run.on.pc = TRUE){

  # Dimensions
  Z <- normHDP_post_output$Z
  G <- normHDP_post_output$G
  J <- normHDP_post_output$J

  # Trace and length of unique samples
  mcmc.mu <- normHDP_post_output$mu_star_1_J_output
  mcmc.phi <- normHDP_post_output$phi_star_1_J_output
  iter.num <- length(mcmc.mu)


  # All possible cluster combinations
  x <- 1:J
  clust.combi <- lapply(2, function(i){
    unique(t(apply(CombSet(x,i,repl = F),1,sort)))
  })[[1]] %>%
    as.data.frame() %>%
    rename(j1 = V1, j2=V2)

  # Register Cores
  if(run.on.pc == FALSE){

    cl <- makeCluster(num.cores,
                      type = "FORK")

  }else{

    cl <- makeCluster(num.cores)
  }


  # For each gene
  loop.result <- pblapply(1:G,
                          cl = cl,
                          FUN = function(g){


                            # Dummy data frame
                            compare.clust.mu <- data.frame(j1 = clust.combi[,1],
                                                           j2 = clust.combi[,2],
                                                           tail.prob.mu = 0,
                                                           mean.lfc.mu = 0)

                            compare.clust.phi <- data.frame(j1 = clust.combi[,1],
                                                            j2 = clust.combi[,2],
                                                            tail.prob.phi = 0,
                                                            mean.lfc.phi = 0)

                            # For each possible cluster pair
                            for(i in 1:nrow(compare.clust.mu)){

                              current.clust = compare.clust.mu$j1[i]
                              other.clust = compare.clust.mu$j2[i]

                              # Absolute LFC for mu
                              lfc.i.DE <- unlist(lapply(1:iter.num, function(iter){

                                abs(log(mcmc.mu[[iter]][current.clust,g])-log(mcmc.mu[[iter]][other.clust,g]))
                              }))

                              # Absolute LFC for phi
                              lfc.i.DD <- unlist(lapply(1:iter.num, function(iter){

                                abs(log(mcmc.phi[[iter]][current.clust,g])-log(mcmc.phi[[iter]][other.clust,g]))
                              }))

                              # Compute the probability that LFC is greater than the threshold
                              compare.clust.mu$tail.prob.mu[i] <- mean(lfc.i.DE >= threshold$mu)
                              compare.clust.phi$tail.prob.phi[i] <- mean(lfc.i.DD >= threshold$phi)

                              # Compute the mean LFC for each pair
                              compare.clust.mu$mean.lfc.mu[i] <- mean(lfc.i.DE)
                              compare.clust.phi$mean.lfc.phi[i] <- mean(lfc.i.DD)
                            }


                            # Select the pair of (j,j') with the largest tail probability
                            compare.clust.mu.max <- compare.clust.mu[which(compare.clust.mu$tail.prob.mu == max(compare.clust.mu$tail.prob.mu)),]
                            compare.clust.phi.max <- compare.clust.phi[which(compare.clust.phi$tail.prob.phi == max(compare.clust.phi$tail.prob.phi)),]

                            # Compute mean LFC and tail probability
                            list(data.frame(gene = g,
                                            mean.lfc.mu = mean(compare.clust.mu.max$mean.lfc.mu),
                                            tail.prob.mu = mean(compare.clust.mu.max$tail.prob.mu)),

                                 data.frame(gene = g,
                                            mean.lfc.phi = mean(compare.clust.phi.max$mean.lfc.phi),
                                            tail.prob.phi = mean(compare.clust.phi.max$tail.prob.phi)))

                          })


  stopCluster(cl)

  # Group loop result and bind by rows
  output.mu <- do.call(rbind,

                       lapply(1:G,
                              function(g) loop.result[[g]][[1]]))

  output.phi <- do.call(rbind,

                        lapply(1:G,
                               function(g) loop.result[[g]][[2]]))


  ##------------------------- Threshold --------------------------

  # Find alpha_M
  if(is.null(alpha_M)){

    efdr.test.mu <- seq(from = 0.1, to = 1, by = 0.01)


    loop.result <- lapply(efdr.test.mu, function(alpha_M){

      sum(ifelse(output.mu$tail.prob.mu > alpha_M, 1, 0)*(1-output.mu$tail.prob.mu))/sum(1-output.mu$tail.prob.mu)
    })

    # Find the alpha_M which controls EFDR to 5 percent
    alpha_M <- efdr.test.mu[which.min(abs(unlist(loop.result) - 0.05))]
  }


  # Find alpha_D
  if(is.null(alpha_D)){

    efdr.test.phi <- seq(from = 0.1, to = 1, by = 0.01)

    loop.result <- lapply(efdr.test.phi, function(alpha_D){

      sum(ifelse(output.phi$tail.prob.phi > alpha_D, 1, 0)*(1-output.phi$tail.prob.phi))/sum(1-output.phi$tail.prob.phi)
    })

    # Find the alpha_M which controls RFDR to 5 percent
    alpha_D <- efdr.test.phi[which.min(abs(unlist(loop.result) - 0.05))]
  }


  ##-- Define class
  output.mu$class.mu <- as.factor(ifelse(output.mu$tail.prob.mu >= alpha_M,
                                         'DE',
                                         'non-DE'))

  output.phi$class.phi <- as.factor(ifelse(output.phi$tail.prob.phi >= alpha_D,
                                           'DD',
                                           'non-DD'))

  ##----------------------- Output ------------------------

  return(list('marker_DE' = output.mu,
              'marker_DD' = output.phi,
              'alpha_M' = alpha_M,
              'alpha_D' = alpha_D))
}



# Function 2: Plot relationship between tail probabilities and log-fold change
global_marker_genes_plot <- function(gmg_output){

  marker_DE <- gmg_output$marker_DE
  marker_DD <- gmg_output$marker_DD
  alpha_M <- gmg_output$alpha_M
  alpha_D <- gmg_output$alpha_D


  # DE
  plot_DE <- ggplot(data = marker_DE, mapping = aes(x = mean.lfc.mu, y = tail.prob.mu))+
    geom_point(size=0.5, mapping = aes(color = class.mu))+
    theme_bw()+
    geom_hline(yintercept = alpha_M, linetype = 'dashed', size=1)+
    scale_color_manual(values=c('Red', 'Grey'))+
    theme(legend.position = "none")+
    xlab('Mean absolute log-fold change of mu')+
    ylab('Tail probabilities')+
    ggtitle('Global DE')+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10))

  # DD
  plot_DD <- ggplot(data = marker_DD, mapping = aes(x = mean.lfc.phi, y = tail.prob.phi))+
    geom_point(size=0.5, mapping = aes(color = class.phi))+
    theme_bw()+
    geom_hline(yintercept = alpha_D, linetype = 'dashed', size=1)+
    scale_color_manual(values=c('Red', 'Grey'))+
    theme(legend.position = "none")+
    xlab('Mean absolute log-fold change of phi')+
    ylab('Tail probabilities')+
    ggtitle('Global DD')+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10))


  # Four types of global marker genes
  marker_DE_DD <- left_join(marker_DE, marker_DD, by = 'gene')

  # Summary graph
  plot_overlappings <- marker_DE_DD %>%

    group_by(class.mu, class.phi) %>%
    summarise(freq=n()) %>%
    ungroup() %>%

    ggplot(mapping = aes(x=class.mu, y = freq))+
    geom_bar(stat='identity', position='dodge', mapping = aes(fill = factor(class.phi)))+
    scale_fill_manual(values=c('Red', 'Grey'))+
    ylab('Number of genes')+
    xlab('Change in mean expression')+
    ggtitle('Overlapping Global Marker Genes')+
    theme_bw()+
    labs(fill = "Change in Dispersion")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10))

  ggarrange(plot_DE,
            plot_DD,
            plot_overlappings,
            nrow = 1)
}



# Function 3: Heat-map of global marker genes
global_marker_genes_heatmap <- function(gmg_output,
                                        normHDP_post_output){

  # Output from dependent functions
  marker_DE <- gmg_output$marker_DE
  marker_DD <- gmg_output$marker_DD
  alpha_M <- gmg_output$alpha_M
  alpha_D <- gmg_output$alpha_D

  mu_star_1_J_output <- normHDP_post_output$mu_star_1_J_output
  phi_star_1_J_output <- normHDP_post_output$phi_star_1_J_output

  sample.length <- length(mu_star_1_J_output)

  # Posterior mean
  mu_posterior_mean <- Reduce('+', mu_star_1_J_output)/sample.length
  phi_posterior_mean <- Reduce('+', phi_star_1_J_output)/sample.length

  # List of Unique clusters
  unique.Z <- 1:normHDP_post_output$J

  # Number of genes
  G <- normHDP_post_output$G

  ##-- Create empty matrix
  output.mu <- matrix(0, nrow = normHDP_post_output$J, ncol = G)
  output.phi <- matrix(0, nrow = normHDP_post_output$J, ncol = G)
  rownames(output.mu) <- unique.Z
  rownames(output.phi) <- unique.Z

  ##-- Rank genes by tail probabilities
  marker_DE <- arrange(marker_DE, tail.prob.mu)
  marker_DD <- arrange(marker_DD, tail.prob.phi)

  ##-- Reorder genes by tail probabilities
  for(clust in unique.Z){

    output.mu[clust,] <- mu_posterior_mean[clust,marker_DE$gene]
    output.phi[clust,] <- phi_posterior_mean[clust, marker_DD$gene]
  }

  par(mfrow=c(1,2),
      mar = c(5.1, 4.1, 4.1, 2.1))

  ##-- Plot with all genes for DE
  image.plot(1:nrow(log(output.mu)),
             1:ncol(log(output.mu)),
             log(output.mu),
             xlab='clusters',
             ylab='genes',
             legend.width = 0.5,
             main = 'DE',
             cex.main = 1.5,
             cex.lab = 1.5,
             xaxt = 'n')

  abline(h=G-length(which(marker_DE$tail.prob.mu >= alpha_M)), lwd=3, col = 'red', lty = 2)

  axis(side = 1,
       at = seq(from = 1, to = length(unique.Z), by = 5),
       labels = paste(seq(from = 1, to = length(unique.Z), by = 5)))


  ##-- Plot with all genes for DD
  image.plot(1:nrow(log(output.phi)),
             1:ncol(log(output.phi)),
             log(output.phi),
             xlab='clusters',
             ylab='genes',
             legend.width = 0.5,
             main = 'DD',
             cex.main = 1.5,
             cex.lab = 1.5,
             xaxt = 'n')

  abline(h=G-length(which(marker_DD$tail.prob.phi >= alpha_D)), lwd=3, col = 'red', lty = 2)

  axis(side = 1,
       at = seq(from = 1, to = length(unique.Z), by = 5),
       labels = paste(seq(from = 1, to = length(unique.Z), by = 5)))


  ##-- Plot for global DE and DD genes
  loop.result <- lapply(1:G, function(g){

    ##-- Mean across all clusters for each gene
    mean.g.mu <- mean(log(output.mu[,g]))
    mean.g.phi <- mean(log(output.phi[,g]))

    list(matrix(log(output.mu[,g]) - mean.g.mu,
                ncol = 1),

         matrix(log(output.phi[,g]) - mean.g.phi,
                ncol = 1))
  })

  loop.result.log.mu <- lapply(1:G, function(g) loop.result[[g]][[1]])
  loop.result.log.phi <- lapply(1:G, function(g) loop.result[[g]][[2]])

  log.mu <- do.call(cbind, loop.result.log.mu)
  log.phi <- do.call(cbind, loop.result.log.phi)

  ##-- Find log mu and phi for global marker genes
  ##-- Note the ordering of genes in marker_DE and marker_DD are in ascending order
  log.mu <- log.mu[,(length(which(marker_DE$tail.prob.mu < alpha_M))):nrow(marker_DE)]
  log.phi <- log.phi[,(length(which(marker_DD$tail.prob.phi < alpha_D))):nrow(marker_DD)]

  ##-- Plot for DE
  image.plot(1:nrow(log.mu),
             1:ncol(log.mu),
             log.mu,
             xlab='clusters',
             ylab='genes',
             legend.width = 0.5,
             main = 'DE',
             cex.main = 1.5,
             cex.lab = 1.5,
             xaxt = 'n')

  axis(side = 1,
       at = seq(from = 1, to = length(unique.Z), by = 5),
       labels = paste(seq(from = 1, to = length(unique.Z), by = 5)))

  ##-- Plot for DD
  image.plot(1:nrow(log.phi),
             1:ncol(log.phi),
             log.phi,
             xlab='clusters',
             ylab='genes',
             legend.width = 0.5,
             main = 'DD',
             cex.main = 1.5,
             cex.lab = 1.5)

  axis(side = 1,
       at = seq(from = 1, to = length(unique.Z), by = 5),
       labels = paste(seq(from = 1, to = length(unique.Z), by = 5)))

}

