#' Global marker genes
#' 
#' @param normHDP_post_output Output from the function normHDP_post_mcmc.
#' @param threshold Threshold value for absolute log-fold change. Values are given for both mean expression and dispersion.
#' @param alpha_M Threshold for the tail probabilities corresponding to the mean expressions, genes with tail probabilities
#' greater than alpha_M are classified as global differentially expressed genes. By default, this is set to control the
#' expected false discovery rate to 5 percent.
#' @param alpha_D Threshold for the tail probabilities corresponding to dispersions, genes with tail probabilities
#' greater than alpha_D are classified as global differentially dispersed genes. By default, this is set to control the
#' expected false discovery rate to 5 percent.
#' @return The output contains 4 items: For marker_DE and marker_DD, each is a data frame with 4 variables; gene index,
#' tail probability, absolute log-fold change and a binary variable to state whether the gene is a global marker gene.
#' The final 2 items are the thresholds alpha_M and alpha_D.
#' @export
global_marker_genes <- function(normHDP_post_output, 
                                threshold,
                                alpha_M = NULL,
                                alpha_D = NULL){
  
  ##-- Input from normHDP_post
  Z <- normHDP_post_output$Z
  G <- normHDP_post_output$G
  
  mcmc.mu <- normHDP_post_output$mu_star_1_J_output
  mcmc.phi <- normHDP_post_output$phi_star_1_J_output
  iter.num <- length(mcmc.mu)
  
  ##-- Find all unique clusters
  unique.cluster <- sort(unique(unlist(Z)))
  
  ##-- Output for mu
  output.mu <- data.frame(gene = numeric(),
                          mean.lfc.mu = numeric(),
                          tail.prob.mu = numeric())
  
  output.phi <- data.frame(gene = numeric(),
                           mean.lfc.phi = numeric(),
                           tail.prob.phi = numeric())
  
  
  ##-- All possible cluster combinations
  x <- unique.cluster
  clust.combi <- lapply(2, function(i){
    unique(t(apply(CombSet(x,i,repl = F),1,sort)))
  })[[1]] %>%
    as.data.frame() %>%
    rename(j1 = V1, j2=V2)
  
  ##-- For each gene g
  loop.result <- lapply(1:G, function(g){
    
    ##-- Empty data frame
    compare.clust.mu <- data.frame(j1 = clust.combi[,1],
                                   j2 = clust.combi[,2],
                                   tail.prob.mu = 0,
                                   mean.lfc.mu = 0)
    
    compare.clust.phi <- data.frame(j1 = clust.combi[,1],
                                    j2 = clust.combi[,2],
                                    tail.prob.phi = 0,
                                    mean.lfc.phi = 0)
    
    for(i in 1:nrow(compare.clust.mu)){
      
      current.clust = compare.clust.mu$j1[i]
      other.clust = compare.clust.mu$j2[i]
      
      ##-- Absolute LFC for mu
      lfc.i.DE <- unlist(lapply(1:iter.num, function(iter){
        
        abs(log(mcmc.mu[[iter]][current.clust,g])-log(mcmc.mu[[iter]][other.clust,g]))
      }))
      
      ##-- Absolute LFC for phi
      lfc.i.DD <- unlist(lapply(1:iter.num, function(iter){
        
        abs(log(mcmc.phi[[iter]][current.clust,g])-log(mcmc.phi[[iter]][other.clust,g]))
      }))
      
      ##-- Compute the probability that LFC is greater than the threshold
      compare.clust.mu$tail.prob.mu[i] <- mean(lfc.i.DE >= threshold$mu)
      compare.clust.phi$tail.prob.phi[i] <- mean(lfc.i.DD >= threshold$phi)
      
      ##-- Compute the mean LFC for each pair
      compare.clust.mu$mean.lfc.mu[i] <- mean(lfc.i.DE)
      compare.clust.phi$mean.lfc.phi[i] <- mean(lfc.i.DD)
    }
    
    ##-- Compute the maximum across all combination pairs
    compare.clust.mu.max <- compare.clust.mu[which(compare.clust.mu$tail.prob.mu == max(compare.clust.mu$tail.prob.mu)),]
    compare.clust.phi.max <- compare.clust.phi[which(compare.clust.phi$tail.prob.phi == min(compare.clust.phi$tail.prob.phi)),]
    
    ##-- Compute the final product
    list(data.frame(gene = g,
                    mean.lfc.mu = mean(compare.clust.mu.max$mean.lfc.mu),
                    tail.prob.mu = max(compare.clust.mu$tail.prob.mu)),
         
         data.frame(gene = g,
                    mean.lfc.phi = mean(compare.clust.phi.max$mean.lfc.phi),
                    tail.prob.phi = max(compare.clust.phi$tail.prob.phi)))
  })
  
  loop.result.DE <- lapply(1:G, function(g) loop.result[[g]][[1]])
  loop.result.DD <- lapply(1:G, function(g) loop.result[[g]][[2]])
  
  ##-- Bind by rows
  loop.result.DE <- do.call(rbind, loop.result.DE)
  loop.result.DD <- do.call(rbind, loop.result.DD)
  
  ##-- Update for every unique cluster
  output.mu <- rbind(output.mu, data.frame(gene = 1:G,
                                           mean.lfc.mu = loop.result.DE$mean.lfc.mu,
                                           tail.prob.mu = loop.result.DE$tail.prob.mu))
  
  output.phi <- rbind(output.phi, data.frame(gene = 1:G,
                                             mean.lfc.phi = loop.result.DD$mean.lfc.phi,
                                             tail.prob.phi = loop.result.DD$tail.prob.phi))
  
  ##------------------------- Determine local marker genes --------------------------
  
  ##-------------- Compute threshold to find local marker genes
  if(is.null(alpha_M)){
    
    efdr.test.mu <- seq(from = 0.1, to = 1, by = 0.01)
    
    ##-- Find alpha_M
    loop.result <- lapply(efdr.test.mu, function(alpha_M){
      
      sum(ifelse(output.mu$tail.prob.mu > alpha_M, 1, 0)*(1-output.mu$tail.prob.mu))/sum(1-output.mu$tail.prob.mu)
    })
    
    ##-- Find the alpha_M which controls EFDR to 5 percent
    alpha_M <- efdr.test.mu[which.min(abs(unlist(loop.result) - 0.05))]
  }
  
  
  ##-- Find alpha_D
  if(is.null(alpha_D)){
    
    efdr.test.phi <- seq(from = 0.1, to = 1, by = 0.01)
    
    loop.result <- lapply(efdr.test.phi, function(alpha_D){
      
      sum(ifelse(output.phi$tail.prob.phi > alpha_D, 1, 0)*(1-output.phi$tail.prob.phi))/sum(1-output.phi$tail.prob.phi)
    })
    
    ##-- Find the alpha_M which controls RFDR to 5 percent 
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


#' Plots for global marker genes
#' 
#' @param gmg_output Output from function global_marker_genes.
#' @return The output contains 3 plots. The first plot shows the relationship between absolute log-fold change
#' and tail probabilities corresponding to mean expressions. Genes above the horizontal line are classified as
#' global differentially expressed genes, and vice versa. Plot 2 has the same structure as plot 1, but for dispersions.
#' Plot 3 summarizes the number of genes with single and both global characteristics.
#' @export
global_marker_genes_plot <- function(gmg_output){
  
  ##-- Output from global marker genes function
  marker_DE <- gmg_output$marker_DE
  marker_DD <- gmg_output$marker_DD
  alpha_M <- gmg_output$alpha_M
  alpha_D <- gmg_output$alpha_D
  
  ##-- DE
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
  
  ##-- DD
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
  
  ##-- Overlapping of DE and DD genes
  
  ##-- Combine marker_DE with marker_DD
  marker_DE_DD <- left_join(marker_DE, marker_DD, by = 'gene')
  
  ##-- barchart
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


#' Heat maps of global marker genes
#' 
#' @param gmg_output Output from function global_marker_genes.
#' @param normHDP_post_output Output from function normHDP_post_mcmc.
#' @return The output contains 4 plots. The first plot show the posterior mean of mean expressions with genes reordered
#' by tail probabilities. Genes above the horizontal dashed line are the global differentially expressed genes. Plot 2
#' has the same structure as plot 1, but for dispersions. Plot 3 only contains global DE genes and mean expressions
#' shown in the heat map are normalized to have zero mean across all clusters for each gene. Plot 4 has the same structure as plot 3.
#' @export
global_marker_genes_heatmap <- function(gmg_output, 
                                        normHDP_post_output){
  
  ##-- Output from global marker genes function
  marker_DE <- gmg_output$marker_DE
  marker_DD <- gmg_output$marker_DD
  alpha_M <- gmg_output$alpha_M
  alpha_D <- gmg_output$alpha_D
  
  ##-- Output from normHDP_post function
  mu_star_1_J_output <- normHDP_post_output$mu_star_1_J_output
  phi_star_1_J_output <- normHDP_post_output$phi_star_1_J_output
  
  sample.length <- length(mu_star_1_J_output)
  
  ##-- Compute posterior mean of mu and phi
  mu_posterior_mean <- Reduce('+', mu_star_1_J_output)/sample.length
  phi_posterior_mean <- Reduce('+', phi_star_1_J_output)/sample.length
  
  ##-- A vector of unique clusters
  unique.Z <- sort(unique(unlist(normHDP_post_output$Z)))
  
  ##-- Number of genes
  G <- length(unique(marker_DE$gene))
  
  ##-- Create empty matrix
  output.mu <- matrix(0, nrow = length(unique.Z), ncol = G)
  output.phi <- matrix(0, nrow = length(unique.Z), ncol = G)
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
       at = 1:length(unique.Z),
       labels = paste(1:length(unique.Z)))
  
  
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
       at = 1:length(unique.Z),
       labels = paste(1:length(unique.Z)))
  
  
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
       at = 1:length(unique.Z),
       labels = paste(1:length(unique.Z)))
  
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
       at = 1:length(unique.Z),
       labels = paste(1:length(unique.Z)))
  
}
