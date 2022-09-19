
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
