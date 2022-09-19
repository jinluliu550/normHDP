
#' Heatmap for local marker genes
#'
#' @param lmg_output Output from function local_marker_genes.
#' @param normHDP_post_output Output from function normHDP_post_mcmc.
#' @return For mean expressions, for each cluster, show a heat-map for locally differentially expressed genes. Genes are reordered by
#' ascending tail probabilities. Graphs for dispersions have the same structure.
#' @export
local_marker_genes_heatmap <- function(lmg_output,
                                       normHDP_post_output){

  ##-- Tail probabilities and mean lfc
  marker_DE <- lmg_output$marker_DE
  marker_DD <- lmg_output$marker_DD

  ##-- Thresholds
  alpha_M <- lmg_output$alpha_M
  alpha_D <- lmg_output$alpha_D

  ##-- Point estimate of clusters
  Z <- normHDP_post_output$Z

  ##-- MCMC samples of mu and phi
  mu_star_1_J_output <- normHDP_post_output$mu_star_1_J_output
  phi_star_1_J_output <- normHDP_post_output$phi_star_1_J_output

  sample.length <- length(mu_star_1_J_output)

  ##-- Compute posterior mean of mu and phi
  mu_posterior_mean <- Reduce('+', mu_star_1_J_output)/sample.length
  phi_posterior_mean <- Reduce('+', phi_star_1_J_output)/sample.length

  ##-- Unique clusters
  unique.Z <- sort(unique(unlist(Z)))

  loop.result <- lapply(unique.Z, function(clust){

    ##-- Find the ones correspond to the target cluster and rearrange by tail probabilities
    marker_DE_clust <- marker_DE %>%
      filter(cluster == clust, tail.prob >= alpha_M) %>%
      arrange(tail.prob)

    marker_DD_clust <- marker_DD %>%
      filter(cluster == clust, tail.prob >= alpha_D) %>%
      arrange(tail.prob)

    ##-- log mu arranged by tail probabilities
    log.DE <- log(mu_posterior_mean)[,marker_DE_clust$gene]
    log.DD <- log(phi_posterior_mean[,marker_DD_clust$gene])

    list(log.DE, log.DD)
  })

  ##-- Extract log mu and log phi from the loop
  log.DE <- lapply(unique.Z, function(clust) loop.result[[clust]][[1]])
  log.DD <- lapply(unique.Z, function(clust) loop.result[[clust]][[2]])

  ##-- Specify dimensions
  col.number <- ceiling(length(unique.Z)/2)
  par(mfrow=c(2, col.number))

  ##-- Plot for each cluster for DE
  for(clust in unique.Z){

    image.plot(1:nrow(log.DE[[clust]]),
               1:ncol(log.DE[[clust]]),
               log.DE[[clust]],
               xlab='clusters',
               ylab='genes',
               legend.width = 0.5,
               main = paste('cluster', clust, sep = ' '),
               cex.lab = 1.5,
               cex.main = 1.5)
  }

  ##-- Plot for each cluster for DD
  par(mfrow=c(2, col.number))

  ##-- Plot for each cluster for DE
  for(clust in unique.Z){

    image.plot(1:nrow(log.DD[[clust]]),
               1:ncol(log.DD[[clust]]),
               log.DD[[clust]],
               xlab='clusters',
               ylab='genes',
               legend.width = 0.5,
               main = paste('cluster', clust, sep = ' '),
               cex.lab = 1.5,
               cex.main = 1.5)
  }
}
