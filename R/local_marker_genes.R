#' Local marker genes
#'
#' Distinguish between local marker and non-marker genes.
#' @param normHDP_post_output Output from function normHDP_mcmc_post.
#' @param threshold Threshold for the log-fold change between mean expressions or dispersions.
#' @param alpha_M Genes with tail probabilities of mean expressions greater than alpha_M are classified as local differentially expressed genes.
#' By default, this is set to control expected false discovery rate to 5 percent.
#' @param alpha_D Genes with tail probabilities of dispersions greater than alpha_D are classified as local differentially dispersed genes.
#' By default, this is set to control expected false discovery rate to 5 percent.
#' @return Output is a list of 4 items. Item 1: marker_DE - a data frame with 5 variables; including cluster, gene, tail probability,
#' log-fold change and a binary variable to indicate whether the gene is a local marker gene. Item 2: marker_DD - same structure as marker_DE.
#' Item 3 and 4 are alpha_M and alpha_D.
#' @export
local_marker_genes <- function(normHDP_post_output,
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

  ##-- Create an empty dataset
  marker.DE.df <- data.frame(gene = numeric(),
                             cluster = numeric(),
                             tail.prob = numeric(),
                             mean.lfc = numeric())

  marker.DD.df <- data.frame(gene = numeric(),
                             cluster = numeric(),
                             tail.prob = numeric(),
                             mean.lfc = numeric())

  ##-- For each of the cluster
  for(clust in unique.cluster){

    loop.result <- lapply(1:G, function(g){

      ##-- Table to store output for DE and DD
      tail.prob.DE <- data.frame(current.clust = clust,
                                 remaining.clust = unique.cluster[-which(unique.cluster==clust)],
                                 tail.prob = 0,
                                 abs_lfc = 0)

      tail.prob.DD <- data.frame(current.clust = clust,
                                 remaining.clust = unique.cluster[-which(unique.cluster==clust)],
                                 tail.prob = 0,
                                 abs_lfc = 0)

      ##-- For each of the cluster pair
      for(i in 1:nrow(tail.prob.DE)){

        current.clust <- clust
        other.clust <- tail.prob.DE$remaining.clust[i]

        ##-- Absolute LFC for mu
        lfc.i.DE <- unlist(lapply(1:iter.num, function(iter){

          abs(log(mcmc.mu[[iter]][current.clust,g])-log(mcmc.mu[[iter]][other.clust,g]))
        }))

        ##-- Absolute LFC for phi
        lfc.i.DD <- unlist(lapply(1:iter.num, function(iter){

          abs(log(mcmc.phi[[iter]][current.clust,g])-log(mcmc.phi[[iter]][other.clust,g]))
        }))

        ##-- Compute the probability that LFC is greater than the threshold
        tail.prob.DE$tail.prob[i] <- mean(lfc.i.DE >= threshold$mu)
        tail.prob.DD$tail.prob[i] <- mean(lfc.i.DD >= threshold$phi)

        ##-- Compute the mean LFC for each pair
        tail.prob.DE$abs_lfc[i] <- mean(lfc.i.DE)
        tail.prob.DD$abs_lfc[i] <- mean(lfc.i.DD)
      }

      ##-- pair which corresponds to the minimum tail probability
      tail.prob.DE.min <- tail.prob.DE[which(tail.prob.DE$tail.prob == min(tail.prob.DE$tail.prob)),]
      tail.prob.DD.min <- tail.prob.DD[which(tail.prob.DD$tail.prob == min(tail.prob.DD$tail.prob)),]

      list(data.frame(gene = g,
                      cluster = clust,
                      tail.prob = min(tail.prob.DE$tail.prob),
                      mean.lfc = mean(tail.prob.DE.min$abs_lfc)),

           data.frame(gene = g,
                      cluster = clust,
                      tail.prob = min(tail.prob.DD$tail.prob),
                      mean.lfc = mean(tail.prob.DD.min$abs_lfc)))

    })

    ##-- lfc for DE and DD
    loop.result.DE <- lapply(1:G, function(g) loop.result[[g]][[1]])
    loop.result.DD <- lapply(1:G, function(g) loop.result[[g]][[2]])

    ##-- Bind by rows
    loop.result.DE <- do.call(rbind, loop.result.DE)
    loop.result.DD <- do.call(rbind, loop.result.DD)

    ##-- Update for every unique cluster
    marker.DE.df <- rbind(marker.DE.df, data.frame(gene = 1:G,
                                                   cluster = clust,
                                                   tail.prob = loop.result.DE$tail.prob,
                                                   mean.lfc = loop.result.DE$mean.lfc))

    marker.DD.df <- rbind(marker.DE.df, data.frame(gene = 1:G,
                                                   cluster = clust,
                                                   tail.prob = loop.result.DD$tail.prob,
                                                   mean.lfc = loop.result.DD$mean.lfc))
  }

  ##-- Convert clust to factor
  marker.DE.df$cluster <- as.factor(marker.DE.df$cluster)
  marker.DD.df$cluster <- as.factor(marker.DD.df$cluster)

  ##------------------------- Determine local marker genes --------------------------

  ##-------------- Compute threshold to find local marker genes
  if(is.null(alpha_M)){

    efdr.test.mu <- seq(from = 0.1, to = 1, by = 0.01)

    ##-- Find alpha_M
    loop.result <- lapply(efdr.test.mu, function(alpha_M){

      sum(ifelse(marker.DE.df$tail.prob > alpha_M, 1, 0)*(1-marker.DE.df$tail.prob))/sum(1-marker.DE.df$tail.prob)
    })

    ##-- Find the alpha_M which controls EFDR to 5 percent
    alpha_M <- efdr.test.mu[which.min(abs(unlist(loop.result) - 0.05))]
  }

  ##-- Find alpha_D
  if(is.null(alpha_D)){

    efdr.test.phi <- seq(from = 0.1, to = 1, by = 0.01)

    loop.result <- lapply(efdr.test.phi, function(alpha_D){

      sum(ifelse(marker.DD.df$tail.prob > alpha_D, 1, 0)*(1-marker.DD.df$tail.prob))/sum(1-marker.DD.df$tail.prob)
    })

    ##-- Find the alpha_M which controls RFDR to 5 percent
    alpha_D <- efdr.test.phi[which.min(abs(unlist(loop.result) - 0.05))]
  }

  ##-- Define class
  marker.DE.df$class <- as.factor(ifelse(marker.DE.df$tail.prob >= alpha_M,
                                         'marker',
                                         'non-marker'))

  marker.DD.df$class <- as.factor(ifelse(marker.DD.df$tail.prob >= alpha_D,
                                         'marker',
                                         'non-marker'))

  ##------------------------------- Return output ------------------------------

  return(list('marker_DE' = marker.DE.df,
              'marker_DD' = marker.DD.df,
              'alpha_M' = alpha_M,
              'alpha_D' = alpha_D))
}

