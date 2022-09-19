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



