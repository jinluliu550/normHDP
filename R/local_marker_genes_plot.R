
#' Plots for local marker genes
#'
#' @param lmg_output Output from function local_marker_genes
#' @return Output contains 4 plots. Plot 1: relationship between absolute log-fold change and tail probabilities for mean
#' expressions with threshold shown. Plot 2: A bar chart to show number of local differentially expressed genes for each cluster. Plot 3 and 4
#' are similar to the structure of the previous 2 plots, but to show relationship and counts for dispersions.
#' @export
local_marker_genes_plot <- function(lmg_output){

  ##-- Input from local marker gene computation
  marker_DE <- lmg_output$marker_DE
  marker_DD <- lmg_output$marker_DD
  alpha_M <- lmg_output$alpha_M
  alpha_D <- lmg_output$alpha_D

  ##-- Plot for DE
  plot_DE_1 <- ggplot(marker_DE) +
    geom_point(mapping = aes(x=mean.lfc, y=tail.prob, color=class), size=0.5)+
    facet_wrap(~cluster, nrow = 2)+
    xlab('Mean absolute log-fold change')+
    ylab('Tail probabilities')+
    theme_bw()+
    geom_hline(yintercept = alpha_M, linetype = "dashed", size = 1)+
    scale_color_manual(values=c('Red', 'Grey'))+
    theme(legend.position = "none")+
    ggtitle('Local Marker Genes (DE)')

  ##-- Barchart
  plot_DE_2 <- marker_DE %>%
    group_by(cluster) %>%
    summarize(marker.count = length(which(class=='marker'))) %>%
    ggplot()+
    geom_bar(mapping = aes(x = cluster, y = marker.count), fill = "#0073C2FF", stat = "identity")+
    theme_pubclean()+
    ylab('Total number of marker genes')+
    xlab('cluster')+
    ggtitle('Local Marker Genes (DE)')

  ##-- Plot for DD
  plot_DD_1 <- ggplot(marker_DD) +
    geom_point(mapping = aes(x=mean.lfc, y=tail.prob, color=class), size=0.5)+
    facet_wrap(~cluster, nrow = 2)+
    xlab('Mean absolute log-fold change')+
    ylab('Tail probabilities')+
    theme_bw()+
    geom_hline(yintercept = alpha_D, linetype = "dashed", size = 1)+
    scale_color_manual(values=c('Red', 'Grey'))+
    theme(legend.position = "none")+
    ggtitle('Local Marker Genes (DD)')

  ##-- Barchart
  plot_DD_2 <- marker_DD %>%
    group_by(cluster) %>%
    summarize(marker.count = length(which(class=='marker'))) %>%
    ggplot()+
    geom_bar(mapping = aes(x = cluster, y = marker.count), fill = "#0073C2FF", stat = "identity")+
    theme_pubclean()+
    ylab('Total number of marker genes')+
    xlab('cluster')+
    ggtitle('Local Marker Genes (DD)')

  ##-- Show plot in R
  print(plot_DE_1)
  print(plot_DE_2)
  print(plot_DD_1)
  print(plot_DD_2)
}

