
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
            plot_overlappings)
}
