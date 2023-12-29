# Compare the relationship between mu and phi for each group

mu_and_phi_plot <- function(normHDP_post_output,
                            row.n){

  # Dimensions
  J <- nrow(normHDP_post_output$mu_star_1_J_output[[1]])
  G <- ncol(normHDP_post_output$mu_star_1_J_output[[1]])

  alpha_phi2 <- mean(normHDP_post_output$alpha_phi2_output)

  # Posterior mean of mu and phi
  mean_mu <- Reduce("+", normHDP_post_output$mu_star_1_J_output)/length(normHDP_post_output$mu_star_1_J_output)
  mean_phi <- Reduce("+", normHDP_post_output$phi_star_1_J_output)/length(normHDP_post_output$phi_star_1_J_output)

  b_posterior_mean <- colMeans(do.call(rbind, normHDP_post_output$b_output))


  # Sequence of log mu
  log_mu_seq <- seq(from = min(log(mean_mu)),
                    to = max(log(mean_mu)),
                    length.out = 200)

  if(length(b_posterior_mean) == 2){

    log_phi_fit <- b_posterior_mean[1] + b_posterior_mean[2]*log_mu_seq

  }else{

    log_phi_fit <- b_posterior_mean[1] + b_posterior_mean[2]*log_mu_seq + b_posterior_mean[3]*(log_mu_seq)^2
  }



  log_phi_lb <- log_phi_fit - 1.96*sqrt(alpha_phi2)
  log_phi_ub <- log_phi_fit + 1.96*sqrt(alpha_phi2)

  df00 <- data.frame(mu = log_mu_seq,
                     fit.phi = log_phi_fit,
                     fit.phi.lb = log_phi_lb,
                     fit.phi.ub = log_phi_ub)

  df0 <- data.frame(mean_mu = as.vector(mean_mu),
                    mean_phi = as.vector(mean_phi),
                    cluster = as.factor(rep(1:J, G)))

  df0 %>%
    ggplot()+
    geom_point(mapping = aes(x = log(mean_mu),
                             y = log(mean_phi)),
               size = 0.1,
               col = 'grey')+
    facet_wrap(~cluster,
               nrow = row.n)+
    theme_bw()+
    xlab('log mean expressions')+
    ylab('log dispersions')+

    geom_line(data = df00,
              mapping = aes(x = mu,
                            y = fit.phi),
              col = 'red',
              size = 0.1)+

    geom_line(data = df00,
              mapping = aes(x = mu,
                            y = fit.phi.lb),
              col = 'red',
              size = 0.1,
              linetype = 'dashed')+

    geom_line(data = df00,
              mapping = aes(x = mu,
                            y = fit.phi.ub),
              col = 'red',
              size = 0.1,
              linetype = 'dashed')


}
