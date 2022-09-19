
#' Plot for posterior predictive checks with multiple datasets
#'
#' @param ppc_output Output from the function ppc.
#' @param data.name Name of datasets.
#' @return For each dataset, plot mean of log shifted counts against standard deviation of log shifted counts,
#' and plot log of mean shifted counts against dropout probabilities.
#' @export
ppc_plot <- function(ppc_output,
                     data.name){

  ##-- Number of datasets
  D <- length(data.name)

  ##-- Compare relationships
  Y_statistics <- ppc_output$Y_statistics
  rep_Y_statistics <- ppc_output$rep_Y_statistics

  ##-- First plot results for observed data
  t.max <- max(rep_Y_statistics$t)

  par(mfrow=c(1,2*D))

  ##-- Relationship between mean and standard deviation of log shifted counts
  for(d in 1:D){

    ##-- Plot for the first replicated Y - keep as the bottom layer
    loess.Y <- loess(sd.log.shifted.counts ~ mean.log.shifted.counts,
                     data = rep_Y_statistics %>%
                       filter(dataset == d, t == 1),
                     span = 0.50)

    ##-- Fitted Y
    predict.Y <- predict(loess.Y)

    ##-- x
    x <- rep_Y_statistics %>%
      filter(dataset == d, t == 1) %>%
      select(mean.log.shifted.counts) %>%
      pull()

    ##-- ordering x
    x.order <- round(order(x))

    ##-- plot
    plot(x = x[x.order],
         y = predict.Y[x.order],
         xlab = 'Mean of log-shifted counts',
         ylab = 'Standard deviation of log-shifted counts',
         main = data.name[d],
         type = 'l',
         col = 'grey',
         lwd = 0.2,
         ylim = c(0,1),
         cex.lab = 1.5)

    ##-- For replicated Y
    for(t0 in 2:t.max){

      loess.Y <- loess(sd.log.shifted.counts ~ mean.log.shifted.counts,
                       data = rep_Y_statistics %>%
                         filter(dataset == d, t == t0),
                       span = 0.50)

      x <- rep_Y_statistics %>%
        filter(dataset == d, t == t0) %>%
        select(mean.log.shifted.counts) %>%
        pull()

      predict.Y <- predict(loess.Y)

      x.order <- round(order(x))

      lines(x[x.order],
            predict.Y[x.order],
            col = 'grey',
            lwd = 0.2)
    }

    ##-- For observed Y
    loess.Y <- loess(sd.log.shifted.counts ~ mean.log.shifted.counts,
                     data = Y_statistics %>%
                       filter(dataset == d),
                     span = 0.50)

    x <- Y_statistics %>%
      filter(dataset == d) %>%
      select(mean.log.shifted.counts) %>%
      pull()

    predict.Y <- predict(loess.Y)

    x.order <- round(order(x))

    lines(x[x.order],
          predict.Y[x.order],
          col = 'red',
          lwd = 3)


  }

  ##-- Relationship between log of mean counts and dropout probabilities
  for(d in 1:D){

    ##-- For replicated Y: number 1

    loess.Y <- loess(dropout.probability ~ log.mean.counts,
                     data = rep_Y_statistics %>%
                       filter(dataset == d, t == 1, !is.infinite(log.mean.counts)),
                     span = 0.50)

    x <- rep_Y_statistics %>%
      filter(dataset == d, t == 1, !is.infinite(log.mean.counts)) %>%
      select(log.mean.counts) %>%
      pull()

    x.order <- round(order(x))
    predict.Y <- predict(loess.Y)

    plot(x = x[x.order],
         y = predict.Y[x.order],
         xlab = 'log of mean gene counts',
         ylab = 'Dropout probabilities',
         main = data.name[d],
         type = 'l',
         col = 'grey',
         ylim = c(0,1),
         cex.lab = 1.5,
         lwd = 0.2)

    ##-- For replicated Y: number 2:t.max
    for(t0 in 2:t.max){

      loess.Y <- loess(dropout.probability ~ log.mean.counts,
                       data = rep_Y_statistics %>%
                         filter(dataset == d, t == t0, !is.infinite(log.mean.counts)),
                       span = 0.50)

      x <- rep_Y_statistics %>%
        filter(dataset == d, t == t0, !is.infinite(log.mean.counts)) %>%
        select(log.mean.counts) %>%
        pull()

      x.order <- round(order(x))
      predict.Y <- predict(loess.Y)

      lines(x[x.order],
            predict.Y[x.order],
            col = 'grey',
            lwd = 0.2)
    }


    ##-- For observed Y
    loess.Y <- loess(dropout.probability ~ log.mean.counts,
                     data = Y_statistics %>%
                       filter(dataset == d, !is.infinite(log.mean.counts)),
                     span = 0.50)

    x <- Y_statistics %>%
      filter(dataset == d, !is.infinite(log.mean.counts)) %>%
      select(log.mean.counts) %>%
      pull()

    x.order <- round(order(x))
    predict.Y <- predict(loess.Y)

    lines(x[x.order],
          predict.Y[x.order],
          col = 'red',
          lwd = 3)
  }
}

