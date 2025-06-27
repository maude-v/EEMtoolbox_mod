#' @title adapted_plot_projections
#' @author Maude Vernet
#' @description
#' Plot projections of species abundance
#' @param projections A data frame containing the projections of species
#' abundance.
#' @param title Title of the plot.
#' @param multiplier A numeric value to scale the abundance values, should be
#' the same as the one used when running `EEM`.
#' @param scaled Logical value indicating whether to scale the y-axis.
#' @return A ggplot object with the projections of species abundance.
#' @export
adapted_plot_projections <- function(projections,
                                     title = NA,
                                     multiplier = 1,
                                     scaled = TRUE) {
  abundance_eq <- dplyr::group_by(projections,
                                  time,
                                  species)
  abundance_eq$pop <- abundance_eq$pop * multiplier
  abundance_eq <- dplyr::summarise(abundance_eq,
                                   median_pop = median(pop),
                                   upper = quantile(pop, 0.975),
                                   lower = quantile(pop, 0.025))

  #basic plot
  if(isTRUE(scaled)) {
  plot <- ggplot2::ggplot(abundance_eq) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Species"),
                    color = ggplot2::guide_legend(title = "Species")) +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Abundance") +
    ggplot2::facet_wrap( ~ species, scales = "free") +
    ggplot2::geom_ribbon(ggplot2::aes(x = time,
                                      ymin = lower,
                                      ymax = upper,
                                      color = species,
                                      fill = species),
                         alpha = 0.1,
                         linewidth = 0.2) +
    ggplot2::geom_line(ggplot2::aes(x = time,
                                    y = median_pop,
                                    color = species),
                       linewidth = 0.8) +
    ggplot2::theme_bw()
  } else {
    plot <- ggplot2::ggplot(abundance_eq) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Species"),
                      color = ggplot2::guide_legend(title = "Species")) +
      ggplot2::xlab("Time") +
      ggplot2::ylab("Abundance") +
      ggplot2::facet_wrap( ~ species) +
      ggplot2::geom_ribbon(ggplot2::aes(x = time,
                                        ymin = lower,
                                        ymax = upper,
                                        color = species,
                                        fill = species),
                           alpha = 0.1,
                           linewidth = 0.2) +
      ggplot2::geom_line(ggplot2::aes(x = time,
                                      y = median_pop,
                                      color = species),
                         linewidth = 0.8) +
      ggplot2::theme_bw()
  }

  if (!is.na(title)) {
    plot <- plot +
      ggplot2::ggtitle(title)
  }

  return(plot)
}

