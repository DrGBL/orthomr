#' Plot Mendelian Randomization Results
#'
#' Creates visualization of MR results with confidence intervals from random draws.
#'
#' @param outcome_true Numeric vector. True outcome values (can be NULL for real data).
#' @param exposure_true Numeric vector. True or observed exposure values.
#' @param orthogonal_poly_coef List. Polynomial coefficients from \code{orthopol()$coefficients}.
#' @param exposure_min Numeric. Minimum exposure value for plotting range.
#' @param exposure_max Numeric. Maximum exposure value for plotting range.
#' @param mvmr_results Data frame. MVMR results from \code{run_mvmr()}.
#' @param n_random_draws Integer. Number of random draws for confidence intervals (default: 100).
#' @param n_plot_steps Integer. Number of steps for prediction range (default: 100).
#' @param inner_quantile_min Numeric. Lower bound for inner interval (optional).
#' @param inner_quantile_max Numeric. Upper bound for inner interval (optional).
#'
#' @return A ggplot object or ggarrange object (if outcome_true is provided).
#'
#' @export
#' @examples
#' \dontrun{
#' plot_mr_res(
#'   outcome_true = sim_result$simulated_data$outcome_true,
#'   exposure_true = sim_result$simulated_data$exposure_true,
#'   orthogonal_poly_coef = poly_result$coefficients,
#'   exposure_min = -3,
#'   exposure_max = 3,
#'   mvmr_results = mvmr_results
#' )
#' }
plot_mr_res <- function(outcome_true,
                        exposure_true,
                        orthogonal_poly_coef,
                        exposure_min,
                        exposure_max,
                        mvmr_results,
                        n_random_draws = 100,
                        n_plot_steps = 100,
                        inner_quantile_min = NULL,
                        inner_quantile_max = NULL) {

  exposure_range <- seq(exposure_min, exposure_max, (exposure_max - exposure_min) / n_plot_steps)

  # Average estimate
  mean_coefficients <- obtain_final_coeffs(orthogonal_poly_coef, mvmr_results, pvalue_threshold = 1)
  prediction_df <- data.frame(
    prediction_group = "Average",
    x = exposure_range,
    y = predict_mvmr(exposure_range, mean_coefficients),
    alpha = 1,
    linewidth = 1
  )

  # Random draws
  for (i in 1:n_random_draws) {
    mvmr_draw <- data.frame(
      exposure = mvmr_results$exposure,
      b = rnorm(nrow(mvmr_results), mean = mvmr_results$b, sd = mvmr_results$se),
      nsnp = mvmr_results$nsnp,
      pval = mvmr_results$pval
    )
    coefficients_draw <- obtain_final_coeffs(orthogonal_poly_coef, mvmr_draw, pvalue_threshold = 1)
    prediction_df <- data.frame(
      prediction_group = paste0("draw", i),
      x = exposure_range,
      y = predict_mvmr(exposure_range, coefficients_draw),
      alpha = 0.25,
      linewidth = 1
    ) %>%
      dplyr::bind_rows(., prediction_df)
  }

  plot_predicted <- prediction_df %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, group = prediction_group)) +
    ggplot2::geom_line(
      data = prediction_df %>% dplyr::filter(prediction_group != "Average"),
      colour = "blue", alpha = 0.125
    ) +
    ggplot2::geom_line(
      data = prediction_df %>% dplyr::filter(prediction_group == "Average"),
      colour = "red", alpha = 1, linewidth = 1.15
    ) +
    ggplot2::theme_bw()

  if (!is.null(outcome_true)) {
    plot_true <- data.frame(X = exposure_true, Y = outcome_true) %>%
      ggplot2::ggplot(ggplot2::aes(X, Y)) +
      ggplot2::geom_line() +
      ggplot2::theme_bw()
  }

  if (!is.null(inner_quantile_min)) {
    plot_predicted <- plot_predicted +
      ggplot2::geom_vline(
        xintercept = inner_quantile_min,
        linetype = "dashed",
        colour = "black"
      )

    if (!is.null(outcome_true)) {
      plot_true <- plot_true +
        ggplot2::geom_vline(
          xintercept = inner_quantile_min,
          linetype = "dashed",
          colour = "black"
        )
    }
  }

  if (!is.null(inner_quantile_max)) {
    plot_predicted <- plot_predicted +
      ggplot2::geom_vline(
        xintercept = inner_quantile_max,
        linetype = "dashed",
        colour = "black"
      )
    if (!is.null(outcome_true)) {
      plot_true <- plot_true +
        ggplot2::geom_vline(
          xintercept = inner_quantile_max,
          linetype = "dashed",
          colour = "black"
        )
    }
  }

  if (!is.null(outcome_true)) {
    combined_plot <- ggpubr::ggarrange(plot_predicted, plot_true, ncol = 2)
  } else {
    combined_plot <- plot_predicted
  }

  return(combined_plot)
}
