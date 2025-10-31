#' Plot Mendelian Randomization Results
#'
#' Creates visualization of MR results with confidence intervals from random draws.
#'
#' @param Y_true Numeric vector. True Y values (can be NULL for real data).
#' @param X_true Numeric vector. True or observed X values.
#' @param polyX List. Polynomial coefficients from \code{orthopol()$coef}.
#' @param xmin Numeric. Minimum X value for plotting range.
#' @param xmax Numeric. Maximum X value for plotting range.
#' @param mvmr_res Data frame. MVMR results from \code{run_mvmr()}.
#' @param n_draws Integer. Number of random draws for confidence intervals (default: 100).
#' @param n_steps Integer. Number of steps for prediction range (default: 100).
#' @param inner_min Numeric. Lower bound for inner interval (optional).
#' @param inner_max Numeric. Upper bound for inner interval (optional).
#'
#' @return A ggplot object or ggarrange object (if Y_true is provided).
#'
#' @export
#' @examples
#' \dontrun{
#' plot_mr_res(
#'   Y_true = sim_data$sims$Y_true,
#'   X_true = sim_data$sims$X_true,
#'   polyX = poly_result$coef,
#'   xmin = -3,
#'   xmax = 3,
#'   mvmr_res = mvmr_results
#' )
#' }
plot_mr_res <- function(Y_true,
                        X_true,
                        polyX,
                        xmin,
                        xmax,
                        mvmr_res,
                        n_draws = 100,
                        n_steps = 100,
                        inner_min = NULL,
                        inner_max = NULL) {

  Xrange <- seq(xmin, xmax, (xmax - xmin) / n_steps)

  # Average estimate
  mean_coeffs <- obtain_final_coeffs(polyX, mvmr_res, pval = 1)
  pred_df <- data.frame(
    pred_group = "Average",
    x = Xrange,
    y = predict_mvmr(Xrange, mean_coeffs),
    alpha = 1,
    linewidth = 1
  )

  # Random draws
  for (i in 1:n_draws) {
    mvmr_draw <- data.frame(
      exposure = mvmr_res$exposure,
      b = rnorm(nrow(mvmr_res), mean = mvmr_res$b, sd = mvmr_res$se),
      nsnp = mvmr_res$nsnp,
      pval = mvmr_res$pval
    )
    coeffs_draw <- obtain_final_coeffs(polyX, mvmr_draw, pval = 1)
    pred_df <- data.frame(
      pred_group = paste0("draw", i),
      x = Xrange,
      y = predict_mvmr(Xrange, coeffs_draw),
      alpha = 0.25,
      linewidth = 1
    ) %>%
      dplyr::bind_rows(., pred_df)
  }

  ptmp_1 <- pred_df %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, group = pred_group)) +
    ggplot2::geom_line(
      data = pred_df %>% dplyr::filter(pred_group != "Average"),
      colour = "blue", alpha = 0.125
    ) +
    ggplot2::geom_line(
      data = pred_df %>% dplyr::filter(pred_group == "Average"),
      colour = "red", alpha = 1, linewidth = 1.15
    ) +
    ggplot2::theme_bw()

  if (!is.null(Y_true)) {
    ptmp_2 <- data.frame(X = X_true, Y = Y_true) %>%
      ggplot2::ggplot(ggplot2::aes(X, Y)) +
      ggplot2::geom_line() +
      ggplot2::theme_bw()
  }

  if (!is.null(inner_min)) {
    ptmp_1 <- ptmp_1 +
      ggplot2::geom_vline(
        xintercept = inner_min,
        linetype = "dashed",
        colour = "black"
      )

    if (!is.null(Y_true)) {
      ptmp_2 <- ptmp_2 +
        ggplot2::geom_vline(
          xintercept = inner_min,
          linetype = "dashed",
          colour = "black"
        )
    }
  }

  if (!is.null(inner_max)) {
    ptmp_1 <- ptmp_1 +
      ggplot2::geom_vline(
        xintercept = inner_max,
        linetype = "dashed",
        colour = "black"
      )
    if (!is.null(Y_true)) {
      ptmp_2 <- ptmp_2 +
        ggplot2::geom_vline(
          xintercept = inner_max,
          linetype = "dashed",
          colour = "black"
        )
    }
  }

  if (!is.null(Y_true)) {
    ptmp <- ggpubr::ggarrange(ptmp_1, ptmp_2, ncol = 2)
  } else {
    ptmp <- ptmp_1
  }

  return(ptmp)
}
