#' Harmonize and Run Complete MR Analysis
#'
#' Wrapper function that performs the complete MR workflow: harmonization,
#' degree selection, MVMR analysis, and visualization.
#'
#' @param mvmr_exp Data frame. MVMR exposure data in TwoSampleMR format.
#' @param mvmr_out Data frame. MVMR outcome data in TwoSampleMR format.
#' @param test List. Simulation output from \code{sim_exp_out_non_linear()}
#'   (optional, for simulation studies).
#' @param X Numeric vector. Exposure values for prediction (optional, for real data).
#' @param inner_int Numeric. Inner interval proportion for plotting (default: 0.95).
#' @param intercept_choice Logical. Whether to include intercept in MVMR (default: TRUE).
#' @param polyX_coef List. Polynomial coefficients from \code{orthopol()$coef}.
#'
#' @return A list containing:
#'   \item{inst_strength}{Instrument strength statistics}
#'   \item{mvmr_res}{MVMR results}
#'   \item{poly_coeff}{Final polynomial coefficients}
#'   \item{final_plots}{Visualization plots}
#'   \item{pairwise_comparison_full}{Pairwise comparison plots (full range, simulation only)}
#'   \item{pairwise_comparison_inner}{Pairwise comparison plots (inner range, simulation only)}
#'
#' @export
#' @examples
#' \dontrun{
#' # For simulation
#' sim_results <- mr_sim_res(
#'   mvmr_exp = mvmr_exp,
#'   mvmr_out = mvmr_out,
#'   test = sim_data,
#'   polyX_coef = poly_result$coef
#' )
#'
#' # For real data
#' real_results <- mr_sim_res(
#'   mvmr_exp = mvmr_exp,
#'   mvmr_out = mvmr_out,
#'   X = observed_exposure,
#'   polyX_coef = poly_result$coef
#' )
#' }
mr_sim_res <- function(mvmr_exp,
                       mvmr_out,
                       test = NULL,
                       X = NULL,
                       inner_int = 0.95,
                       intercept_choice = TRUE,
                       polyX_coef) {

  gwas_harm_choice <- choose_n_degrees(mvmr_exp, mvmr_out)
  gwas_harm <- gwas_harm_choice[[1]]
  strength_mvmr <- gwas_harm_choice[[2]]

  mvmr_res <- run_mvmr(gwas_harm, intercept = intercept_choice)

  XYcoeff <- obtain_final_coeffs(polyX_coef, mvmr_res, pval = 1, set_higher_to_zero = FALSE)

  if (!is.null(test)) {
    Y <- test[["sims"]]$Y
    X <- test[["sims"]]$X
    Y_true <- test[["sims"]]$Y_true
    X_true <- test[["sims"]]$X_true

    dfres <- data.frame(
      X_true = X_true,
      X = X,
      Y_true = Y_true,
      Y = Y
    )

    quants <- quantile(dfres$X, c((1 - inner_int) / 2, 1 - (1 - inner_int) / 2))

    plot_sim <- plot_mr_res(
      Y_true,
      X_true,
      polyX = polyX_coef,
      xmin = min(dfres$X),
      xmax = max(dfres$X),
      inner_min = quants[1],
      inner_max = quants[2],
      mvmr_res,
      n_draws = 100
    )

    # Error plots without inner_interval
    predY <- predict_mvmr(dfres$X, XYcoeff)
    dfres <- dfres %>%
      dplyr::mutate(Y_predicted = predY) %>%
      dplyr::mutate(error = Y_true - Y_predicted)

    pair_plot_full <- GGally::ggpairs(dfres)

    # Error plots with inner_interval
    dfres <- dfres %>%
      dplyr::filter(X >= quants[1] & X <= quants[2]) %>%
      dplyr::filter(X_true >= quants[1] & X_true <= quants[2])
    predY <- predict_mvmr(dfres$X, XYcoeff)

    dfres <- dfres %>%
      dplyr::mutate(Y_predicted = predY) %>%
      dplyr::mutate(error = Y_true - Y_predicted)

    pair_plot_inner <- GGally::ggpairs(dfres)

    return(list(
      inst_strength = strength_mvmr,
      mvmr_res = mvmr_res,
      poly_coeff = XYcoeff,
      final_plots = plot_sim,
      pairwise_comparison_full = pair_plot_full,
      pairwise_comparison_inner = pair_plot_inner
    ))
  }

  if (!is.null(X)) {
    quants <- quantile(X, c((1 - inner_int) / 2, 1 - (1 - inner_int) / 2))
    plot_sim <- plot_mr_res(
      Y_true = NULL,
      X_true = X,
      polyX = polyX_coef,
      xmin = min(X),
      xmax = max(X),
      inner_min = quants[1],
      inner_max = quants[2],
      mvmr_res = mvmr_res,
      n_draws = 100
    )

    return(list(
      inst_strength = strength_mvmr,
      mvmr_res = mvmr_res,
      poly_coeff = XYcoeff,
      final_plots = plot_sim
    ))
  }

  return(list(
    inst_strength = strength_mvmr,
    mvmr_res = mvmr_res,
    poly_coeff = XYcoeff
  ))
}
