#' Run Complete MR Analysis Workflow
#'
#' Wrapper function that performs the complete MR workflow: harmonization,
#' degree selection, MVMR analysis, and visualization.
#'
#' @param mvmr_exposure_data Data frame. MVMR exposure data in TwoSampleMR format.
#' @param mvmr_outcome_data Data frame. MVMR outcome data in TwoSampleMR format.
#' @param simulation_data List. Simulation output from \code{sim_exp_out_non_linear()}
#'   (optional, for simulation studies).
#' @param exposure_values Numeric vector. Exposure values for prediction (optional, for real data).
#' @param inner_interval Numeric. Inner interval proportion for plotting (default: 0.95).
#' @param include_intercept Logical. Whether to include intercept in MVMR (default: TRUE).
#' @param orthogonal_poly_coef List. Polynomial coefficients from \code{orthopol()$coefficients}.
#'
#' @return A list containing:
#'   \item{instrument_strength}{Instrument strength statistics}
#'   \item{mvmr_results}{MVMR results}
#'   \item{polynomial_coefficients}{Final polynomial coefficients}
#'   \item{plots}{Visualization plots}
#'   \item{pairwise_comparison_full}{Pairwise comparison plots (full range, simulation only)}
#'   \item{pairwise_comparison_inner}{Pairwise comparison plots (inner range, simulation only)}
#'
#' @export
#' @examples
#' \dontrun{
#' # For simulation
#' sim_results <- mr_sim_res(
#'   mvmr_exposure_data = mvmr_exposure_data,
#'   mvmr_outcome_data = mvmr_outcome_data,
#'   simulation_data = sim_result,
#'   orthogonal_poly_coef = poly_result$coefficients
#' )
#'
#' # For real data
#' real_results <- mr_sim_res(
#'   mvmr_exposure_data = mvmr_exposure_data,
#'   mvmr_outcome_data = mvmr_outcome_data,
#'   exposure_values = observed_exposure,
#'   orthogonal_poly_coef = poly_result$coefficients
#' )
#' }
mr_sim_res <- function(mvmr_exposure_data,
                       mvmr_outcome_data,
                       simulation_data = NULL,
                       exposure_values = NULL,
                       inner_interval = 0.95,
                       include_intercept = TRUE,
                       orthogonal_poly_coef) {

  degree_selection_result <- choose_n_degrees(mvmr_exposure_data, mvmr_outcome_data)
  harmonized_data <- degree_selection_result[[1]]
  instrument_strength <- degree_selection_result[[2]]

  mvmr_results <- run_mvmr(harmonized_data, include_intercept = include_intercept)

  final_coefficients <- obtain_final_coeffs(
    orthogonal_poly_coef,
    mvmr_results,
    pvalue_threshold = 1,
    set_higher_to_zero = FALSE
  )

  if (!is.null(simulation_data)) {
    outcome_observed <- simulation_data[["simulated_data"]]$outcome_observed
    exposure_observed <- simulation_data[["simulated_data"]]$exposure_observed
    outcome_true <- simulation_data[["simulated_data"]]$outcome_true
    exposure_true <- simulation_data[["simulated_data"]]$exposure_true

    results_df <- data.frame(
      exposure_true = exposure_true,
      exposure_observed = exposure_observed,
      outcome_true = outcome_true,
      outcome_observed = outcome_observed
    )

    quantiles <- quantile(
      results_df$exposure_observed,
      c((1 - inner_interval) / 2, 1 - (1 - inner_interval) / 2)
    )

    mr_plot <- plot_mr_res(
      outcome_true,
      exposure_true,
      orthogonal_poly_coef = orthogonal_poly_coef,
      exposure_min = min(results_df$exposure_observed),
      exposure_max = max(results_df$exposure_observed),
      inner_quantile_min = quantiles[1],
      inner_quantile_max = quantiles[2],
      mvmr_results,
      n_random_draws = 100
    )

    # Error plots without inner interval
    predicted_outcome <- predict_mvmr(results_df$exposure_observed, final_coefficients)
    results_df <- results_df %>%
      dplyr::mutate(outcome_predicted = predicted_outcome) %>%
      dplyr::mutate(error = outcome_true - outcome_predicted)

    pairwise_plot_full <- GGally::ggpairs(results_df)

    # Error plots with inner interval
    results_df <- results_df %>%
      dplyr::filter(exposure_observed >= quantiles[1] & exposure_observed <= quantiles[2]) %>%
      dplyr::filter(exposure_true >= quantiles[1] & exposure_true <= quantiles[2])
    predicted_outcome <- predict_mvmr(results_df$exposure_observed, final_coefficients)

    results_df <- results_df %>%
      dplyr::mutate(outcome_predicted = predicted_outcome) %>%
      dplyr::mutate(error = outcome_true - outcome_predicted)

    pairwise_plot_inner <- GGally::ggpairs(results_df)

    return(list(
      instrument_strength = instrument_strength,
      mvmr_results = mvmr_results,
      polynomial_coefficients = final_coefficients,
      plots = mr_plot,
      pairwise_comparison_full = pairwise_plot_full,
      pairwise_comparison_inner = pairwise_plot_inner
    ))
  }

  if (!is.null(exposure_values)) {
    quantiles <- quantile(exposure_values, c((1 - inner_interval) / 2, 1 - (1 - inner_interval) / 2))
    mr_plot <- plot_mr_res(
      outcome_true = NULL,
      exposure_true = exposure_values,
      orthogonal_poly_coef = orthogonal_poly_coef,
      exposure_min = min(exposure_values),
      exposure_max = max(exposure_values),
      inner_quantile_min = quantiles[1],
      inner_quantile_max = quantiles[2],
      mvmr_results = mvmr_results,
      n_random_draws = 100
    )

    return(list(
      instrument_strength = instrument_strength,
      mvmr_results = mvmr_results,
      polynomial_coefficients = final_coefficients,
      plots = mr_plot
    ))
  }

  return(list(
    instrument_strength = instrument_strength,
    mvmr_results = mvmr_results,
    polynomial_coefficients = final_coefficients
  ))
}
