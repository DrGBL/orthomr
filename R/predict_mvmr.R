#' Predict Using MVMR Polynomial Coefficients
#'
#' Given exposure values and polynomial coefficients from MVMR, predicts
#' the outcome values.
#'
#' @param exposure_values Numeric vector. Exposure values.
#' @param polynomial_coefficients List. Polynomial coefficients from \code{obtain_final_coeffs()}.
#'
#' @return Numeric vector of predicted outcome values.
#'
#' @export
#' @examples
#' \dontrun{
#' exposure_new <- seq(-3, 3, 0.1)
#' outcome_predicted <- predict_mvmr(exposure_new, final_coefficients)
#' }
predict_mvmr <- function(exposure_values, polynomial_coefficients) {
  coefficients_working <- polynomial_coefficients

  if ("Intercept" %in% names(polynomial_coefficients)) {
    coefficients_working[["Intercept"]] <- NULL
  }

  predicted_outcome <- rep(0, length(exposure_values))
  for (i in 1:length(coefficients_working)) {
    for (j in 1:length(coefficients_working[[i]])) {
      if (!is.na(coefficients_working[[i]][j])) {
        predicted_outcome <- predicted_outcome + coefficients_working[[i]][j] * exposure_values^(j - 1)
      }
    }
  }

  if ("Intercept" %in% names(polynomial_coefficients)) {
    predicted_outcome <- predicted_outcome + polynomial_coefficients[["Intercept"]]
  }

  return(predicted_outcome)
}
