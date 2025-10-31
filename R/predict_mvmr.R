#' Predict Using MVMR Polynomial Coefficients
#'
#' Given exposure values X and polynomial coefficients from MVMR, predicts
#' the outcome Y.
#'
#' @param X Numeric vector. Exposure values.
#' @param XYcoeff List. Polynomial coefficients from \code{obtain_final_coeffs()}.
#'
#' @return Numeric vector of predicted Y values.
#'
#' @export
#' @examples
#' \dontrun{
#' X_new <- seq(-3, 3, 0.1)
#' Y_pred <- predict_mvmr(X_new, final_coef)
#' }
predict_mvmr <- function(X, XYcoeff) {
  XYcoeff_tmp <- XYcoeff

  if ("Intercept" %in% names(XYcoeff)) {
    XYcoeff_tmp[["Intercept"]] <- NULL
  }

  predY <- rep(0, length(X))
  for (i in 1:length(XYcoeff_tmp)) {
    for (j in 1:length(XYcoeff_tmp[[i]])) {
      if (!is.na(XYcoeff_tmp[[i]][j])) {
        predY <- predY + XYcoeff_tmp[[i]][j] * X^(j - 1)
      }
    }
  }

  if ("Intercept" %in% names(XYcoeff)) {
    predY <- predY + XYcoeff[["Intercept"]]
  }

  return(predY)
}
