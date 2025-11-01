#' Create Orthogonal Polynomials
#'
#' Generates orthogonal polynomials up to a specified degree and numerically
#' approximates their coefficients.
#'
#' @param degree Integer. The maximum degree of the polynomial.
#' @param exposure_values Numeric vector. The exposure values.
#' @param force_zero_intercept Logical. If TRUE (default), forces the constant term to zero.
#'
#' @return A list containing:
#'   \item{values}{Matrix of orthogonal polynomial values}
#'   \item{coefficients}{List of coefficient vectors for each polynomial degree}
#'
#' @export
#' @examples
#' \dontrun{
#' exposure_values <- rnorm(100)
#' poly_result <- orthopol(degree = 3, exposure_values = exposure_values)
#' }
orthopol <- function(degree, exposure_values, force_zero_intercept = TRUE) {
  poly_matrix <- poly(exposure_values, degree = degree, raw = FALSE)

  coefficient_list <- list()
  powers_df <- c()

  for (i in 1:degree) {
    powers_df <- dplyr::bind_cols(powers_df, data.frame(tmp = exposure_values^i)) %>%
      dplyr::rename(!!paste0("X", i) := tmp)
    model <- lm(poly_matrix[, i] ~ as.matrix(powers_df))
    coefficient_list[[i]] <- model$coefficients
  }

  if (force_zero_intercept == TRUE) {
    for (i in 1:degree) {
      poly_matrix[, i] <- poly_matrix[, i] - coefficient_list[[i]][1]
    }
  }

  return(list(values = poly_matrix, coefficients = coefficient_list))
}
