#' Create Orthogonal Polynomials
#'
#' Generates orthogonal polynomials up to a specified degree and numerically
#' approximates their coefficients.
#'
#' @param degree Integer. The maximum degree of the polynomial.
#' @param X Numeric vector. The exposure values.
#' @param force_x_zero Logical. If TRUE (default), forces the constant term to zero.
#'
#' @return A list containing:
#'   \item{values}{Matrix of orthogonal polynomial values}
#'   \item{coef}{List of coefficient vectors for each polynomial degree}
#'
#' @export
#' @examples
#' \dontrun{
#' X <- rnorm(100)
#' poly_result <- orthopol(degree = 3, X = X)
#' }
orthopol <- function(degree, X, force_x_zero = TRUE) {
  polyX <- poly(X, degree = degree, raw = FALSE)

  coef_list <- list()
  def_df <- c()

  for (i in 1:degree) {
    def_df <- dplyr::bind_cols(def_df, data.frame(tmp = X^i)) %>%
      dplyr::rename(!!paste0("X", i) := tmp)
    mod <- lm(polyX[, i] ~ as.matrix(def_df))
    coef_list[[i]] <- mod$coefficients
  }

  if (force_x_zero == TRUE) {
    for (i in 1:degree) {
      polyX[, i] <- polyX[, i] - coef_list[[i]][1]
    }
  }

  return(list(values = polyX, coef = coef_list))
}
