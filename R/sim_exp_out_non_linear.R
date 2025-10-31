#' Simulate Exposure and Outcome with Non-Linear Relationship
#'
#' Simulates genetic variants, exposure, and outcome data with a non-linear
#' relationship between exposure and outcome, including confounders and
#' measurement error.
#'
#' @param n_s Integer. Sample size (default: 10000).
#' @param n_ins Integer. Number of genetic instruments (default: 30).
#' @param beta_shape1 Numeric. Beta distribution shape parameter 1 for allele
#'   frequencies (default: 0.1).
#' @param beta_shape2 Numeric. Beta distribution shape parameter 2 for allele
#'   frequencies (default: 0.3).
#' @param prop_dom_ins Numeric. Proportion of dominant instruments (default: 0.1).
#' @param prop_non_lin_ins Numeric. Proportion of non-linear instruments (default: 0.2).
#' @param val_non_lin Numeric. Value for non-linear instrument coding (default: -1).
#' @param func_conf_x Function. Polynomial function relating confounder to X
#'   (default: polynomial with coef c(0,9,-6)).
#' @param func_conf_y Function. Polynomial function relating confounder to Y
#'   (default: polynomial with coef c(0,-2,4)).
#' @param conf_mean Numeric. Confounder mean (default: 0).
#' @param conf_sd Numeric. Confounder standard deviation (default: 1).
#' @param prop_var_err Numeric. Proportion of X variance from measurement error
#'   (default: 0.2).
#' @param prop_var_exp_X Numeric. Proportion of Y variance explained by X (default: 0.5).
#' @param a Numeric vector. Coefficients for polynomial from X to Y (currently unused).
#' @param pleiotropy Logical. Whether to include pleiotropic effects (default: FALSE).
#' @param prop_pleio Numeric. Proportion of variants with pleiotropic effects
#'   (default: 0.25).
#' @param func Function. Polynomial function relating X to Y (default: polynomial
#'   with coef c(1,0,5,2)).
#'
#' @return A list containing:
#'   \item{sims}{Data frame with X_true, X, Y_true, Y}
#'   \item{geno}{Matrix of genotypes}
#'   \item{af_ins}{Vector of allele frequencies}
#'
#' @export
#' @examples
#' \dontrun{
#' library(polynom)
#' sim_data <- sim_exp_out_non_linear(
#'   n_s = 5000,
#'   func = as.function(polynomial(coef = c(1, 10, 3, 0)))
#' )
#' }
sim_exp_out_non_linear <- function(n_s = 10000,
                                   n_ins = 30,
                                   beta_shape1 = 0.1,
                                   beta_shape2 = 0.3,
                                   prop_dom_ins = 0.1,
                                   prop_non_lin_ins = 0.2,
                                   val_non_lin = -1,
                                   func_conf_x = as.function(polynom::polynomial(coef = c(0, 9, -6))),
                                   func_conf_y = as.function(polynom::polynomial(coef = c(0, -2, 4))),
                                   conf_mean = 0,
                                   conf_sd = 1,
                                   prop_var_err = 0.2,
                                   prop_var_exp_X = 0.5,
                                   a = c(10, 100),
                                   pleiotropy = FALSE,
                                   prop_pleio = 0.25,
                                   func = as.function(polynom::polynomial(coef = c(1, 0, 5, 2)))) {

  # Allele frequencies of each instrument
  af_ins <- rbeta(n_ins, 1, 3)

  # Simulate genotypes
  n_dom <- floor(n_ins * prop_dom_ins)
  n_non_lin <- floor(prop_non_lin_ins * (n_ins - n_dom))
  G_1_tmp <- matrix(rbinom(n_s * (n_ins - n_dom), 2, af_ins[1:(n_ins - n_dom)]),
                    nrow = n_s, ncol = (n_ins - n_dom), byrow = TRUE)
  G_1_1 <- G_1_tmp[, 1:n_non_lin]
  G_1_2 <- G_1_tmp[, (n_non_lin + 1):(n_ins - n_dom)]
  G_1_2[which(G_1_2 == 2)] <- val_non_lin
  G_1 <- cbind(G_1_1, G_1_2)
  G_2 <- matrix(rbinom(n_s * n_dom, 2, af_ins[(n_ins - n_dom + 1):n_ins]),
                nrow = n_s, ncol = n_dom, byrow = TRUE)
  G_2[which(G_2 == 1)] <- 0
  G_2[which(G_2 == 2)] <- 1
  G <- cbind(G_1, G_2)

  # Confounder
  U <- rnorm(1, conf_mean, conf_sd)

  # Genetic effects
  beta <- rnorm(n_ins, mean = 0, sd = sqrt((af_ins * (1 - af_ins))^(-0.25)))

  # Generate X
  X_true <- G %*% beta + func_conf_x(U)
  X_true <- scale(X_true)

  # Add measurement error to X_true
  tot_var_x <- var(X_true)
  X <- X_true + rnorm(nrow(G), 0, sqrt(tot_var_x - (1 - prop_var_err) * tot_var_x))

  # Generate Y (from X_true, not X!)
  Y_true <- func(X_true)

  # Add error to Y_true
  tot_var_y <- var(Y_true)
  Y <- Y_true + rnorm(nrow(G), 0, sqrt(tot_var_y - prop_var_exp_X * tot_var_y)) + func_conf_y(U)

  if (pleiotropy == TRUE) {
    n_pleio <- floor(ncol(G) * prop_pleio)
    beta_pleio <- rnorm(n_pleio, mean = 0,
                        sd = sqrt((af_ins[1:n_pleio] * (1 - af_ins[1:n_pleio]))^(-0.25)))
    Y <- Y + G[, 1:n_pleio] %*% beta_pleio
  }

  return(list(
    sims = data.frame(X_true = X_true, X = X, Y_true = Y_true, Y = Y),
    geno = cbind(G_1_tmp, G_2),
    af_ins = af_ins
  ))
}
