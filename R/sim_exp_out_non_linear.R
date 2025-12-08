#' Simulate Exposure and Outcome with Non-Linear Relationship
#'
#' Simulates genetic variants, exposure, and outcome data with a non-linear
#' relationship between exposure and outcome, including confounders and
#' measurement error.
#'
#' @param sample_size Integer. Sample size (default: 10000).
#' @param n_instruments Integer. Number of genetic instruments (default: 30).
#' @param allele_freq_shape1 Numeric. Beta distribution shape parameter 1 for allele
#'   frequencies (default: 0.1).
#' @param allele_freq_shape2 Numeric. Beta distribution shape parameter 2 for allele
#'   frequencies (default: 0.3).
#' @param prop_dominant_instruments Numeric. Proportion of dominant instruments (default: 0.1).
#' @param prop_nonlinear_instruments Numeric. Proportion of non-linear instruments (default: 0.2).
#' @param nonlinear_coding_value Numeric. Value for non-linear instrument coding (default: -1).
#' @param confounder_to_exposure_func Function. Polynomial function relating confounder to exposure
#'   (default: polynomial with coef c(0,9,-6)).
#' @param confounder_to_outcome_func Function. Polynomial function relating confounder to outcome
#'   (default: polynomial with coef c(0,-2,4)).
#' @param confounder_mean Numeric. Confounder mean (default: 0).
#' @param confounder_sd Numeric. Confounder standard deviation (default: 1).
#' @param prop_measurement_error Numeric. Proportion of exposure variance from measurement error
#'   (default: 0.2).
#' @param prop_variance_explained Numeric. Proportion of outcome variance explained by exposure (default: 0.5).
#' @param polynomial_coef_unused Numeric vector. Legacy parameter, currently unused.
#' @param include_pleiotropy Logical. Whether to include pleiotropic effects (default: FALSE).
#' @param prop_pleiotropic Numeric. Proportion of variants with pleiotropic effects
#'   (default: 0.25).
#' @param exposure_to_outcome_func Function. Polynomial function relating exposure to outcome
#'   (default: polynomial with coef c(1,0,5,2)).
#'
#' @return A list containing:
#'   \item{simulated_data}{Data frame with exposure_true, exposure_observed, outcome_true, outcome_observed}
#'   \item{genotypes}{Matrix of genotypes}
#'   \item{allele_frequencies}{Vector of allele frequencies}
#'
#' @export
#' @examples
#' \dontrun{
#' library(polynom)
#' sim_result <- sim_exp_out_non_linear(
#'   sample_size = 5000,
#'   exposure_to_outcome_func = as.function(polynomial(coef = c(1, 10, 3, 0)))
#' )
#' }
sim_exp_out_non_linear <- function(sample_size = 10000,
                                   n_instruments = 30,
                                   allele_freq_shape1 = 0.1,
                                   allele_freq_shape2 = 0.3,
                                   prop_dominant_instruments = 0.1,
                                   prop_nonlinear_instruments = 0.2,
                                   nonlinear_coding_value = -1,
                                   confounder_to_exposure_func = as.function(polynom::polynomial(coef = c(0, 9, -6))),
                                   confounder_to_outcome_func = as.function(polynom::polynomial(coef = c(0, -2, 4))),
                                   confounder_mean = 0,
                                   confounder_sd = 1,
                                   prop_measurement_error = 0.2,
                                   prop_variance_explained = 0.5,
                                   polynomial_coef_unused = c(10, 100),
                                   include_pleiotropy = FALSE,
                                   prop_pleiotropic = 0.25,
                                   exposure_to_outcome_func = as.function(polynom::polynomial(coef = c(1, 0, 5, 2)))) {

  # Allele frequencies of each instrument
  allele_frequencies <- rbeta(n_instruments, 1, 3)

  # Simulate genotypes
  n_dominant <- floor(n_instruments * prop_dominant_instruments)
  n_nonlinear <- floor(prop_nonlinear_instruments * (n_instruments - n_dominant))

  genotypes_additive_all <- matrix(
    rbinom(sample_size * (n_instruments - n_dominant), 2, allele_frequencies[1:(n_instruments - n_dominant)]),
    nrow = sample_size,
    ncol = (n_instruments - n_dominant),
    byrow = TRUE
  )

  genotypes_nonlinear <- genotypes_additive_all[, 1:n_nonlinear]
  genotypes_additive <- genotypes_additive_all[, (n_nonlinear + 1):(n_instruments - n_dominant)]
  genotypes_additive[which(genotypes_additive == 2)] <- nonlinear_coding_value
  genotypes_combined_additive <- cbind(genotypes_nonlinear, genotypes_additive)

  genotypes_dominant <- matrix(
    rbinom(sample_size * n_dominant, 2, allele_frequencies[(n_instruments - n_dominant + 1):n_instruments]),
    nrow = sample_size,
    ncol = n_dominant,
    byrow = TRUE
  )
  genotypes_dominant[which(genotypes_dominant == 1)] <- 0
  genotypes_dominant[which(genotypes_dominant == 2)] <- 1

  genotypes <- cbind(genotypes_combined_additive, genotypes_dominant)

  # Confounder
  confounder <- rnorm(1, confounder_mean, confounder_sd)

  # Genetic effects
  genetic_effects <- rnorm(
    n_instruments,
    mean = 0,
    sd = sqrt((allele_frequencies * (1 - allele_frequencies))^(-0.25))
  )

  # Generate true exposure
  exposure_true <- genotypes %*% genetic_effects + confounder_to_exposure_func(confounder)
  exposure_true <- scale(exposure_true)

  # Add measurement error to exposure
  total_variance_exposure <- var(exposure_true)
  exposure_observed <- exposure_true + rnorm(
    nrow(genotypes),
    0,
    sqrt(total_variance_exposure - (1 - prop_measurement_error) * total_variance_exposure)
  )

  # Generate true outcome (from true exposure, not observed!)
  outcome_true <- exposure_to_outcome_func(exposure_true)
  outcome_true <- scale(outcome_true,center=TRUE,scale=FALSE)

  # Add error to outcome
  total_variance_outcome <- var(outcome_true)
  outcome_observed <- outcome_true +
    rnorm(nrow(genotypes), 0, sqrt(total_variance_outcome - prop_variance_explained * total_variance_outcome)) +
    confounder_to_outcome_func(confounder)

  if (include_pleiotropy == TRUE) {
    n_pleiotropic <- floor(ncol(genotypes) * prop_pleiotropic)
    pleiotropic_effects <- rnorm(
      n_pleiotropic,
      mean = 0,
      sd = sqrt((allele_frequencies[1:n_pleiotropic] * (1 - allele_frequencies[1:n_pleiotropic]))^(-0.25))
    )
    outcome_observed <- outcome_observed + genotypes[, 1:n_pleiotropic] %*% pleiotropic_effects
  }

  return(list(
    simulated_data = data.frame(
      exposure_true = exposure_true,
      exposure_observed = exposure_observed,
      outcome_true = outcome_true,
      outcome_observed = outcome_observed
    ),
    genotypes = cbind(genotypes_additive_all, genotypes_dominant),
    allele_frequencies = allele_frequencies
  ))
}
