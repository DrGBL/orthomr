#' Run Multivariable Mendelian Randomization
#'
#' Performs MVMR analysis using harmonized GWAS data.
#'
#' @param harmonized_data List. Harmonized GWAS data from \code{TwoSampleMR::mv_harmonise_data()}
#'   or \code{choose_n_degrees()}.
#' @param include_intercept Logical. Whether to include an intercept term (default: TRUE).
#' @param pvalue_threshold Numeric. P-value threshold for counting significant SNPs (default: 5e-8).
#'
#' @return A data frame with MVMR results containing columns:
#'   \item{exposure}{Exposure name}
#'   \item{b}{Effect estimate}
#'   \item{se}{Standard error}
#'   \item{pval}{P-value}
#'   \item{nsnp}{Number of significant SNPs}
#'
#' @export
#' @examples
#' \dontrun{
#' mvmr_results <- run_mvmr(harmonized_data, include_intercept = TRUE)
#' }
run_mvmr <- function(harmonized_data, include_intercept = TRUE, pvalue_threshold = 5e-8) {

  if (include_intercept == TRUE) {
    model_summary <- summary(lm(
      harmonized_data$outcome_beta ~ harmonized_data$exposure_beta,
      weights = 1 / harmonized_data$outcome_se^(2)
    ))
    mvmr_results <- data.frame(
      exposure = rownames(model_summary$coefficients),
      b = model_summary$coefficients[, 1],
      se = model_summary$coefficients[, 2],
      pval = model_summary$coefficients[, 4],
      nsnp = NA
    ) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "harmonized_data\\$exposure_beta", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\(", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\)", ""))
  } else {
    model_summary <- summary(lm(
      harmonized_data$outcome_beta ~ harmonized_data$exposure_beta - 1,
      weights = 1 / harmonized_data$outcome_se^(2)
    ))
    mvmr_results <- data.frame(
      exposure = rownames(model_summary$coefficients),
      b = model_summary$coefficients[, 1],
      se = model_summary$coefficients[, 2],
      pval = model_summary$coefficients[, 4],
      nsnp = NA
    ) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "harmonized_data\\$exposure_beta", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\(", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\)", ""))
  }

  for (i in 1:nrow(mvmr_results)) {
    if (mvmr_results$exposure[i] == "Intercept") {
      mvmr_results$nsnp[i] <- nrow(harmonized_data$exposure_pval)
    } else {
      if (include_intercept == TRUE) {
        mvmr_results$nsnp[i] <- length(which(harmonized_data$exposure_pval[, i - 1] < pvalue_threshold))
      } else {
        mvmr_results$nsnp[i] <- length(which(harmonized_data$exposure_pval[, i] < pvalue_threshold))
      }
    }
  }

  return(mvmr_results)
}
