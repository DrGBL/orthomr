#' Run Multivariable Mendelian Randomization
#'
#' Performs MVMR analysis using harmonized GWAS data.
#'
#' @param gwas_harm List. Harmonized GWAS data from \code{TwoSampleMR::mv_harmonise_data()}
#'   or \code{choose_n_degrees()}.
#' @param intercept Logical. Whether to include an intercept term (default: TRUE).
#' @param pvalue Numeric. P-value threshold for counting significant SNPs (default: 5e-8).
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
#' mvmr_results <- run_mvmr(gwas_harm, intercept = TRUE)
#' }
run_mvmr <- function(gwas_harm, intercept = TRUE, pvalue = 5e-8) {

  if (intercept == TRUE) {
    summary <- summary(lm(gwas_harm$outcome_beta ~ gwas_harm$exposure_beta,
                          weights = 1 / gwas_harm$outcome_se^(2)))
    mvmr_res <- data.frame(
      exposure = rownames(summary$coefficients),
      b = summary$coefficients[, 1],
      se = summary$coefficients[, 2],
      pval = summary$coefficients[, 4],
      nsnp = NA
    ) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "gwas_harm\\$exposure_beta", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\(", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\)", ""))
  } else {
    summary <- summary(lm(gwas_harm$outcome_beta ~ gwas_harm$exposure_beta - 1,
                          weights = 1 / gwas_harm$outcome_se^(2)))
    mvmr_res <- data.frame(
      exposure = rownames(summary$coefficients),
      b = summary$coefficients[, 1],
      se = summary$coefficients[, 2],
      pval = summary$coefficients[, 4],
      nsnp = NA
    ) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "gwas_harm\\$exposure_beta", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\(", "")) %>%
      dplyr::mutate(exposure = stringr::str_replace(exposure, "\\)", ""))
  }

  for (i in 1:nrow(mvmr_res)) {
    if (mvmr_res$exposure[i] == "Intercept") {
      mvmr_res$nsnp[i] <- nrow(gwas_harm$exposure_pval)
    } else {
      if (intercept == TRUE) {
        mvmr_res$nsnp[i] <- length(which(gwas_harm$exposure_pval[, i - 1] < pvalue))
      } else {
        mvmr_res$nsnp[i] <- length(which(gwas_harm$exposure_pval[, i] < pvalue))
      }
    }
  }

  return(mvmr_res)
}
