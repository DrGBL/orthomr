#' Choose Number of Polynomial Degrees
#'
#' Progressively takes smaller p-value cutoffs for instrument selection until
#' all polynomial degrees have F-statistics above a threshold.
#'
#' @param mvmr_exposure_data Data frame. MVMR exposure data in TwoSampleMR format.
#' @param mvmr_outcome_data Data frame. MVMR outcome data in TwoSampleMR format.
#' @param log10_pvalue_range Numeric vector of length 2. Range of log10 p-values to test
#'   (default: c(7, 300)).
#' @param f_statistic_threshold Numeric. F-statistic threshold (default: 10).
#'
#' @return A list containing:
#'   \item{harmonized_data}{Harmonized GWAS data}
#'   \item{instrument_strength}{Instrument strength statistics}
#'   Returns NULL if no valid configuration is found.
#'
#' @export
#' @examples
#' \dontrun{
#' result <- choose_n_degrees(
#'   mvmr_exposure_data,
#'   mvmr_outcome_data,
#'   f_statistic_threshold = 10
#' )
#' }
choose_n_degrees <- function(mvmr_exposure_data,
                             mvmr_outcome_data,
                             log10_pvalue_range = c(7, 300),
                             f_statistic_threshold = 10) {
  exposure_data_working <- mvmr_exposure_data

  max_degree <- exposure_data_working %>%
    dplyr::select(exposure) %>%
    dplyr::distinct() %>%
    dplyr::mutate(exposure = as.numeric(stringr::str_replace(exposure, "deg_", ""))) %>%
    dplyr::pull(exposure) %>%
    max(.)

  for (j in max_degree:1) {
    exposure_data_working <- exposure_data_working %>%
      dplyr::filter(id.exposure != paste0("deg_", as.character(j + 1)))

    for (i in seq(log10_pvalue_range[1], log10_pvalue_range[2], 1)) {
      # Suppress messages and warnings safely
      harmonized_gwas <- suppressWarnings(suppressMessages(
        TwoSampleMR::mv_harmonise_data(
          exposure_data_working %>%
            dplyr::group_by(SNP) %>%
            dplyr::mutate(min_pval = min(pval.exposure)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(min_pval < exp(-i)) %>%
            dplyr::group_by(id.exposure) %>%
            dplyr::mutate(min_pval = min(pval.exposure)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(min_pval < exp(-i)),
          mvmr_outcome_data
        )
      ))

      if (nrow(harmonized_gwas$expname) < j) {
        break()
      }

      # Strength instruments
      formatted_data <- MVMR::format_mvmr(
        BXGs = harmonized_gwas$exposure_beta,
        BYG = harmonized_gwas$outcome_beta,
        seBXGs = harmonized_gwas$exposure_se,
        seBYG = harmonized_gwas$outcome_se,
        RSID = paste("rs", c(1:nrow(harmonized_gwas$exposure_beta)))
      )
      instrument_strength <- MVMR::strength_mvmr(r_input = formatted_data, gencov = 0)

      if (length(which(instrument_strength > f_statistic_threshold)) == length(instrument_strength[1, ])) {
        return(list(harmonized_gwas, instrument_strength))
      }
    }
  }
  return()
}
