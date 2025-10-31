#' Choose Number of Polynomial Degrees
#'
#' Progressively takes smaller p-value cutoffs for instrument selection until
#' all polynomial degrees have F-statistics above a threshold.
#'
#' @param mvmr_exp Data frame. MVMR exposure data in TwoSampleMR format.
#' @param mvmr_out Data frame. MVMR outcome data in TwoSampleMR format.
#' @param log10p_range Numeric vector of length 2. Range of log10 p-values to test
#'   (default: c(7, 300)).
#' @param F_thresh Numeric. F-statistic threshold (default: 10).
#'
#' @return A list containing:
#'   \item{gwas_harm}{Harmonized GWAS data}
#'   \item{strength}{Instrument strength statistics}
#'   Returns NULL if no valid configuration is found.
#'
#' @export
#' @examples
#' \dontrun{
#' result <- choose_n_degrees(mvmr_exp, mvmr_out, F_thresh = 10)
#' }
choose_n_degrees <- function(mvmr_exp,
                             mvmr_out,
                             log10p_range = c(7, 300),
                             F_thresh = 10) {
  mvmr_exp_tmp <- mvmr_exp

  max_deg <- mvmr_exp_tmp %>%
    dplyr::select(exposure) %>%
    dplyr::distinct() %>%
    dplyr::mutate(exposure = as.numeric(stringr::str_replace(exposure, "deg_", ""))) %>%
    dplyr::pull(exposure) %>%
    max(.)

  for (j in max_deg:1) {
    mvmr_exp_tmp <- mvmr_exp_tmp %>%
      dplyr::filter(id.exposure != paste0("deg_", as.character(j + 1)))

    for (i in seq(log10p_range[1], log10p_range[2], 1)) {
      sink(nullfile())
      suppressMessages(
        gwas_harm <- TwoSampleMR::mv_harmonise_data(
          mvmr_exp_tmp %>%
            dplyr::group_by(SNP) %>%
            dplyr::mutate(min_pval = min(pval.exposure)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(min_pval < exp(-i)) %>%
            dplyr::group_by(id.exposure) %>%
            dplyr::mutate(min_pval = min(pval.exposure)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(min_pval < exp(-i)),
          mvmr_out
        )
      )

      if (nrow(gwas_harm$expname) < j) {
        break()
      }

      # Strength instruments
      Fdata <- MVMR::format_mvmr(
        BXGs = gwas_harm$exposure_beta,
        BYG = gwas_harm$outcome_beta,
        seBXGs = gwas_harm$exposure_se,
        seBYG = gwas_harm$outcome_se,
        RSID = paste("rs", c(1:nrow(gwas_harm$exposure_beta)))
      )
      stre <- MVMR::strength_mvmr(r_input = Fdata, gencov = 0)
      sink()

      if (length(which(stre > 10)) == length(stre[1, ])) {
        return(list(gwas_harm, stre))
      }
    }
  }
  return()
}
