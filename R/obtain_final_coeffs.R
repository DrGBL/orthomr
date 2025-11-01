#' Obtain Final Polynomial Coefficients
#'
#' Combines orthogonal polynomial coefficients with MVMR coefficients to obtain
#' final polynomial coefficients.
#'
#' @param orthogonal_poly_coef List. Output from \code{orthopol()$coefficients}.
#' @param mvmr_results Data frame. MVMR results with columns: exposure, b, pval, nsnp.
#' @param pvalue_threshold Numeric. P-value threshold for including polynomial terms (default: 0.05).
#' @param set_higher_to_zero Logical. If TRUE, sets all higher-order terms to zero
#'   once a non-significant term is encountered (default: FALSE).
#' @param remove_degree_if_no_snps Logical. If TRUE, removes polynomial degrees
#'   with no significant SNPs (default: FALSE).
#'
#' @return A list of final coefficients for each polynomial degree.
#'
#' @export
#' @examples
#' \dontrun{
#' final_coef <- obtain_final_coeffs(
#'   poly_result$coefficients,
#'   mvmr_results,
#'   pvalue_threshold = 0.05
#' )
#' }
obtain_final_coeffs <- function(orthogonal_poly_coef,
                                mvmr_results,
                                pvalue_threshold,
                                set_higher_to_zero = FALSE,
                                remove_degree_if_no_snps = FALSE) {
  final_coefficients <- orthogonal_poly_coef

  if (mvmr_results$exposure[1] == "Intercept") {
    intercept_data <- mvmr_results[1, ]
    mvmr_results_filtered <- mvmr_results[-1, ]
  } else {
    mvmr_results_filtered <- mvmr_results
  }

  if (remove_degree_if_no_snps) {
    mvmr_results_filtered <- mvmr_results_filtered %>%
      dplyr::mutate(degree_int = as.numeric(stringr::str_extract(exposure, "[0-9]*$"))) %>%
      dplyr::arrange(degree_int) %>%
      dplyr::mutate(keep_deg = FALSE)
    for (i in 1:nrow(mvmr_results_filtered)) {
      if (i == 1) {
        if (mvmr_results_filtered$nsnp[i] > 0) {
          mvmr_results_filtered$keep_deg[i] <- TRUE
        }
      } else {
        if (mvmr_results_filtered$nsnp[i] > 0 & mvmr_results_filtered$keep_deg[i - 1] == TRUE) {
          mvmr_results_filtered$keep_deg[i] <- TRUE
        }
      }
    }
    mvmr_results_filtered <- mvmr_results_filtered %>% dplyr::filter(keep_deg == TRUE)
    final_coefficients <- final_coefficients[1:nrow(mvmr_results_filtered)]
  }

  if (set_higher_to_zero == TRUE) {
    stop_coefficient <- FALSE
    for (i in 1:nrow(mvmr_results_filtered)) {
      exposure_name <- paste0("deg_", i)
      exposure_index <- which(mvmr_results_filtered$exposure == exposure_name)

      if (mvmr_results_filtered$pval[exposure_index] < pvalue_threshold &
          stop_coefficient == FALSE) {
        final_coefficients[[i]] <- final_coefficients[[i]] *
          mvmr_results_filtered$b[exposure_index]
      } else {
        final_coefficients[[i]] <- final_coefficients[[i]] * 0
        stop_coefficient <- TRUE
      }
    }
  } else {
    for (i in 1:nrow(mvmr_results_filtered)) {
      exposure_name <- paste0("deg_", i)
      exposure_index <- which(mvmr_results_filtered$exposure == exposure_name)

      if (mvmr_results_filtered$pval[exposure_index] < pvalue_threshold) {
        final_coefficients[[i]] <- final_coefficients[[i]] *
          mvmr_results_filtered$b[exposure_index]
      } else {
        final_coefficients[[i]] <- final_coefficients[[i]] * 0
      }
    }
  }

  if (mvmr_results$exposure[1] == "Intercept") {
    final_coefficients[["Intercept"]] <- intercept_data$b
  }

  return(final_coefficients)
}
