#' Obtain Final Coefficients
#'
#' Combines orthogonal polynomial coefficients with MVMR coefficients to obtain
#' final polynomial coefficients.
#'
#' @param polyX List. Output from \code{orthopol()$coef}.
#' @param mvmr_res Data frame. MVMR results with columns: exposure, b, pval, nsnp.
#' @param pval Numeric. P-value threshold for including polynomial terms (default: 0.05).
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
#' final_coef <- obtain_final_coeffs(polyX$coef, mvmr_res, pval = 0.05)
#' }
obtain_final_coeffs <- function(polyX, mvmr_res, pval,
                                set_higher_to_zero = FALSE,
                                remove_degree_if_no_snps = FALSE) {
  final_coefs <- polyX

  if (mvmr_res$exposure[1] == "Intercept") {
    intercept_df <- mvmr_res[1, ]
    mvmr_res_tmp <- mvmr_res[-1, ]
  } else {
    mvmr_res_tmp <- mvmr_res
  }

  if (remove_degree_if_no_snps) {
    mvmr_res_tmp <- mvmr_res_tmp %>%
      dplyr::mutate(degree_int = as.numeric(stringr::str_extract(exposure, "[0-9]*$"))) %>%
      dplyr::arrange(degree_int) %>%
      dplyr::mutate(keep_deg = FALSE)
    for (i in 1:nrow(mvmr_res_tmp)) {
      if (i == 1) {
        if (mvmr_res_tmp$nsnp[i] > 0) {
          mvmr_res_tmp$keep_deg[i] <- TRUE
        }
      } else {
        if (mvmr_res_tmp$nsnp[i] > 0 & mvmr_res_tmp$keep_deg[i - 1] == TRUE) {
          mvmr_res_tmp$keep_deg[i] <- TRUE
        }
      }
    }
    mvmr_res_tmp <- mvmr_res_tmp %>% dplyr::filter(keep_deg == TRUE)
    final_coefs <- final_coefs[1:nrow(mvmr_res_tmp)]
  }

  if (set_higher_to_zero == TRUE) {
    stop_coeff <- FALSE
    for (i in 1:nrow(mvmr_res_tmp)) {
      if (mvmr_res_tmp$pval[which(mvmr_res_tmp$exposure == paste0("deg_", i))] < pval &
          stop_coeff == FALSE) {
        final_coefs[[i]] <- final_coefs[[i]] *
          mvmr_res_tmp$b[which(mvmr_res_tmp$exposure == paste0("deg_", i))]
      } else {
        final_coefs[[i]] <- final_coefs[[i]] * 0
        stop_coeff <- TRUE
      }
    }
  } else {
    for (i in 1:nrow(mvmr_res_tmp)) {
      if (mvmr_res_tmp$pval[which(mvmr_res_tmp$exposure == paste0("deg_", i))] < pval) {
        final_coefs[[i]] <- final_coefs[[i]] *
          mvmr_res_tmp$b[which(mvmr_res_tmp$exposure == paste0("deg_", i))]
      } else {
        final_coefs[[i]] <- final_coefs[[i]] * 0
      }
    }
  }

  if (mvmr_res$exposure[1] == "Intercept") {
    final_coefs[["Intercept"]] <- intercept_df$b
  }

  return(final_coefs)
}
