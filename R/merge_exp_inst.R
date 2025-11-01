#' Merge Exposure Instruments
#'
#' Loops through GWAS exposure list, extracts significant instruments from all
#' GWAS, merges them, and removes duplicates.
#'
#' @param gwas_exposure_list List of data frames. Each data frame should contain columns:
#'   instrument_id, beta, se, pval, allele_freq.
#'
#' @return A data frame with merged instruments in TwoSampleMR format with columns:
#'   SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure,
#'   beta.exposure, se.exposure, pval.exposure, eaf.exposure.
#'
#' @export
#' @examples
#' \dontrun{
#' gwas_list <- list(gwas1, gwas2, gwas3)
#' merged_instruments <- merge_exp_inst(gwas_list)
#' }
merge_exp_inst <- function(gwas_exposure_list) {
  merged_instruments <- c()

  for (i in 1:length(gwas_exposure_list)) {
    significant_instruments <- gwas_exposure_list[[i]] %>%
      dplyr::filter(pval < 5e-8) %>%
      dplyr::pull(ins)

    merged_instruments <- gwas_exposure_list[[i]] %>%
      dplyr::filter(ins %in% significant_instruments) %>%
      dplyr::rename(SNP = ins) %>%
      dplyr::mutate(exposure = paste0("deg_", i)) %>%
      dplyr::mutate(id.exposure = paste0("deg_", i)) %>%
      dplyr::mutate(effect_allele.exposure = "A") %>%
      dplyr::mutate(other_allele.exposure = "C") %>%
      dplyr::rename(beta.exposure = beta) %>%
      dplyr::rename(se.exposure = se) %>%
      dplyr::rename(pval.exposure = pval) %>%
      dplyr::rename(eaf.exposure = af) %>%
      dplyr::bind_rows(merged_instruments, .)

    for (j in 1:length(gwas_exposure_list)) {
      if (i != j) {
        merged_instruments <- gwas_exposure_list[[j]] %>%
          dplyr::filter(ins %in% significant_instruments) %>%
          dplyr::rename(SNP = ins) %>%
          dplyr::mutate(exposure = paste0("deg_", j)) %>%
          dplyr::mutate(id.exposure = paste0("deg_", j)) %>%
          dplyr::mutate(effect_allele.exposure = "A") %>%
          dplyr::mutate(other_allele.exposure = "C") %>%
          dplyr::rename(beta.exposure = beta) %>%
          dplyr::rename(se.exposure = se) %>%
          dplyr::rename(pval.exposure = pval) %>%
          dplyr::rename(eaf.exposure = af) %>%
          dplyr::bind_rows(merged_instruments, .)
      }
    }
  }
  merged_instruments <- merged_instruments %>% dplyr::distinct()
  return(merged_instruments)
}
