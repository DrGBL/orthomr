#' Merge Exposure Instruments
#'
#' Loops through GWAS exposure list, extracts significant instruments from all
#' GWAS, merges them, and removes duplicates.
#'
#' @param gwas_exp List of data frames. Each data frame should contain columns:
#'   ins, beta, se, pval, af.
#'
#' @return A data frame with merged instruments in TwoSampleMR format with columns:
#'   SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure,
#'   beta.exposure, se.exposure, pval.exposure, eaf.exposure.
#'
#' @export
#' @examples
#' \dontrun{
#' gwas_list <- list(gwas1, gwas2, gwas3)
#' merged_inst <- merge_exp_inst(gwas_list)
#' }
merge_exp_inst <- function(gwas_exp) {
  inst_pre_mvmr <- c()
  for (i in 1:length(gwas_exp)) {
    sig_ins <- gwas_exp[[i]] %>%
      dplyr::filter(pval < 5e-8) %>%
      dplyr::pull(ins)

    inst_pre_mvmr <- gwas_exp[[i]] %>%
      dplyr::filter(ins %in% sig_ins) %>%
      dplyr::rename(SNP = ins) %>%
      dplyr::mutate(exposure = paste0("deg_", i)) %>%
      dplyr::mutate(id.exposure = paste0("deg_", i)) %>%
      dplyr::mutate(effect_allele.exposure = "A") %>%
      dplyr::mutate(other_allele.exposure = "C") %>%
      dplyr::rename(beta.exposure = beta) %>%
      dplyr::rename(se.exposure = se) %>%
      dplyr::rename(pval.exposure = pval) %>%
      dplyr::rename(eaf.exposure = af) %>%
      dplyr::bind_rows(inst_pre_mvmr, .)

    for (j in 1:length(gwas_exp)) {
      if (i != j) {
        inst_pre_mvmr <- gwas_exp[[j]] %>%
          dplyr::filter(ins %in% sig_ins) %>%
          dplyr::rename(SNP = ins) %>%
          dplyr::mutate(exposure = paste0("deg_", j)) %>%
          dplyr::mutate(id.exposure = paste0("deg_", j)) %>%
          dplyr::mutate(effect_allele.exposure = "A") %>%
          dplyr::mutate(other_allele.exposure = "C") %>%
          dplyr::rename(beta.exposure = beta) %>%
          dplyr::rename(se.exposure = se) %>%
          dplyr::rename(pval.exposure = pval) %>%
          dplyr::rename(eaf.exposure = af) %>%
          dplyr::bind_rows(inst_pre_mvmr, .)
      }
    }
  }
  inst_pre_mvmr <- inst_pre_mvmr %>% dplyr::distinct()
  return(inst_pre_mvmr)
}
