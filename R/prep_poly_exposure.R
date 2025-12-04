#' Prepare Polynomial Exposure for GWAS
#'
#' Reads phenotype data, creates orthogonal polynomials, and saves GWAS-ready
#' output files along with polynomial coefficients for later use in MVMR.
#'
#' @param phenotype_file Character. Path to phenotype file. Must be a tab-separated
#'   file with either 2 columns (IID, phenotype) or 3 columns (FID, IID, phenotype).
#' @param output_directory Character. Path to output directory where files will be saved.
#' @param degree Integer. Degree of polynomial (must be >= 2).
#' @param output_prefix Character. Optional output prefix for generated files (default: NULL).
#'   If NULL, files are named "polynomial_phenotype_df.tsv.gz" and "polynomial_coef_list.RDS".
#'   If provided, files are named "{prefix}_polynomial_phenotype_df.tsv.gz" and
#'   "{prefix}_polynomial_coef_list.RDS".
#' @param force_zero_intercept Logical. If TRUE, forces the constant term to zero in
#'   orthogonal polynomials (default: FALSE).
#'
#' @return Invisibly returns a list containing:
#'   \item{polynomial_data}{Data frame ready for GWAS with FID, IID, and polynomial columns}
#'   \item{coefficients}{List of polynomial coefficients for use in MVMR}
#'   \item{phenotype_file}{Path to saved phenotype file}
#'   \item{coefficient_file}{Path to saved coefficient file}
#'
#' @details
#' This function is designed to prepare exposure data for non-linear Mendelian
#' randomization. It:
#' \enumerate{
#'   \item Reads phenotype data (with 2 or 3 columns)
#'   \item Standardizes the phenotype values
#'   \item Creates orthogonal polynomials up to specified degree
#'   \item Saves a GWAS-ready file with polynomial values
#'   \item Saves polynomial coefficients for later use in MVMR analysis
#' }
#'
#' The output phenotype file can be used directly with PLINK or other GWAS software.
#' The coefficient file should be loaded after MVMR to obtain final effect estimates.
#'
#' @export
#' @examples
#' \dontrun{
#' # With 2-column phenotype file (IID, phenotype)
#' prep_poly_exposure(
#'   phenotype_file = "vitd.tsv.gz",
#'   output_directory = "/output/path/",
#'   degree = 5,
#'   output_prefix = "vitd"
#' )
#'
#' # With 3-column phenotype file (FID, IID, phenotype)
#' prep_poly_exposure(
#'   phenotype_file = "exposure.tsv",
#'   output_directory = "./results/",
#'   degree = 3
#' )
#' }
prep_poly_exposure <- function(phenotype_file,
                               output_directory,
                               degree,
                               output_prefix = NULL,
                               force_zero_intercept = FALSE) {

  # Input validation
  if (!file.exists(phenotype_file)) {
    stop("Phenotype file does not exist: ", phenotype_file, call. = FALSE)
  }

  if (!dir.exists(output_directory)) {
    stop("Output directory does not exist: ", output_directory, call. = FALSE)
  }

  if (degree < 2) {
    stop("The value of degree must be an integer greater than or equal to 2.", call. = FALSE)
  }

  # Load phenotype data
  message("Reading phenotype file: ", phenotype_file)
  phenotype_data <- readr::read_tsv(phenotype_file, show_col_types = FALSE)

  # Validate number of columns
  if (ncol(phenotype_data) == 1 || ncol(phenotype_data) > 3) {
    stop("Incorrect number of columns in the phenotype file. Must have 2 or 3 columns.",
         call. = FALSE)
  }

  # Create orthogonal polynomials
  message("Computing orthogonal polynomials of degree ", degree)

  if (ncol(phenotype_data) == 2) {
    # 2 columns: IID, phenotype
    exposure_scaled <- scale(dplyr::pull(phenotype_data[, 2]))
    poly_result <- orthopol(
      degree = degree,
      exposure_values = exposure_scaled,
      force_zero_intercept = force_zero_intercept
    )

    polynomial_df <- data.frame(
      FID = dplyr::pull(phenotype_data[, 1]),
      IID = dplyr::pull(phenotype_data[, 1])
    ) %>%
      dplyr::bind_cols(., as.data.frame(poly_result[["values"]]))

  } else {
    # 3 columns: FID, IID, phenotype
    exposure_scaled <- scale(dplyr::pull(phenotype_data[, 3]))
    poly_result <- orthopol(
      degree = degree,
      exposure_values = exposure_scaled,
      force_zero_intercept = force_zero_intercept
    )

    polynomial_df <- data.frame(
      FID = dplyr::pull(phenotype_data[, 1]),
      IID = dplyr::pull(phenotype_data[, 2])
    ) %>%
      dplyr::bind_cols(., as.data.frame(poly_result[["values"]]))
  }

  # Construct output file paths
  if (is.null(output_prefix)) {
    phenotype_output <- file.path(output_directory, "polynomial_phenotype_df.tsv.gz")
    coefficient_output <- file.path(output_directory, "polynomial_coef_list.RDS")
  } else {
    phenotype_output <- file.path(
      output_directory,
      paste0(output_prefix, "_polynomial_phenotype_df.tsv.gz")
    )
    coefficient_output <- file.path(
      output_directory,
      paste0(output_prefix, "_polynomial_coef_list.RDS")
    )
  }

  # Save outputs
  message("Saving outputs...")
  readr::write_tsv(polynomial_df, phenotype_output)
  saveRDS(poly_result[["coefficients"]], coefficient_output)

  message("Polynomial phenotype GWAS-ready dataframe saved to: ", phenotype_output)
  message("Polynomial coefficients (to be used after MVMR step) saved to: ", coefficient_output)

  # Return invisibly
  invisible(list(
    polynomial_data = polynomial_df,
    coefficients = poly_result[["coefficients"]],
    phenotype_file = phenotype_output,
    coefficient_file = coefficient_output
  ))
}
