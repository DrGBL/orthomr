# orthomr

Orthogonal polynomials for non-linear Mendelian randomization.

## Installation

You can install the development version of orthomr from GitHub:
``` r
# install.packages("remotes")
remotes::install_github("yourusername/orthomr", build_vignettes = TRUE)
```

## Example

The example below and the simulates data and GWAS, then uses orthomr.

```r
library(orthomr)
# See the vignette for a complete workflow
browseVignettes("orthomr")
```

To use the package on real data, you need to first create orthogonal polynomials by inputting a 3 columns file (```you_data_df.tsv.gz```) to the function ```prep_poly_exposure```:

Here, we used polynomials of degrees 1 through 5.

```r
prep_poly_exposure(
 phenotype_file = "you_data_df.tsv.gz",
 output_directory = "/output/path/",
 degree = 5,
 output_prefix = "output_prefix"
)
```

This will create a file called ```output_prefix_polynomial_phenotype_df.tsv.gz``` to perform your genome-wide association studies on (with your favourite software).

It also outputs an R object list of the polynomial coefficients and saves it in the same directory as ```output_prefix_polynomial_coef_list.RDS```.

You then need to LD-clump the GWAS so that all instruments for all phenotypes are in linkage equilibrium. A suggested helper file is provided here (```inst_select.sh```.), but this can be done differently by the user. The required input for this function are:

- ```--summ_stats_prefix```: prefix of orthogonal polynomial GWAS summary statistics files (see next line for details).

- ```--path_summ_stats```: path to the directory containing the GWAS summary statistics results. The summary statistics are assumed to all be in the same directory and be called ```path_summ_stats/summ_stats_prefix_chr{chrN}_deg_{degN}.tsv.gz```.

- ```--col_chr COLNAME```: Column name for chromosome in the summary statistics files.

- ```--col_snp COLNAME```: Column name for SNP in the summary statistics files.

- ```--col_bp COLNAME```: Column name for base pair position in the summary statistics files.

- ```--col_a1 COLNAME```: Column name for effect allele in the summary statistics files.

- ```--col_p COLNAME```: Column name for p-value in the summary statistics files.

- ```--path_tmp```: path to a directory where temporary files will be output. Defaults to your current directory otherwise.

- ```--path_ref```: path to the directory where the LD reference cohort files are located. These are assumed to be plink binary files (bim/bed/fam), one triplet per chromosome.

- ```--ref_prefix```: the prefix of the reference cohort plink files. These are assumed to be of the form ```ref_prefix_chr{chrN}```.

- ```--path_clumps```: path to output directory.

- ```--max_degree```: the maximum polynomial degree to work with.

Additional optional arguments are explained in the ```--help``` command.

The summary statistics of the chosen instruments can then be loaded using the TwoSampleMR package multivariable MR framework, using the ```mr_sim_res```  command as follows:

```r
poly_result<-readRDS("output_prefix_polynomial_coef_list.RDS")

real_results <- mr_sim_res(
  mvmr_exposure_data = mvmr_exposure_data,
  mvmr_outcome_data = mvmr_outcome_data,
  exposure_values = observed_exposure,
  orthogonal_poly_coef = poly_result$coefficients
)
```

## License

GPL (>= 3)
