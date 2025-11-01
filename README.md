# orthomr

Orthogonal polynomials for non-linear Mendelian randomization.

## Installation

You can install the development version of orthomr from GitHub:
``` r
# install.packages("devtools")
devtools::install_github("yourusername/orthomr", build_vignettes = TRUE)
```

## Example
```r
library(orthomr)

# Create orthogonal polynomials
exposure_values <- rnorm(100)
poly_result <- orthopol(degree = 3, exposure_values = exposure_values)

# See the vignette for a complete workflow
browseVignettes("orthomr")
```

## License

GPL (>= 3)
