#' The example data is meant to represent one dataset for Scenario II in simulation study,
#' which is explored in Section S3 of Supplementary Materials in the paper.
#' The 'sampledata' file contains 378 rows and 83 variables.
#'
#' @format A data frame with 378 rows and 83 variables:
#' \describe{
#'   \item{First  column}{Yvariable Binary indicating dependent variable for the couple disease status, 1 for disease}
#'   \item{2:37 columns}{X.f_mat chemical exposure variables for female individual}
#'   \item{38:73 columns}{X.m_mat chemical exposure variables for male individual}
#'   \item{74:78 columns}{covariate_f_mat subject-specific covariates such as age or smoking status for female individual}
#'   \item{79:83 columns}{covariate_m_mat subject-specific covariates such as age or smoking status for male individual}
#'    }
#' @source \url{https://doi.org/10.1111/biom.12972}
"sampledata"
