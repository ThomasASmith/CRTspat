#' Prediction mesh (created with ncells= 100) for INLA analysis of AvecNet geography
#' @format list:
#' \itemize{
#' \item \code{prediction}: data.frame containing the prediction points and covariate values
#' \item \code{A}: projection matrix from the observations to the mesh nodes.
#' \item \code{Ap}: projection matrix from the prediction points to the mesh nodes.
#' \item \code{indexs}:  index set for the SPDE model
#' \item \code{spde}: SPDE model
#' })
"inlaMesh100"

#' Prediction mesh (created with ncells= 100) for INLA analysis of test (Rusinga) geography
#' @format list:
#' \itemize{
#' \item \code{prediction}: data.frame containing the prediction points and covariate values
#' \item \code{A}: projection matrix from the observations to the mesh nodes.
#' \item \code{Ap}: projection matrix from the prediction points to the mesh nodes.
#' \item \code{indexs}:  index set for the SPDE model
#' \item \code{spde}: SPDE model
#' })
"inlaMeshTest"



