#' Designs, simulates, and analyses a cluster randomised trial (CRT) with geographical contamination
#'
#' \code{CRTspillover} estimates the required number of clusters and the extent of contamination between arms for a CRT based on an input set of locations.
#' The trial is simulated and the simulated data are analysed.
#'
#' @param alpha confidence level
#' @param power power
#' @param effect required effect size
#' @param ICC Intra-Cluster Correlation obtained from other studies
#' @param pC baseline prevalence
#' @param cont contamination range in km, obtained from other studies
#' @param coordinates dataframe containing coordinates of households. Columns 'x' and 'y' should contain Cartesian (x,y) coordinates. Units are expected to be km.
#' @param h  proposal for the number of coordinates in each cluster
#' @param algo algorithm for cluster boundaries, choose between
#' "TSP": travelling salesman problem heuristic;
#' "NN": nearest neighbor;
#' "kmeans": kmeans algorithm
#' @param reuseTSP indicator of whether a pre-existing path should be used by the TSP algorithm
#' @param sd standard deviation measuring spatial smoothing of proportion
#' @param method statistical method used to analyse trial.
#' Options are 'piecewise_linear','logit','sigmoid','empirical','GEE','MCMC01','MCMC02','MCMC03'
#' @param excludeBuffer exclude any buffer zone (records with buffer=TRUE) from the analysis
#' @param requireBootstrap logical indicator of whether bootstrap confidence intervals are required
#' @param resamples number of bootstrap samples if bootstrap confidence intervals are required
#' @param burnin number of burnin interations of MCMC algorithm (for MCMC methods)
#' @param iter total number of interations of MCMC algorithm (for MCMC methods)
#' @return list object with the following attributes:
#' \itemize{
#' \item \code{design}: list of vector-specific entomological parameters (see
#' description of output from function \code{Design_CRT})
#' \item \code{simulation}: vector of intervention-specific entomological
#' parameters and impact (vectorial capacity) values
#' \item \code{analysis}: list of results of statistical analysis of the simulated trial (see
#' description of output from function \code{Analyse_CRT})
#' }
#' @export
#'
#' @examples
#'
CRTspillover = function(alpha = 0.05, power = 0.8, effect = 0.4, ICC = 0.175,
  pC = 0.4, cont = 0.25, coordinates=AvecNet_coordinates, h = 80, algo = "kmeans",
  reuseTSP=FALSE, sd = 0.4, method='MCMC03', excludeBuffer=FALSE, requireBootstrap=FALSE,
  resamples=1000, iter=10000, burnin=5000){

  design =Design_CRT(alpha = alpha, power = power, effect = effect, ICC = ICC,
              pC = pC, cont = cont, coordinates=coordinates, h = h, algo = algo, reuseTSP=reuseTSP)

  simulation = Simulate_CRT(trial=design$assignments, efficacy=effect, initialPrevalence=pC, sd=sd)

  analysis = Analyse_CRT(trial=simulation, method=method, excludeBuffer=excludeBuffer, requireBootstrap=requireBootstrap,
              alpha=alpha, resamples=resamples, burnin=burnin)

  CRT=list(design=design,simulation=simulation,analysis=analysis)
return(CRT)}


