#' Coordinates from the Avecnet trial
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' }
"AvecNet_coordinates"

#' Test dataset containing the results of running DesignCRT() with the default options on a cut-down subset of the AvecNet dataset
#' @format list object with the following attributes:
#' \itemize{
#' \item \code{arm}: vector of assignments to trial arms
#' \item \code{alpha}: confidence level
#' \item \code{power}: power
#' \item \code{seed}: random number seed
#' \item \code{effect}: Required effect size
#' \item \code{ICC}: Intra-Cluster Correlation obtained from other studies
#' \item \code{DE}: calculated Design Effect
#' \item \code{pC}: baseline prevalence
#' \item \code{cont}: contamination range in km, obtained from other studies
#' \item \code{coordinate_source}: filename for coordinates of households
#' \item \code{h}: proposal for the number of coordinates in each cluster
#' \item \code{algo}: algorithm used for cluster boundaries
#' \item \code{n_ind}: algorithm used for cluster boundaries
#' \item \code{assignments}: data frame containing locations and cluster assignments
#' \item \code{min_c}: minimum number of clusters required
#' }
"test_AvecNet"

#' Dummy
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' }
"test_result"

#' Test dataset containing the results of running Simulate_TrialSite() with the default options
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' }
"test_site"

#' Results of toy GEE analysis for use in testing
#' @format list object with the following attributes:
#' \itemize{
#' \item \code{description}: Description of the trial dataset
#' \item \code{method}: statistical method
#' \item \code{PointEstimates}: point estimates
#' \item \code{IntervalEstimates}: interval estimates
#' \item \code{ModelObject}: object returned by the fitting function
#' }
"test_Analyse_CRT"

#' Test dataset containing the results of running Simulate_CRT() with the default options
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{cluster}: cluster assignment
#' \item \code{arm}: assignment to trial arm
#' \item \code{infectiousness_proxy}: infectiousness proxy
#' \item \code{nearestDiscord}: distance to nearest discordant location
#' \item \code{p0}: value of underlying prevalence surface without intervention
#' \item \code{p1}: value of underlying prevalence surface with intervention
#' \item \code{denom}: number of samples evaluated at the location
#' \item \code{num}: number of positive samples
#' \item \code{index}: row number
#' }
"test_Simulate_CRT"

#' Trial dataframe containing cluster assignments based on Rusinga locations
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{cluster}: cluster assignment
#' }
"testClusters"

#' Trial dataframe containing arm assignments based on Rusinga locations
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{cluster}: cluster assignment
#' \item \code{arm}: assignment to trial arm
#' }
"testArms"

#' Trial dataframe containing buffer assignment based on Rusinga locations
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{cluster}: cluster assignment
#' \item \code{arm}: assignment to trial arm
#' \item \code{nearestDiscord}: distance to nearest discordant location
#' \item \code{buffer}: logical indicator of whether a point is in the buffer
#' \item \code{index}: row number
#' }
"testBuffer"

#' Trial dataframe including simulated outcomes based on Rusinga locations
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{cluster}: cluster assignment
#' \item \code{arm}: assignment to trial arm
#' \item \code{infectiousness_proxy}: infectiousness proxy
#' \item \code{nearestDiscord}: distance to nearest discordant location
#' \item \code{p0}: value of underlying prevalence surface without intervention
#' \item \code{p1}: value of underlying prevalence surface with intervention
#' \item \code{buffer}: logical indicator of whether a point is in the buffer
#' \item \code{denom}: number of samples evaluated at the location
#' \item \code{num}: number of positive samples
#' \item \code{index}: row number
#' }
"testOutcomes"

#' Baseline RDT positivity data from the SolarMal trial
#' @format dataframe:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{RDT_test_result}: test result (0= negative, 1 = positive)
#' \item \code{X}: row number
#' }
"Solarmal_baseline"
