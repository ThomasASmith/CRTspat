#' Design a CRT of a malaria intervention with contamination
#'
#' \code{designCRT} estimates the required number of clusters and the extent of contamination between arms for a CRT
#' with the input set of locations. Outputs are:
#' (i) Estimates of the required numbers of clusters.
#' (ii) A proposal for the cluster and arm assignments to the input coordinates.
#' (A warning is output if the number of locations is too small to allow randomisation of sufficient clusters).
#' (iii) the proportion of households in the input geography falling within the core of the clusters
#' (i.e. outside the contamination range of locations in the opposite arm)
#' @param trial dataframe or CRT object containing Cartesian coordinates of locations in columns 'x' and 'y'.
#' Units are expected to be km.
#' @param alpha confidence level
#' @param desiredPower desired power
#' @param effect required effect size
#' @param ICC Intra-Cluster Correlation
#' @param pC baseline prevalence
#' @param outcome.type options are 'Continuous', 'Count or rate', 'Proportion', 'Dichotomous'
#' @param buffer.width required buffer width in km
#' @param h  proposal for the number of coordinates in each cluster
#' @param algorithm algorithm used to determine cluster boundaries, options are-
#' 'TSP': travelling salesman problem heuristic;
#' 'NN': nearest neighbor;
#' 'kmeans': kmeans
#' @return A CRT object comprising the input data, cluster and arm assignments, trial description and results of power calculations
#' \itemize{
#' \item \code{arm}: arm assignment
#' \item \code{cluster}: cluster assignment
#' \item \code{buffer}: indicator of whether the location is in the buffer between arms
#' \item \code{input.parameters} list of input parameters
#' \item \code{CRT.design.full} list of calculated filename for coordinates of households
#' }
#' @export
#' @examples
#'
#' exampleDesign = designCRT(trial=readdata('test_site.csv'),
#'                 ICC=0.10, effect=0.4, pC=0.35, desiredPower = 0.8, algorithm = 'kmeans',
#'                 buffer.width=0.25, h=100, outcome.type ='Dichotomous')
designCRT <- function(trial, alpha = 0.05, desiredPower = 0.8, effect, ICC,
    pC, buffer.width = 0, h, algorithm = "kmeans", outcome.type = outcome.type) {

    input.parameters <- list(pC = pC, alpha = alpha, desiredPower = desiredPower,
        h = h, algorithm = algorithm, outcome.type = outcome.type, buffer.width = buffer.width,
        effect = effect, ICC = ICC)

    # Step J: specify or compute cluster boundaries
    CRT <- specify.clusters(trial = trial, h = h, algorithm = algorithm, reuseTSP = FALSE)

    # Step K: Random assignment of clusters to arms
    CRT <- randomizeCRT(CRT)

    # augment the trial data frame with distance to nearest discordant
    # coordinate (specifyBuffer assigns a buffer only if a buffer.width> 0 is input)

    CRT <- specify.buffer(trial = CRT, buffer.width = buffer.width)
    CRT <- convert.data.frame.CRT(CRT, input.parameters = input.parameters)
    return(CRT)
}


#' Power and sample size calculations for a CRT
#'
#' \code{calculateCRTpower} carries out power and sample-size calculations for a CRT
#' The output is a CRT object containing values for:
#' - The required numbers of clusters to achieve the specified power.
#' - The design effect based on the input ICC.
#' - The power achievable.
#' @param locations total number of units
#' @param ICC Intracluster Correlation Coefficient
#' @param effect Required effect size
#' @param alpha Significance level
#' @param outcome.type options are 'Continuous', 'Count or rate', 'Proportion', 'Dichotomous'
#' @param desiredPower desired power (expressed as a proportion)
#' @param pC expected average outcome in control arm
#' @param mean_h mean number of units per cluster
#' @param sd_h standard deviation of number of units per cluster
#' @return A CRT object comprising the input coordinates :
#' \itemize{
#' \item \code{input.parameters} list of input parameters
#' \item \code{CRT.design.full} list of calculated filename for coordinates of households
#' }
#' @export
#' @examples
#'
#' examplePower = calculateCRTpower(locations = 3000, ICC=0.10, effect=0.4, alpha = 0.05,
#'     outcome.type = 'Dichotomous', desiredPower = 0.8, pC=0.35, mean_h=100, sd_h=5)
calculateCRTpower <- function(locations, ICC, effect, alpha, outcome.type, desiredPower,
    pC, mean_h, sd_h) {
    # TODO: add non-binomial outcomes

    # Step A: confidence level Step B: power Step C: Required effect
    # size Step D: ICC, obtained from other studies Step E: baseline
    # prevalence Step F: buffer width based on postulated contamination
    # range in km, obtained from other studies Step G\ coordinates of
    # households in study area Step H: proposal for the number of
    # households in each cluster

    # convert power and significance level to normal deviates
    Zsig <- -qnorm(alpha/2)
    Zpow <- qnorm(desiredPower)

    # Power and sample size calculations based on Hemming et al, 2011
    # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102

    pI <- pC * (1 - effect)  # probability in intervened group
    d <- pC - pI  # difference between groups
    sigma2 <- 1/2 * (pI * (1 - pI) + pC * (1 - pC))

    # required individuals per arm in individually randomized trial
    n_ind <- 2 * sigma2 * ((Zsig + Zpow)/d)^2

    # number of clusters assigned to each arm
    k <- floor(locations/(mean_h * 2))

    # coefficient of variation of the cluster sizes
    cv_h <- sd_h/mean_h

    # design effect (Variance Inflation Factor) allowing for varying
    # cluster sizes(Hemming eqn 6)
    DE <- 1 + (cv_h^2 + 1) * (mean_h - 1) * ICC

    # number of individuals required per arm in CRT with equal cluster
    # sizes
    n_crt <- n_ind * DE

    # Step I: Calculations for the required minimum number of clusters
    # for both arms

    # minimum numbers of clusters required assuming varying cluster
    # sizes per arm (Hemming eqn 8)
    min_k <- ceiling(n_ind * (1 + ((cv_h + 1) * mean_h - 1) * ICC)/mean_h)

    clustersRequired <- 2 * min_k

    # power with k clusters per arm power =
    # stats::pnorm(sqrt(k*mean_h/(2*(1+(mean_h-1)*ICC)))*
    # d/sqrt(sigma2) - Zsig) #Hemming eqn 27

    power <- stats::pnorm(sqrt(k * mean_h/(2 * DE)) * d/sqrt(sigma2) - Zsig)  #unequal cluster sizes

    CRT.design.full <- list(locations = locations, mean_h = mean_h, sd_h = sd_h,
        min_k = min_k, clustersRequired = clustersRequired, DE = DE, power = power)
    input.parameters <- list(locations = locations, ICC = ICC, effect = effect,
        alpha = alpha, outcome.type = outcome.type, desiredPower = desiredPower,
        pC = pC, mean_h = mean_h, sd_h = sd_h)
    CRT <- list(CRT.design.full = CRT.design.full, input.parameters = input.parameters)
    class(CRT) <- "CRT"
    return(CRT)
}
