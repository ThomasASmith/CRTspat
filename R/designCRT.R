#' Design a CRT for a given set of locations
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
#' @param yC baseline outcome
#' @param outcome.type options are:
#'  \tabular{ll}{
#' \code{"y"}\tab continous \cr
#' \code{"n"}\tab count \cr
#' \code{"e"}\tab event rate \cr
#' \code{"p"}\tab proportion  \cr
#' \code{"d"}\tab dichotomous \cr
#' }
#' @param buffer.width required buffer width in km
#' @param h  proposal for the number of coordinates in each cluster
#' @param algorithm algorithm used to determine cluster boundaries, options are-
#' \code{"TSP"}: travelling salesman problem heuristic;
#' \code{"NN"}: nearest neighbor;
#' \code{"kmeans"}: kmeans
#' @return A \code{"CRT"} object comprising the input data, cluster and arm assignments, trial description and results of power calculations
#'  \tabular{lll}{
#'  \code{CRT.design.full}   \tab list \tab summary statistics describing the site,
#'  cluster assignments, and randomization.\cr
#'  \code{x} \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y} \tab numeric vector \tab y-coordinates of locations \cr
#'  \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \code{arm} \tab factor \tab assignments to \code{"control"} or \code{"intervention"} for each location \cr
#'  \code{buffer} \tab logical \tab indicator of whether the point is within the buffer \cr
#'  \code{...} \tab other objects included in the input \code{"CRT"} object or data frame \cr
#'  }
#' @export
#' @examples
#' exampleDesign = designCRT(trial=readdata('test_site.csv'),
#'                 ICC=0.10, effect=0.4, yC=0.35, desiredPower = 0.8, algorithm = 'kmeans',
#'                 buffer.width=0.25, h=100, outcome.type ='d')
designCRT <- function(trial, alpha = 0.05, desiredPower = 0.8, effect, ICC,
    yC, buffer.width = 0, h, algorithm = "kmeans", outcome.type = outcome.type) {

    trial <- CRT_as_data.frame(trial)
    input.parameters <- list(yC = yC, alpha = alpha, desiredPower = desiredPower,
        h = h, algorithm = algorithm, outcome.type = outcome.type, buffer.width = buffer.width,
        effect = effect, ICC = ICC)

    # Step J: specify or compute cluster boundaries
    CRT <- specify.clusters(trial = trial, h = h, algorithm = algorithm, reuseTSP = FALSE)

    # Step K: Random assignment of clusters to arms
    CRT <- randomizeCRT(CRT)

    # augment the trial data frame with distance to nearest discordant
    # coordinate (specifyBuffer assigns a buffer only if a buffer.width> 0 is input)

    CRT <- specify.buffer(trial = CRT, buffer.width = buffer.width)
    CRT <- data.frame_as_CRT(CRT, input.parameters = input.parameters)
    return(CRT)
}


#' Power and sample size calculations for a CRT
#'
#' \code{calculateCRTpower} carries out power and sample-size calculations for a CRT
#' The output is a CRT object containing values for:
#' - The required numbers of clusters to achieve the specified power.
#' - The design effect based on the input ICC.
#' - The power achievable.
#' @param locations numeric: total number of units
#' @param ICC numeric: Intracluster Correlation Coefficient
#' @param effect numeric: Required effect size
#' @param alpha numeric: Significance level
#' @param outcome.type:
#'  \tabular{ll}{
#' \code{"y"}\tab continuous \cr
#' \code{"n"}\tab count \cr
#' \code{"e"}\tab event rate \cr
#' \code{"p"}\tab proportion  \cr
#' \code{"d"}\tab dichotomous \cr
#' }
#' @param sigma2 numeric variance of the outcome (for \code{outcome.type = "y"})
#' @param phi overdispersion parameter (for \code{outcome.type = "n"} or \code{outcome.type = "e"})
#' @param N denominator for proportions (for \code{outcome.type = "p"})
#' @param desiredPower desired power (expressed as a proportion)
#' @param yC expected average outcome in control arm
#' @param mean_h mean number of units per cluster
#' @param sd_h standard deviation of number of units per cluster
#' @return A CRT object comprising the input coordinates :
#' \itemize{
#' \item \code{input.parameters} list of input parameters
#' \item \code{CRT.design.full} list of calculated filename for coordinates of households
#' }
#' @export
#' @details for count data, or event rates a quasi–Poisson model is assumed
#'


#' @examples
#' examplePower = calculateCRTpower(locations = 3000, ICC=0.10, effect=0.4, alpha = 0.05,
#'     outcome.type = 'd', desiredPower = 0.8, yC=0.35, mean_h=100, sd_h=5)
calculateCRTpower <- function(locations, ICC, effect, alpha, outcome.type, sigma2 = NULL,
                              phi = 1, N = 1, desiredPower, yC, mean_h, sd_h) {
    # TODO: add non-binomial outcomes

    # Step A: confidence level Step B: power Step C: Required effect
    # size Step D: ICC, obtained from other studies Step E: baseline
    # outcome Step F: buffer width based on postulated contamination
    # range in km, obtained from other studies Step G\ coordinates of
    # households in study area Step H: proposal for the number of
    # households in each cluster

    # convert power and significance level to normal deviates
    Zsig <- -qnorm(alpha/2)
    Zpow <- qnorm(desiredPower)

    # Power and sample size calculations based on Hemming et al, 2011
    # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102

    link <- switch(outcome.type,
              'y' = "identity",
              "n" = "log",
              "e" = "log",
              'p' = "logit",
              'd' = "logit")
    if (identical(link,'identity')) {
        yI <- yC - effect
        d <- yC - yI
        # with normal models, sigma2 is an input variable
    } else if (identical(link,'log')) {
        yI <- yC * (1 - effect)  # probability in intervened group
        d <- yC - yI  # difference between groups
        phi # Poisson variance is equal to mean.
        sigma2 <- yC * phi
    } else if (identical(link,'logit')) {
        yI <- yC * (1 - effect)  # probability in intervened group
        d <- yC - yI  # difference between groups
        sigma2 <- 1/2 * (yI * (1 - yI) + yC * (1 - yC))
    }
    # required individuals per arm in individually randomized trial
    n_ind <- 2 * sigma2 * ((Zsig + Zpow)/d)^2

    # effective cluster sizes (inflating for multiple observations at the same location)
    mean_eff <- mean_h * N
    sd_eff <- sd_h * N

    # number of clusters assigned to each arm
    k <- floor(locations/(mean_eff * 2))

    # coefficient of variation of the cluster sizes
    cv_eff <- sd_eff/mean_eff

    # design effect (Variance Inflation Factor) allowing for varying
    # cluster sizes(Hemming eqn 6)
    DE <- 1 + (cv_eff^2 + 1) * (mean_eff - 1) * ICC

    # number of individuals required per arm in CRT with equal cluster
    # sizes
    n_crt <- n_ind * DE

    # Step I: Calculations for the required minimum number of clusters
    # for both arms

    # minimum numbers of clusters required assuming varying cluster
    # sizes per arm (Hemming eqn 8)
    min_k <- ceiling(n_ind * (1 + ((cv_eff + 1) * mean_eff - 1) * ICC)/mean_eff)

    clustersRequired <- 2 * min_k

    # power with k clusters per arm power =
    # stats::pnorm(sqrt(k*mean_h/(2*(1+(mean_h-1)*ICC)))*
    # d/sqrt(sigma2) - Zsig) #Hemming eqn 27

    power <- stats::pnorm(sqrt(k * mean_eff/(2 * DE)) * d/sqrt(sigma2) - Zsig)  #unequal cluster sizes

    CRT.design.full <- list(locations = locations, mean_h = mean_h, sd_h = sd_h,
        min_k = min_k, clustersRequired = clustersRequired, DE = DE, power = power)
    input.parameters <- list(locations = locations, ICC = ICC, effect = effect,
        alpha = alpha, outcome.type = outcome.type,  sigma2 = sigma2,
        phi = phi, N = N, desiredPower = desiredPower,
        yC = yC, mean_h = mean_h, sd_h = sd_h)
    CRT <- list(CRT.design.full = CRT.design.full, input.parameters = input.parameters)
    class(CRT) <- "CRT"
    return(CRT)
}
