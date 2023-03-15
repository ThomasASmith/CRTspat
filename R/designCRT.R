#' Design a CRT for a given set of locations
#'
#' \code{designCRT} estimates the required number of clusters and the extent of contamination between arms for a CRT
#' with the input set of locations. Outputs are:
#' (i) Estimates of the required numbers of clusters.
#' (ii) A proposal for the cluster and arm assignments to the input coordinates.
#' (A warning is output if the number of locations is too small to allow randomisation of sufficient clusters).
#' (iii) the proportion of households in the input geography falling within the core of the clusters.
#' @param trial dataframe or \code{"CRT"} object containing Cartesian coordinates of locations in columns 'x' and 'y'.
#' Units are expected to be km.
#' @param alpha numeric: confidence level
#' @param desiredPower numeric: desired power
#' @param effect numeric: required effect size
#' @param yC numeric: baseline outcome
#' @param outcome.type character with options -
#' \code{"y"}: continuous;
#' \code{"n"}: count;
#' \code{"e"}: event rate;
#' \code{"p"}: proportion;
#' \code{"d"}: dichotomous.
#' @param sigma2 numeric: variance of the outcome (for \code{outcome.type = "y"})
#' @param phi numeric: overdispersion parameter (for \code{outcome.type = "n"} or \code{outcome.type = "e"})
#' @param N numeric: denominator for proportions (for \code{outcome.type = "p"})
#' @param ICC numeric: Intra-Cluster Correlation
#' @param buffer.width numeric: required buffer width in km
#' @param h numeric: proposal for the number of coordinates in each cluster
#' @param algorithm algorithm used to determine cluster boundaries, options are
#' \code{"TSP"}: travelling salesman problem heuristic;
#' \code{"NN"}: nearest neighbor;
#' \code{"kmeans"}: kmeans.
#' @return A \code{"CRT"} object comprising the input data, cluster and arm assignments, trial description and results of power calculations
#'  \tabular{lll}{
#'  \code{geom.full}   \tab list \tab summary statistics describing the site,
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
#' exampleDesign = designCRT(trial = readdata('test_site.csv'),
#'                 desiredPower = 0.8, effect=0.4, yC=0.35, outcome.type = 'd', ICC = 0.05,
#'                 buffer.width = 0.25, h=100, algorithm = 'kmeans')
designCRT <- function(trial, alpha = 0.05, desiredPower = 0.8, effect, yC, outcome.type ='d', sigma2 = NULL,
                      phi = 1, N = 1,  ICC, buffer.width = 0, h, algorithm = "kmeans") {

    trial <- CRT_as_data.frame(trial)
    design <- list(alpha = alpha, desiredPower = desiredPower, effect = effect, yC = yC,
         outcome.type = outcome.type, h = h, ICC = ICC, buffer.width = buffer.width, algorithm = algorithm)

    # Step J: specify or compute cluster boundaries
    CRT <- specify.clusters(trial = trial, h = h, algorithm = algorithm, reuseTSP = FALSE)

    # Step K: Random assignment of clusters to arms
    CRT <- randomizeCRT(CRT)

    # augment the trial data frame with distance to nearest discordant
    # coordinate (specifyBuffer assigns a buffer only if a buffer.width> 0 is input)

    CRT <- specify.buffer(trial = CRT, buffer.width = buffer.width)
    CRT <- data.frame_as_CRT(CRT, design = design)
    return(CRT)
}


#' Power and sample size calculations for a CRT
#'
#' \code{calculateCRTpower} carries out power and sample-size calculations for a CRT.
#' The output is an object of class \code{"CRT"} containing values for:
#' - The required numbers of clusters to achieve the specified power.
#' - The design effect based on the input ICC.
#' - The power achievable with clusters of the specified size.
#' @param locations numeric: total number of units available for randomization
#' @param alpha numeric: confidence level
#' @param desiredPower numeric: desired power
#' @param effect numeric: required effect size
#' @param yC numeric: baseline outcome
#' @param outcome.type character with options-
#' \code{"y"}: continuous;
#' \code{"n"}: count;
#' \code{"e"}: event rate;
#' \code{"p"}: proportion;
#' \code{"d"}: dichotomous.
#' @param sigma2 numeric: variance of the outcome (for \code{outcome.type = "y"})
#' @param phi numeric: overdispersion parameter (for \code{outcome.type = "n"} or \code{outcome.type = "e"})
#' @param N numeric: denominator for proportions (for \code{outcome.type = "p"})
#' @param ICC numeric: Intra-Cluster Correlation
#' @param mean_h mean number of units per cluster
#' @param sd_h standard deviation of number of units per cluster
#' @return An object of class \code{"CRT"} comprising:
#' \itemize{
#' \item \code{design} list of input parameters
#' \item \code{geom.full} list containing the calculations of power, required number of clusters and design effect
#' }
#' @export
#' @details
#' Power and sample size calculations are for a two-arm trial, and use the formulae of
#' [Hemming et al, 2011](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102). For counts
#' or event rate data a quasi–Poisson model is assumed.
#'
#' @examples
#' examplePower = calculateCRTpower(locations = 3000, ICC=0.10, effect=0.4, alpha = 0.05,
#'     outcome.type = 'd', desiredPower = 0.8, yC=0.35, mean_h=100, sd_h=5)
calculateCRTpower <- function(locations, alpha, desiredPower, effect, yC, outcome.type, sigma2 = NULL,
                              phi = 1, N = 1, ICC, mean_h, sd_h) {
    # Step A: confidence level Step B: power Step C: Required effect
    # size Step D: ICC, obtained from other studies Step E: baseline
    # outcome Step F: buffer width based on postulated contamination
    # range in km, obtained from other studies Step G\ coordinates of
    # households in study area Step H: proposal for the number of
    # households in each cluster

    # convert power and significance level to normal deviates
    Zsig <- -qnorm(alpha/2)
    Zpow <- qnorm(desiredPower)

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

    geom.full <- list(locations = locations, mean_h = mean_h, sd_h = sd_h,
        min_k = min_k, clustersRequired = clustersRequired, DE = DE, power = power)
    design <- list(locations = locations, alpha = alpha, desiredPower = desiredPower, effect = effect,
                             yC = yC, outcome.type = outcome.type,  sigma2 = sigma2,
                             phi = phi, N = N, ICC = ICC, mean_h = mean_h, sd_h = sd_h)
    CRT <- list(geom.full = geom.full, design = design)
    class(CRT) <- "CRT"
    return(CRT)
}
