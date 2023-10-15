#' Power and sample size calculations for a CRT
#'
#' \code{CRTpower} carries out power and sample size calculations for CRTs.
#'
#' @param trial dataframe or \code{'CRTsp'} object containing Cartesian coordinates of locations in columns 'x' and 'y'.
#' @param locations numeric: total number of units available for randomization (required if \code{trial} is not specified)
#' @param alpha numeric: confidence level
#' @param desiredPower numeric: desired power
#' @param effect numeric: required effect size
#' @param yC numeric: baseline value of outcome
#' @param outcome_type character with options -
#' \code{'y'}: continuous;
#' \code{'n'}: count;
#' \code{'e'}: event rate;
#' \code{'p'}: proportion;
#' \code{'d'}: dichotomous.
#' @param sigma2 numeric: variance of the outcome (required for \code{outcome_type = 'y'})
#' @param phi numeric: overdispersion parameter (for \code{outcome_type = 'n'} or \code{outcome_type = 'e'})
#' @param N numeric: mean of the denominator for proportions (for \code{outcome_type = 'p'})
#' @param ICC numeric: Intra-Cluster Correlation
#' @param k integer: number of clusters in each arm (required if \code{trial} is not specified)
#' @param sd_h standard deviation of number of units per cluster (required if \code{trial} is not specified)
#' @returns A list of class \code{'CRTsp'} object comprising the input data, cluster and arm assignments,
#' trial description and results of power calculations
#' @export
#' @details
#' Power and sample size calculations are for a two-arm trial using the formulae of
#' [Hemming et al, 2011](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102) which use a
#' normal approximation for the inter-cluster variation. For counts
#' or event rate data a quasi–Poisson model is assumed. The functions do not consider any loss
#' in power due to contamination, loss to follow-up etc.
#' If geolocations are not input power and sample size calculations are based on the scalar input parameters.\cr\cr
#' If a trial dataframe or \code{'CRTsp'} object containing a pre-existing randomization is input then the
#' numbers and sizes of clusters are in the input data are used to estimate the power. If buffer zones have been specified
#' then separate calculations are made for the core area and for the full site\cr
#' The output is an object of class \code{'CRTsp'} containing any input trial data.frame and values for:
#' - The required numbers of clusters to achieve the specified power.
#' - The design effect based on the input ICC.
#' - Calculations of the nominal power (ignoring any bias caused by contamination effects)\cr
#' @examples
#' {# Example without input geolocations
#' examplePower1 = CRTpower(locations = 3000, ICC=0.10, effect=0.4, alpha = 0.05,
#'     outcome_type = 'd', desiredPower = 0.8, yC=0.35, k = 20, sd_h=5)
#' summary(examplePower1)
#' # Example with input geolocations and randomisation
#' examplePower2 = CRTpower(trial = readdata('example_site.csv'), desiredPower = 0.8,
#'     effect=0.4, yC=0.35, outcome_type = 'd', ICC = 0.05, k = 20)
#' summary(examplePower2)
#' }
CRTpower <- function(trial = NULL, locations = NULL, alpha = 0.05, desiredPower = 0.8,
    effect = NULL, yC = NULL, outcome_type = "d", sigma2 = NULL, phi = 1,
    N = 1, ICC = NULL, k = NULL, sd_h = 0) {

    CRT <- CRTsp(trial)

    # populate a design list with a data about the input trial (if
    # available) and the input parameters

    # sigma2 is only required for continuous data, so its absence
    # should otherwise not crash the program
    if (!identical(outcome_type, "y")) sigma2 <- NA

    design <- ifelse(is.null(CRT$design$locations), list(), CRT$design)
    design$locations <- ifelse((nrow(CRT$trial) == 0), locations, nrow(CRT$trial))
    parnames <- c("alpha", "desiredPower", "effect", "yC", "outcome_type",
        "sigma2", "phi", "N", "ICC", "k", "sd_h")
    # Identify which variables to retrieve from the pre-existing design
    from_old <- lapply(mget(parnames), FUN = is.null)
    design[parnames] <- ifelse(from_old, design[parnames], mget(parnames))
    missing <- lapply(design[parnames], FUN = is.null)
    if (TRUE %in% missing) {
        stop("*** Value(s) must be supplied for: ", parnames[missing == TRUE],
            " ***")
    }
    CRT <- CRTsp(CRT, design = design)
}


# Characteristics of a trial design. The input is a data frame or
# CRTsp object. The output list conforms to the requirements for a
# CRTsp object
get_geom <- function(trial = NULL, design = NULL) {

    sd_distance <- clustersRequired <- DE <- power <- NULL
    if(!is.null(design)) {
        locations <- design$locations
        sd_h <- design$sd_h
        k <- design$k
        mean_h = locations/(2 * k)
    } else {
        mean_h <- k <- sd_h <- locations <- NULL
    }
    geom <- list(locations = locations,
                 sd_h = sd_h,
                 k= k,
                 records = 0,
                 mean_h = mean_h,
                 DE = NULL,
                 power = NULL,
                 clustersRequired = NULL)

    # overwrite values from the design with those from the data frame if these are present
    if (!is.null(trial) & nrow(trial) > 0) {
        coordinates <- data.frame(cbind(x = trial$x, y = trial$y))
        geom$records <- nrow(trial)
        geom$locations <- nrow(dplyr::distinct(coordinates))

        if (!is.null(trial$cluster)) {
            # reassign the cluster levels in case some are not
            # represented in this geom (otherwise nlevels() counts
            # clusters that are not present)
            trial$cluster <- as.factor(as.character(trial$cluster))

            geom$k <- floor(nlevels(trial$cluster)/2)
            # mean number of locations randomized in each cluster
            geom$mean_h <- mean(table(trial$cluster))

            # standard deviation of locations assigned to each
            # cluster
            geom$sd_h <- stats::sd(table(trial$cluster))

        }
        if (!is.null(trial$arm)) {
            geom$sd_distance <- stats::sd(trial$nearestDiscord)
            arms <- unique(cbind(trial$cluster, trial$arm))[, 2]  #assignments
        }
    }
    # cluster size
    geom$k <- ifelse(is.null(geom$k), NA, round(geom$k))
    geom$mean_h <- geom$locations/(2 * geom$k)

    if (!is.null(design$effect)) {
        if (is.null(geom$locations)) {
            stop("*** Number of locations is a required input ***")
        }
        if (is.null(geom$k)) {
            stop("*** Number of clusters is a required input ***")
        }
        if (identical(geom$sd_h, 0)) {
            message("*** Assuming all clusters are the same size ***")
        }
        effect <- design$effect
        yC <- design$yC


        # Step A: confidence level Step B: power Step C: Required
        # effect size Step D: ICC, obtained from other studies Step E:
        # baseline outcome Step F: buffer width based on postulated
        # contamination range in km, obtained from other studies Step
        # G\ coordinates of households in study area Step H: proposal
        # for the number of households in each cluster

        link <- switch(design$outcome_type, y = "identity", n = "log", e = "log",
            p = "logit", d = "logit")
        if (identical(link, "identity")) {
            yI <- yC - design$effect
            d <- yC - yI
            sigma2 <- design$sigma2
            # with normal models, sigma2 is an input variable
        } else if (identical(link, "log")) {
            yI <- yC * (1 - design$effect)  # probability in intervened group
            d <- yC - yI  # difference between groups
            # Poisson variance is equal to mean.
            sigma2 <- yC * design$phi
        } else if (identical(link, "logit")) {
            yI <- yC * (1 - design$effect)  # probability in intervened group
            d <- yC - yI  # difference between groups
            # This is the variance for a Bernoulli. The cluster sizes are
            # inflated for the binomial case (below)
            sigma2 <- 1/2 * (yI * (1 - yI) + yC * (1 - yC))
        }


        # convert power and significance level to Zvalues
        Zsig <- -qnorm(design$alpha/2)
        Zpow <- qnorm(design$desiredPower)

        # corresponding t values (not used) df <- k - 1 t_sig <- - qt(p
        #                               = design$alpha/2, df = df) t_pow <- qt(p =
        #                               design$desiredPower, df = df)

        # required individuals per arm in individually randomized trial
        n_ind <- 2 * sigma2 * ((Zsig + Zpow)/d)^2

        # effective cluster sizes (inflating for multiple observations
        # at the same location)
        mean_eff <- geom$mean_h * design$N
        sd_eff <- geom$sd_h * design$N

        # coefficient of variation of the cluster sizes
        cv_eff <- sd_eff/mean_eff

        # design effect (Variance Inflation Factor) allowing for
        # varying cluster sizes(Hemming eqn 6)
        geom$DE <- 1 + (cv_eff^2 + 1) * (mean_eff - 1) * design$ICC

        # number of individuals required per arm in CRT with equal
        # cluster sizes
        n_crt <- n_ind * geom$DE

        # Step I: Calculations for the required minimum number of
        # clusters for both arms

        # minimum total numbers of clusters required assuming varying
        # cluster sizes per arm (Hemming eqn 8)
        geom$clustersRequired <- ceiling(2 * (n_ind * (1 + ((cv_eff + 1) *
            mean_eff - 1) * design$ICC)/mean_eff))

        # power with k clusters per arm and unequal cluster sizes
        geom$power <- with(geom, stats::pnorm(
            sqrt(k * mean_eff/(2 * geom$DE)) * d/sqrt(sigma2) - Zsig))

    }
    return(geom)
}

