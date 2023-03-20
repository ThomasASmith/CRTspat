#' Create co-ordinates for a simulated CRT
#'
#' \code{simulateSite} creates a set of Cartesian co-ordinates for use as the locations in a simulated trial site
#' @param geoscale standard deviation of random displacement from each settlement cluster center
#' @param locations number of locations in population
#' @param kappa intensity of Poisson process of settlement cluster centers
#' @param mu mean  number of points per settlement cluster
#' @returns A list of class \code{"CRTspat"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site\tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  }
#' @details \code{simulateSite} simulates a human settlement pattern using the Thomas algorithm (\code{rThomas} function
#' in [\code{spatstat}](http://spatstat.org/) allowing the user to defined the density of locations and degree of spatial clustering.
#' The results are output as a set of Cartesian coordinates centred at the origin.
#' @export
#' @examples
#' # Generate a simulated area with 10,000 locations
#' example_area = simulateSite(geoscale = 1, locations=10000, kappa=3, mu=40)
simulateSite <- function(geoscale, locations, kappa, mu) {
    scaling = geoscale * 10
    # Poisson point pattern with Thomas algorithm
    p <- spatstat.random::rThomas(kappa, geoscale, mu, win = spatstat.geom::owin(c(0, scaling), c(0, scaling)))
    # expected number of points: kappa*mu*scaling^2

    # create locations and specify co-ordinates
    hhID <- c(1:locations)
    x <- p$x[seq(1:locations)]
    y <- p$y[seq(1:locations)]
    coordinates <- data.frame(x = x - mean(x), y = y - mean(y))
    CRT <- as_CRTspat(coordinates, design = NULL)
    return(CRT)
}


#' Simulation of cluster randomized trial with contamination
#'
#' \code{simulateCRT} generates simulated data for a cluster randomized trial (CRT) with geographic contamination between arms.
#'
#' @param trial an object of class \code{"CRTspat"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}). Each location may also be
#'   assigned a \code{propensity} (see details).
#' @param effect numeric. The simulated effect size (defaults to 0)
#' @param outcome0 numeric. The anticipated value of the outcome in the absence of intervention
#' @param generateBaseline logical. If \code{TRUE} then baseline data will be simulated
#' @param matchedPair logical. If \code{TRUE} then the function tries to carry out randomization
#' using pair-matching on the baseline data (see details)
#' @param scale measurement scale of the outcome. Options are: 'proportion' (the default); 'count'; 'continuous'.
#' @param baselineNumerator optional name of numerator variable for baseline data
#' @param baselineDenominator optional name of denominator variable for baseline data
#' @param denominator optional name of denominator variable for the outcome
#' @param ICC_inp numeric. Target intra cluster correlation, provided as input when baseline data are to be simulated
#' @param sd numeric. standard deviation of the normal kernel measuring spatial smoothing leading to contamination
#' @param theta_inp numeric. input contamination range
#' @param tol numeric. tolerance of output ICC
#' @param vr numeric. ratio of location variance to cluster variance (for continuous outcomes)
#' @returns A list of class \code{"CRTspat"} containing the following components:
#' \tabular{llll}{
#' \code{geom.full}\tab list: \tab summary statistics describing the site
#' cluster assignments, and randomization \tab\cr
#' \code{design}\tab list: \tab values of input parameters to the design \tab\cr
#' \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#' \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#' \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#' \tab\code{cluster} \tab factor: \tab assignments to cluster of each location  \cr
#' \tab\code{arm} \tab factor: \tab assignments to \code{control} or \code{intervention} for each location \cr
#' \tab\code{nearestDiscord} \tab numeric vector: \tab Euclidean distance to nearest discordant location (km) \cr
#' \tab\code{propensity} \tab numeric vector: \tab propensity of location \cr
#' \tab\code{base_denom} \tab numeric vector: \tab denominator for baseline \cr
#' \tab\code{base_num} \tab numeric vector: \tab numerator for baseline \cr
#' \tab\code{denom} \tab numeric vector: \tab denominator for the outcome \cr
#' \tab\code{num} \tab numeric vector: \tab numerator for the outcome \cr
#' \tab\code{...} \tab\tab other objects included in the input \code{"CRTspat"} object
#' or \code{data.frame}\cr
#' }
#' @details Synthetic data are generated by sampling around the values of
#' variable \code{propensity}, which is a numerical vector
#' (taking positive values) of length equal to the number of locations.
#' There are three ways in which \code{propensity} can arise:
#' \enumerate{
#' \item \code{propensity} can be provided as part of the input \code{trial} object.
#' \item Baseline numerators and denominators (values of \code{baselineNumerator}
#' and \code{baselineDenominator} may be provided.
#' \code{propensity} is then generated as the numerator:denominator ratio
#' for each location in the input object
#' \item Otherwise \code{propensity} is generated using a 2D Normal
#' kernel density. The [\code{OOR::StoSOO}](https://rdrr.io/cran/OOR/man/StoSOO.html)
#' is used to achieve an intra-cluster correlation coefficient (ICC) that approximates
#' the value of \code{'ICC_inp'} by searching for an appropriate value of the kernel bandwidth.
#' }
#' \code{num[i]}, the synthetic outcome for location \code{i}
#' is simulated with expectation: \cr
#' \deqn{\code{E(num[i]) = outcome0[i] * propensity[i] * denom[i] * (1 - effect*I[i])/mean(outcome0[] * propensity[])}} \cr
#' The sampling distribution of \code{num[i]} depends on the value of \code{scale} as follows: \cr
#' \tabular{llll}{
#' \code{scale}=’continuous’ \tab Values of \code{num} are sampled from a
#' Normal distributions with means \code{E(num[i])}
#' and variance determined by \code{vr} and the fitting to \code{ICC_inp}.\cr
#' \code{scale}=’count’ \tab Simulated events are allocated to locations via multivariate hypergeometric distributions
#' parameterised with \code{E(num[i])}.\cr
#' \code{scale}=’proportion’\tab Simulated events are allocated to locations via multinomial distributions
#' parameterised with \code{E(num[i])}.\cr
#' }
#' \code{denominator} may specify a vector of numeric (non-zero) values
#' in the input \code{"CRTspat"} or \code{data.frame} which is returned
#' as variable \code{denom}. It acts as a scale-factor for continuous outcomes, rate-multiplier
#' for counts, or denominator for proportions. For discrete data all values of \code{denom}
#' must be > 0.5 and are rounded to the nearest integer in calculations of \code{num}.\cr\cr
#' By default, \code{denom} is generated as a vector of ones, leading to simulation of
#' dichotomous outcomes if \code{scale}=’proportion’.\cr
#'
#' If baseline numerators and denominators are provided then the output vectors
#' \code{base_denom} and  \code{base_num} are set to the input values. If baseline numerators and denominators
#' are not provided then the synthetic baseline data are generated by sampling around \code{propensity} in the same
#' way as the outcome data, but with the effect size set to zero.
#'
#' If \code{matchedPair} is \code{TRUE} then pair-matching on the baseline data will be used in randomization providing
#' there are an even number of clusters. If there are an odd number of clusters then matched pairs are not generated and
#' an unmatched randomization is output.
#'
#' Either \code{sd} or \code{theta_inp} must be provided. If both are provided then
#' the value of \code{sd} is overwritten
#' by the standard deviation implicit in the value of \code{theta_inp}.
#' Contamination is simulated as arising from a diffusion-like process.
#'
#' For further details see [Multerer (2021)](https://edoc.unibas.ch/85228/)
#' @export
#'
#' @examples
#' example_simulated_CRT =  simulateCRT(trial=readdata('test_Arms.csv'),
#'                                      effect=0.25,
#'                                      ICC_inp=0.05,
#'                                      outcome0=0.5,
#'                                      matchedPair = FALSE,
#'                                      scale='proportion',
#'                                      sd=0.6,
#'                                      tol=0.05)
simulateCRT <- function(trial = NULL, effect = 0, outcome0 = NULL, generateBaseline = TRUE, matchedPair = TRUE, scale = "proportion",
    baselineNumerator = "base_num", baselineDenominator = "base_denom", denominator = NULL, ICC_inp = NULL, sd = NULL, theta_inp = NULL,
    tol = 1e-04, vr = 0.5) {

    # Written by Tom Smith, July 2017. Adapted by Lea Multerer, September 2017
    cat("\n=====================    SIMULATION OF CLUSTER RANDOMISED TRIAL    =================\n")
    bw <- NULL
    if (identical(class(trial), "data.frame")) {
        CRT <- list(trial = trial, design = NULL)
        class(CRT) <- "CRTspat"
    } else {
        CRT <- trial
        trial <- CRT$trial
    }

    trial$arm <- as.factor(trial$arm)
    trial$cluster <- as.factor(trial$cluster)

    # set the denominator variable to be 'denom'
    if (is.null(denominator))
        denominator <- "denom"
    trial$denom <- trial[[denominator]]
    if (denominator != "denom")
        trial[[denominator]] <- NULL

    # use contamination range if this is available
    if (!is.null(theta_inp)) {
        sd <- theta_inp/(sqrt(2) * qnorm(0.95))
    }
    if (is.null(sd)) {
        print("Error: contamination range or s.d. of spatial kernel must be provided")
        trial <- NULL
        return(trial)
    }
    # trial needs to be ordered for GEE analyses (estimation of ICC)
    trial <- trial[order(trial$cluster), ]

    # For the smoothing step compute contributions to the relative effect size from other locations as a function of
    # distance to the other locations

    euclid <- distance_matrix(trial$x, trial$y)

    # generate baseline data if required and exposure proxy if this is not provided

    if (!"propensity" %in% colnames(trial) & baselineNumerator %in% colnames(trial) & baselineDenominator %in% colnames(trial)) {
        trial$propensity <- trial[[baselineNumerator]]/trial[[baselineDenominator]]
        # TODO: perhaps replace this with estimation via INLA with a spatial model

    } else if ("propensity" %in% colnames(trial)) {
        # create a baseline dataset using a pre-existing exposure proxy
        trial <- syntheticBaseline(bw = NULL, trial = trial, sd = sd, euclid = euclid, outcome0 = outcome0)

    } else if (generateBaseline) {

        # determine the required smoothing bandwidth by fitting to the pre-specified ICC
        # compute approximate diagonal of clusters

        approx_diag <- sqrt((max(trial$x) - min(trial$x))^2 + (max(trial$y) - min(trial$y))^2)/sqrt(length(unique(trial$cluster)))
        cat("Estimating the smoothing required to achieve the target ICC of", ICC_inp, "\n")

        loss <- 999
        nb_iter <- 20
        if (identical(Sys.getenv("TESTTHAT"), "true")) nb_iter <- 5
        # in testing, the number of iterations is reduced giving very approximate output
        while (loss > tol) {
            ICC.loss <- OOR::StoSOO(par = NA, fn = ICCdeviation, lower = -5, upper = 5, nb_iter = nb_iter,
                                    trial = trial, ICC_inp = ICC_inp, approx_diag = approx_diag, sd = sd,
                                    scale = scale, euclid = euclid, effect = effect, outcome0 = outcome0)
            loss <- ICC.loss$value
        }
        logbw <- ICC.loss$par
        # overprint the output that was recording progress
        cat("\r                                                         \n")

        # set the seed so that the same result is obtained for a specific bandwidth
        bw <- exp(logbw)
        set.seed(round(bw * 1e+06))
        # create a baseline dataset using the optimized bandwidth
        trial <- syntheticBaseline(bw = bw, trial = trial, sd = sd, euclid = euclid, outcome0 = outcome0)

    }

    trial <- get_assignments(trial = trial, scale = scale, euclid = euclid, sd = sd, effect = effect,
                             outcome0 = outcome0, denominator = denominator)
    CRT <- updateCRT(CRT = CRT, trial = trial)
    return(CRT)
}

# Assign expected outcome to each location assuming a fixed effect size.
get_assignments <- function(trial, scale, euclid, sd, effect, outcome0, denominator) {

    # Indicator of whether the source is intervened is (as.numeric(trial$arm[i]) - 1 smoothedIntervened is the value of
    # propensity decremented by the effect of intervention and smoothed to allow for mosquito movement

    if (sd > 0) {
        smoothedIntervened <- gauss(sd, euclid) %*% (trial$propensity * (1 - effect * (as.numeric(trial$arm) - 1)))
    } else {
        smoothedIntervened <- trial$propensity * (1 - effect * (as.numeric(trial$arm) - 1))
    }
    # distances to nearest discordant locations
    discord <- outer(trial$arm, trial$arm, "!=")  #returns true & false.
    euclidd <- ifelse(discord, euclid, 99999.9)

    # for the control arm return the minimal distance with a minus sign
    trial$nearestDiscord <- ifelse(trial$arm == "control", -apply(euclidd, MARGIN = 2, min), apply(euclidd, MARGIN = 2, min))

    trial <- distributePositives(trial = trial, outcome0 = outcome0, smoothed = smoothedIntervened, effect = effect, scale = scale,
        denominator = "denom", numerator = "num")
    return(trial)
}

# allocate the positives to locations
distributePositives <- function(trial, outcome0, smoothed, scale, effect, denominator, numerator) {
    expected_ratio <- num <- rowno <- sumnum <- NULL
    if (!(denominator %in% colnames(trial)))
        trial[[denominator]] <- 1

    # the denominator must be an integer; this changes the value if a non-integral value is input
    trial[[denominator]] <- round(trial[[denominator]], digits = 0)

    # compute the total positives expected given the input effect size
    npositives <- round(outcome0 * sum(trial[[denominator]]) * (1 - 0.5 * effect))

    # scale to input value of initial prevalence by assigning required number of infections with probabilities proportional
    # to smoothedIntervened multiplied by the denominator

    expected_allocation <- smoothed * trial[[denominator]]/sum(smoothed * trial[[denominator]])

    trial$expected_ratio <- expected_allocation/trial[[denominator]]
    trial$rowno <- seq(1:nrow(trial))

    # expand the vector of locations to allow for denominators > 1
    triallong <- trial %>%
        tidyr::uncount(trial[[denominator]])

    # To generate count data, records in triallong can be sampled multiple times. To generate proportions each record can
    # only be sampled once.
    replacement <- identical(scale, "log")

    # sample generates a multinomial sample and outputs the indices of the locations assigned
    positives <- sample(x = nrow(triallong), size = npositives, replace = FALSE, prob = triallong$expected_ratio)
    triallong$num <- 0
    triallong$num[positives] <- 1

    # summarise the numerator values into the original set of locations
    numdf <- dplyr::group_by(triallong, rowno) %>%
        dplyr::summarise(sumnum = sum(num))
    numdf[[numerator]] <- numdf$sumnum

    # remove any superseded numerator variable
    trial[[numerator]] <- NULL

    # use left_join to merge into the original data frame (records with zero denominator do not appear in numdf)
    trial <- trial %>%
        dplyr::left_join(numdf, by = "rowno")

    # remove temporary variables and replace missing numerators with zero (because the multinomial sampling algorithm leaves
    # NA values where no events are assigned)
    trial <- subset(trial, select = -c(rowno, expected_ratio, sumnum))
    if (sum(is.na(trial[[numerator]])) > 0) {
        cat("** Warning: some records have zero denominator after rounding **\n")
        cat("You may want to remove these records or rescale the denominators \n")
        trial[[numerator]][is.na(trial[[numerator]])] <- 0
    }
    return(trial)
}

# deviation of ICC from target as a function of bandwidth
ICCdeviation <- function(logbw, trial, ICC_inp, approx_diag, sd, scale, euclid, effect, outcome0) {
    cluster <- NULL
    # set the seed so that a reproducible result is obtained for a specific bandwidth
    if (!is.null(logbw)) {
        bw <- exp(logbw)
        set.seed(round(bw * 1e+06))
    }

    trial <- syntheticBaseline(bw = bw, trial = trial, sd = sd, euclid = euclid, outcome0 = outcome0)
    trial <- get_assignments(trial = trial, scale = scale, euclid = euclid, sd = sd, effect = effect,
                             outcome0 = outcome0, denominator = "denom")
    trial$neg <- trial$denom - trial$num
    fit <- geepack::geeglm(cbind(num, neg) ~ arm, id = cluster, corstr = "exchangeable", data = trial, family = binomial(link = "logit"))
    summary_fit <- summary(fit)
    # Intracluster correlation
    ICC <- noLabels(summary_fit$corr[1])  #with corstr = 'exchangeable', alpha is the ICC
    cat("\rbandwidth: ", bw, "  ICC=", ICC, "        \r")
    loss <- (ICC - ICC_inp)^2
    return(loss)
}

# assign initial pattern if it does not exist
syntheticBaseline <- function(bw, trial, sd, euclid, outcome0) {
    if (!is.null(bw)) {
        trial$propensity <- KDESmoother(trial$x, trial$y, kernnumber = 200, bandwidth = bw, low = 0, high = 1)
    }
    # Smooth the exposure proxy to allow for mosquito movement.  the s.d. in each dimension of the 2 d gaussian is
    # sd/sqrt(2) smoothedBaseline is the amount received by the each cluster from the contributions (propensity) of each
    # source
    if (sd > 0) {
        smoothedBaseline <- gauss(sd, euclid) %*% trial$propensity
    } else {
        smoothedBaseline <- trial$propensity
    }

    trial <- distributePositives(trial = trial, outcome0 = outcome0, smoothed = smoothedBaseline, effect = 0, scale = scale,
        denominator = "base_denom", numerator = "base_num")
    return(trial)
}

# compute a euclidian distance matrix
distance_matrix <- function(x, y) {
    # generates square matrices of differences
    xdist <- outer(x, x, "-")
    ydist <- outer(y, y, "-")
    euclid <- sqrt(xdist * xdist + ydist * ydist)
}


# add lognormal noise: not sure this function is needed X is the input vector comprising a sample from a smoothed
# distribution varXY is the required variance
add_noise <- function(X, varXY) {
    muY <- 1
    varY <- (varXY - var(X))/(1 + var(X))
    mu <- log(muY/sqrt(1 + varY/muY))
    var <- log(1 + varY/muY)
    Y <- stats::rlnorm(length(XY), meanlog = mu, sdlog = sqrt(var))
    XY <- X * Y
    return(XY)
}


# contribution of i to j as a function of the Gaussian process used in simulating contamination
gauss <- function(sd, euclid) {
    # definition of a gauss function
    f <- (1/(2 * pi * sd^2)) * exp(-(euclid^2)/(2 * (sd^2)))
    totalf <- rowSums(f)  #sums over rows of matrix f
    # Careful here, totalf sums over very small numbers, consider replacing
    return(f/totalf)
}

# generate a random pattern of vectorial capacity with smoothing
KDESmoother <- function(x, y, kernnumber, bandwidth, low, high) {

    # force bandwidth to be scalar
    bandwidth <- bandwidth[1]

    sam <- sample(1:length(x), kernnumber, replace = F)

    xdist <- outer(x, x[sam], "-")
    ydist <- outer(y, y[sam], "-")
    euclid <- sqrt(xdist * xdist + ydist * ydist)  #is not a square matrix

    f <- (1/(2 * pi * bandwidth^2)) * exp(-(euclid^2)/(2 * (bandwidth^2)))
    totalf <- (1/ncol(euclid)) * rowSums(f)  #sums over rows of matrix f

    smoother <- low + totalf * ((high - low)/max(totalf))

    return(smoother)
}

