#' Create co-ordinates for a simulated CRT
#'
#' \code{simulateSite} creates a set of Cartesian co-ordinates for use as the locations in a simulated trial site
#' @param scale standard deviation of random displacement from each settlement cluster center
#' @param locations number of locations in population
#' @param kappa intensity of Poisson process of settlement cluster centers
#' @param mu mean  number of points per settlement cluster
#' @returns A list with class \code{CRT} containing the following components:
#'  \tabular{llll}{
#'  \code{CRT.design.full}:   \tab list \tab summary statistics describing the site\cr
#'  \code{x}: \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y}: \tab numeric vector \tab y-coordinates of locations \cr
#' }
#' @details \code{simulateSite} simulates a human settlement pattern using the Thomas algorithm (\code{rThomas} function
#' in [\code{spatstat}](http://spatstat.org/) allowing the user to defined the density of locations and degree of spatial clustering.
#' The results are output as a set of Cartesian coordinates centred at the origin.
#' @export
#' @examples
#' # Generate a simulated area with 10,000 locations
#' example_area = simulateSite(scale = 2, locations=10000, kappa=3, mu=40)
simulateSite <- function(scale, locations, kappa, mu) {
    scaling = scale * 20
    # Poisson point pattern with Thomas algorithm
    p <- spatstat.random::rThomas(kappa, scale, mu, win = spatstat.geom::owin(c(0, scaling), c(0, scaling)))
    # expected number of points: kappa*mu*scaling^2

    # create locations and specify co-ordinates
    hhID <- c(1:locations)
    x <- p$x[seq(1:locations)]
    y <- p$y[seq(1:locations)]
    coordinates <- data.frame(x = x, y = y)
    CRT <- convert.data.frame.CRT(coordinates, input.parameters = NULL)
    return(CRT)
}


#' Simulation of cluster randomized trial with contamination
#'
#' \code{simulateCRT} generates simulated data for a cluster randomized trial (CRT) with geographic contamination between arms and a binary outcome. Contamination is simulated as arising from a diffusion-like process.
#'
#' @param trial a \code{CRT} object or \code{data.frame} containing locations in (x,y) coordinates; cluster
#'   assignments (factor \code{cluster}); and arm assignments (factor \code{arm})
#' @param efficacy simulated efficacy (defaults to 0)
#' @param initialPrevalence prevalence in control arm (assumed equal to initial proportion)
#' @param generateBaseline logical indicator of whether baseline data should be simulated
#' @param matchedPair logical: indicator of whether pair-matching on the baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (if present)
#' @param baselineDenominator name of denominator variable for baseline data (if present)
#' @param denominator name of denominator variable for the outcome (if present)
#' @param ICC_inp Intra Cluster Correlation, provided as input when baseline data are to be simulated
#' @param sd standard deviation of the normal kernel measuring spatial smoothing leading to contamination
#' @param theta_inp input contamination range
#' @param tol tolerance of output ICC
#' @returns A list of class \code{CRT} containing the following components:
#' \tabular{llll}{
#' \code{input.parameters}:   \tab list \tab values of input parameters to the function \cr
#' \code{CRT.design.full}:   \tab list \tab summary statistics describing the site,
#'  cluster assignments, and randomization \cr
#' \code{x}: \tab numeric vector \tab x-coordinates of locations \cr
#' \code{y}: \tab numeric vector \tab y-coordinates of locations \cr
#' \code{cluster}: \tab factor \tab assignments to cluster of each location  \cr
#' \code{arm}: \tab factor \tab assignments to \code{control} or \code{intervention} for each location \cr
#' \code{nearestDiscord}: \tab numeric vector \tab Euclidean distance to nearest discordant location (km) \cr
#' \code{propensity}: propensity of location \cr
#' \code{base_denom}: denominator for baseline \cr
#' \code{base_num}: numerator for baseline \cr
#' \code{denom}: denominator for the outcome \cr
#' \code{num}: numerator for the outcome \cr
#' \code{...};   \tab other objects included in the input \code{CRT} object or \code{data.frame}  \cr
#' }
#' @details For details see [Multerer (2021)](https://edoc.unibas.ch/85228/)
#' @export
#'
#' @examples
#' example_simulated_CRT =  simulateCRT(trial=readdata('test_Arms.csv'),
#'                                      efficacy=0.25,
#'                                      ICC_inp=0.05,
#'                                      initialPrevalence=0.5,
#'                                      matchedPair = TRUE,
#'                                      sd=0.6,
#'                                      tol=0.05)
simulateCRT <- function(trial = NULL, efficacy = 0, initialPrevalence = NULL, generateBaseline = TRUE, matchedPair = TRUE,
    baselineNumerator = "base_num", baselineDenominator = "base_denom", denominator = NULL, ICC_inp = NULL, sd = NULL, theta_inp = NULL,
    tol = 1e-04) {

    # Written by Tom Smith, July 2017. Adapted by Lea Multerer, September 2017
    cat("\n=====================    SIMULATION OF CLUSTER RANDOMISED TRIAL    =================\n")
    bw <- NULL
    # several operations require the input data as data.frame. Descriptors will be replaced
    trial <- convertCRT.data.frame(CRT = trial)

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
        # TODO: replace this with estimation via INLA with a spatial model

    } else if ("propensity" %in% colnames(trial)) {
        # create a baseline dataset using a pre-existing exposure proxy
        trial <- syntheticBaseline(bw = NULL, trial = trial, sd = sd, euclid = euclid, initialPrevalence = initialPrevalence)

    } else if (generateBaseline) {

        # determine the required smoothing bandwidth by fitting to the pre-specified ICC compute approximate diagonal of
        # clusters
        approx_diag <- sqrt((max(trial$x) - min(trial$x))^2 + (max(trial$y) - min(trial$y))^2)/sqrt(length(unique(trial$cluster)))
        cat("Estimating the smoothing required to achieve the target ICC of", ICC_inp, "\n")

        loss <- 999
        while (loss > tol) {
            ICC.loss <- OOR::StoSOO(par = NA, fn = ICCdeviation, lower = -5, upper = 5, nb_iter = 20, trial = trial, ICC_inp = ICC_inp,
                approx_diag = approx_diag, sd = sd, euclid = euclid, efficacy = efficacy, initialPrevalence = initialPrevalence)
            loss <- ICC.loss$value
        }
        logbw <- ICC.loss$par
        # overprint the output that was recording progress
        cat("\r                                                         \n")

        # set the seed so that the same result is obtained for a specific bandwidth
        bw <- exp(logbw)
        set.seed(round(bw * 1e+06))
        # create a baseline dataset using the optimized bandwidth
        trial <- syntheticBaseline(bw = bw, trial = trial, sd = sd, euclid = euclid, initialPrevalence = initialPrevalence)

    }

    trial <- assignPositives(trial = trial, euclid = euclid, sd = sd, efficacy = efficacy, initialPrevalence = initialPrevalence,
        denominator = denominator)
    class(trial) <- "CRT"
    return(trial)
}

# Assign expected proportions to each location assuming a fixed efficacy.
assignPositives <- function(trial, euclid, sd, efficacy, initialPrevalence, denominator) {

    # Indicator of whether the source is intervened is (as.numeric(trial$arm[i]) - 1 smoothedIntervened is the value of
    # propensity decremented by the effect of intervention and smoothed to allow for mosquito movement

    if (sd > 0) {
        smoothedIntervened <- gauss(sd, euclid) %*% (trial$propensity * (1 - efficacy * (as.numeric(trial$arm) - 1)))
    } else {
        smoothedIntervened <- trial$propensity * (1 - efficacy * (as.numeric(trial$arm) - 1))
    }
    # distances to nearest discordant locations
    discord <- outer(trial$arm, trial$arm, "!=")  #returns true & false.
    euclidd <- ifelse(discord, euclid, 99999.9)

    # for the control arm return the minimal distance with a minus sign
    trial$nearestDiscord <- ifelse(trial$arm == "control", -apply(euclidd, MARGIN = 2, min), apply(euclidd, MARGIN = 2, min))

    trial <- distributePositives(trial = trial, initialPrevalence = initialPrevalence, smoothed = smoothedIntervened, efficacy = efficacy,
        denominator = "denom", numerator = "num")
    return(trial)
}

# allocate the positives to locations
distributePositives <- function(trial, initialPrevalence, smoothed, efficacy, denominator, numerator) {
    expected_proportion <- num <- rowno <- sumnum <- NULL
    if (!(denominator %in% colnames(trial)))
        trial[[denominator]] <- 1

    # the denominator must be an integer; this changes the value if a non-integral value is input
    trial[[denominator]] <- round(trial[[denominator]], digits = 0)

    # compute the total positives expected given the input efficacy
    npositives <- round(initialPrevalence * sum(trial[[denominator]]) * (1 - 0.5 * efficacy))

    # scale to input value of initial prevalence by assigning required number of infections with probabilities
    # proportional to smoothedIntervened multiplied by the denominator

    expected_allocation <- smoothed * trial[[denominator]]/sum(smoothed * trial[[denominator]])

    trial$expected_proportion <- expected_allocation/trial[[denominator]]
    trial$rowno <- seq(1:nrow(trial))

    # expand the vector of locations to allow for denominators > 1
    triallong <- trial %>%
        tidyr::uncount(trial[[denominator]])

    # sample generates a multinomial sample and outputs the indices of the locations assigned
    positives <- sample(x = nrow(triallong), size = npositives, replace = FALSE, prob = triallong$expected_proportion)
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

    # remove temporary variables and replace missing numerators with zero (because the multinomial sampling algorithm
    # leaves NA values where no events are assigned)
    trial <- subset(trial, select = -c(rowno, expected_proportion, sumnum))
    if (sum(is.na(trial[[numerator]])) > 0) {
        cat("** Warning: some records have zero denominator after rounding **\n")
        cat("You may want to remove these records or rescale the denominators \n")
        trial[[numerator]][is.na(trial[[numerator]])] <- 0
    }
    return(trial)
}

# deviation of ICC from target as a function of bandwidth
ICCdeviation <- function(logbw, trial, ICC_inp, approx_diag, sd, euclid, efficacy, initialPrevalence) {
    cluster <- NULL
    # set the seed so that a reproducible result is obtained for a specific bandwidth
    if (!is.null(logbw)) {
        bw <- exp(logbw)
        set.seed(round(bw * 1e+06))
    }

    trial <- syntheticBaseline(bw = bw, trial = trial, sd = sd, euclid = euclid, initialPrevalence = initialPrevalence)
    trial <- assignPositives(trial = trial, euclid = euclid, sd = sd, efficacy = efficacy, initialPrevalence = initialPrevalence,
        denominator = "denom")
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
syntheticBaseline <- function(bw, trial, sd, euclid, initialPrevalence) {
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

    trial <- distributePositives(trial = trial, initialPrevalence = initialPrevalence, smoothed = smoothedBaseline, efficacy = 0,
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


