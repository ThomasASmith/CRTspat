#' Simulation of cluster randomized trial with contamination
#'
#' \code{Simulate_CRT} generates simulated data for a cluster randomized trial (CRT) with geographic contamination between arms and a binary outcome. Contamination is simulated as arising from a diffusion-like process.
#' For details see Multerer (PhD thesis).
#'
#' @param trial a dataframe containing locations (x,y), cluster assignments, and arm assignments
#' @param efficacy simulated efficacy (defaults to 0)
#' @param initialPrevalence prevalence in control arm (assumed equal to initial proportion)
#' @param generateBaseline logical indicator of whether baseline data should be simulated
#' @param matchedPair logical: indicator of whether pair-matching on the baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (if present)
#' @param baselineDenominator name of denominator variable for baseline data (if present)
#' @param ICC_inp Intra Cluster Correlation, provided as input when baseline data are to be simulated
#' @param sd standard deviation of the normal kernel measuring spatial smoothing leading to contamination
#' @param theta_inp input contamination range
#' @param tol tolerance of output ICC
#' @return A trial simulation object comprising a data frame containing the following numerical quantities:
#' \itemize{
#' \item \code{x}: x-coordinates of location
#' \item \code{y}: y-coordinates of location
#' \item \code{cluster}: assignment to cluster
#' \item \code{arm}: assignment to trial arm
#' \item \code{infectiousness_proxy}: infectiousness proxy
#' \item \code{base_denom}: denominator for baseline
#' \item \code{base_num}: numerator for baseline
#' \item \code{nearestDiscord}: distance to nearest discordant location
#' \item \code{denom}: denominator for the outcome
#' \item \code{num}: numerator for the outcome
#' }
#' @export
#'
#' @examples
#' example_simulated_CRT =  Simulate_CRT(trial=CRTspillover::testArms,,
#'                                      efficacy=0.25,
#'                                      ICC_inp=0.05,
#'                                      initialPrevalence=0.5,
#'                                      matchedPair = TRUE,
#'                                      sd=0.6,
#'                                      tol=0.05)
Simulate_CRT <- function(trial = NULL, efficacy = 0, initialPrevalence = NULL,
                         generateBaseline = TRUE, matchedPair = TRUE, baselineNumerator = "base_num",
                         baselineDenominator = "base_denom", ICC_inp = NULL, sd = NULL, theta_inp = NULL,
                         tol = 1e-04) {

# Written by Tom Smith, July 2017. Adapted by Lea Multerer, September 2017
  cat("\n=====================    SIMULATION OF CLUSTER RANDOMISED TRIAL    =================\n")
  bw <- NULL
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

  # For the smoothing step compute contributions to the relative
  # effect size from other households as a function of distance to
  # the other households

  euclid <- distance_matrix(trial$x, trial$y)

  # generate baseline data if required and exposure proxy if this
  # is not provided

  if (!"infectiousness_proxy" %in% colnames(trial) & baselineNumerator %in%
      colnames(trial) & baselineDenominator %in% colnames(trial)) {
    trial$infectiousness_proxy <- trial[[baselineNumerator]]/trial[[baselineDenominator]]
  } else if ("infectiousness_proxy" %in% colnames(trial)) {
    # create a baseline dataset using a pre-existing exposure
    # proxy
    trial <- syntheticBaseline(bw = NULL, trial = trial, sd = sd, euclid = euclid,
                               initialPrevalence = initialPrevalence)
  } else if (generateBaseline) {
    # determine the required smoothing bandwidth by fitting to
    # the pre-specified ICC

    # compute approximate diagonal of clusters
    approx_diag <- sqrt((max(trial$x) - min(trial$x))^2 + (max(trial$y) -
                         min(trial$y))^2)/sqrt(length(unique(trial$cluster)))

    cat("Estimating the smoothing required to achieve the target ICC of",
        ICC_inp, "\n")

    loss <- 999
    while (loss > tol) {
      ICC.loss <- OOR::StoSOO(par = NA, fn = ICCdeviation, lower = -5,
                              upper = 5, nb_iter = 20, trial = trial, ICC_inp = ICC_inp,
                              approx_diag = approx_diag, sd = sd, euclid = euclid, efficacy = efficacy,
                              initialPrevalence = initialPrevalence)
      loss <- ICC.loss$value
    }
    logbw <- ICC.loss$par
    # overprint the output that was recording progress
    cat("\r                                                         \n")

    # set the seed so that the same result is obtained for a
    # specific bandwidth
    bw <- exp(logbw)
    set.seed(round(bw * 1e+06))
    # create a baseline dataset using the optimized bandwidth
    trial <- syntheticBaseline(bw = bw, trial = trial, sd = sd, euclid = euclid,
                               initialPrevalence = initialPrevalence)
  }

  trial <- assignPositives(trial = trial, euclid = euclid, sd = sd, efficacy = efficacy,
                           initialPrevalence = initialPrevalence)

  return(trial)
}

# Assign expected proportions to each location assuming a fixed efficacy.
assignPositives <- function(trial, euclid, sd, efficacy, initialPrevalence) {

  # Indicator of whether the source is intervened is
  # (as.numeric(trial$arm[i]) - 1 smoothedIntervened is the value
  # of infectiousness_proxy decremented by the effect of
  # intervention and smoothed to allow for mosquito movement

  if (sd > 0) {
    smoothedIntervened <- gauss(sd, euclid) %*% (trial$infectiousness_proxy *
                                                   (1 - efficacy * (as.numeric(trial$arm) - 1)))
  } else {
    smoothedIntervened <- trial$infectiousness_proxy * (1 - efficacy *
                                                          (as.numeric(trial$arm) - 1))
  }
  # distances to nearest discordant households
  discord <- outer(trial$arm, trial$arm, "!=")  #returns true & false.
  euclidd <- ifelse(discord, euclid, 99999.9)

  # for the control arm return the minimal distance with a minus sign
  trial$nearestDiscord <- ifelse(trial$arm == "control", -apply(euclidd,MARGIN = 2, min),
                                                          apply(euclidd, MARGIN = 2, min))

  if (!("denom" %in% colnames(trial))) trial$denom <- 1

  # the denominator must be an integer; this changes the denominator if a non-integral value is input
  trial$denom <- round(trial$denom, digits = 0)

  # compute the total positives expected given the input efficacy
  npositives <- round(initialPrevalence * sum(trial$denom) * (1 - 0.5 * efficacy))

  trial <- distributePositives(trial = trial, npositives = npositives, smoothed = smoothedIntervened)
  return(trial)
}

# allocate the positives to locations
distributePositives <- function(trial, npositives, smoothed){
  denom <- expected_proportion <- num <- rowno <- NULL

  # scale to input value of initial prevalence by assigning required number of infections with
  # probabilities proportional to smoothedIntervened multiplied by the denominator

  expected_allocation <- smoothed * trial$denom/sum(smoothed * trial$denom)

  trial$expected_proportion <- expected_allocation/trial$denom
  trial$rowno <- seq(1:length(trial$denom))

  # expand the vector of locations to allow for denominators > 1
  triallong <- trial %>% tidyr::uncount(denom)

  # sample generates a multinomial sample and outputs the indices
  # of the locations assigned
  positives <- sample(x = nrow(triallong), size = npositives, replace = FALSE,
                      prob = triallong$expected_proportion)
  triallong$num <- 0
  triallong$num[positives] <- 1

  # summarise the numerator values into the original set of locations
  num <- dplyr::group_by(triallong, rowno) %>%
    dplyr::summarise(num = sum(num))

  # use left_join to merge into the original data frame (records
  # with zero denominator do not appear in num)
  trial <- trial %>% dplyr::left_join(num, by="rowno")

  # remove temporary variables and replace numerators with zero
  # where denominator is zero
  trial <- subset(trial, select = -c(rowno, expected_proportion))
  if (sum(is.na(trial$num)) > 0) {
    cat("** Warning: some records have zero denominator after rounding **\n")
    cat("You may want to remove these records or rescale the denominators \n")
    trial$num[is.na(trial$num)] <- 0
  }
  return(trial)
}

# deviation of ICC from target as a function of bandwidth
ICCdeviation <- function(logbw, trial, ICC_inp, approx_diag, sd, euclid,
                         efficacy, initialPrevalence) {
  cluster <- NULL
  # set the seed so that a reproducible result is obtained for a specific bandwidth
  if (!is.null(logbw)) {
    bw <- exp(logbw)
    set.seed(round(bw * 1e+06))
  }

  trial <- syntheticBaseline(bw = bw, trial = trial, sd = sd, euclid = euclid,
                             initialPrevalence = initialPrevalence)
  trial <- assignPositives(trial = trial, euclid = euclid, sd = sd, efficacy = efficacy,
                           initialPrevalence = initialPrevalence)
  trial$neg <- trial$denom - trial$num
  fit <- geepack::geeglm(cbind(num, neg) ~ arm, id = cluster, corstr = "exchangeable",
                         data = trial, family = binomial(link = "logit"))
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
    trial$infectiousness_proxy <- KDESmoother(trial$x, trial$y, kernnumber = 200,
                                              bandwidth = bw, low = 0, high = 1)
  }
  # Smooth the exposure proxy to allow for mosquito movement Note
  # that the s.d. in each dimension of the 2 d gaussian is
  # sd/sqrt(2) smoothedBaseline is the amount received by the each
  # cluster from the contributions (infectiousness_proxy) of each
  # source
  if (sd > 0) {
    smoothedBaseline <- gauss(sd, euclid) %*% trial$infectiousness_proxy
  } else {
    smoothedBaseline <- trial$infectiousness_proxy
  }

  # scale to input value of initial prevalence by assigning
  # required number of infections with probabilities proportionate
  # to infectiousness_proxy

  # number of positives required to match the specified prevalence
  npositives <- round(initialPrevalence * nrow(trial))

  positives <- sample(x = nrow(trial), size = npositives, replace = FALSE,
                      prob = smoothedBaseline)

  trial$base_denom <- 1
  trial$base_num <- 0
  trial$base_num[positives] <- 1
  return(trial)
}

# compute a euclidian distance matrix
distance_matrix <- function(x, y) {
  #generates square matrices of differences
  xdist <- outer(x, x, "-")
  ydist <- outer(y, y, "-")
  euclid <- sqrt(xdist * xdist + ydist * ydist)
}


# add lognormal noise:
# not sure this function is needed
# X is the input vector comprising a sample from a smoothed distribution
# varXY is the required variance
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
  # Careful here, totalf sums over very small numbers, consider
  # replacing
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


