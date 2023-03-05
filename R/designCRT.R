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

# Characteristics of a trial design. The input is a data frame. The output list
# conforms to the requirements for a CRT object
describeTrial <- function(trial, input.parameters = NULL) {

    # set the class to data.frame (removing the descriptors) so that
    # the function works both for data.frame and CRT input
    trial <- convertCRT.data.frame(trial)

    sd_distance <- mean_h <- sd_h <- min_k <- clustersRequired <- DE <- power <- NULL
    CRT.design <- list(locations = nrow(trial))

    if (!is.null(trial$cluster)) {

        # mean number of locations randomized in each cluster
        mean_h <- mean(table(trial$cluster))

        # standard deviation of locations randomized in each cluster
        sd_h <- stats::sd(table(trial$cluster))

    }
    if (!is.null(trial$arm)) {

        arm <- unique(cbind(trial$cluster, trial$arm))[, 2]  #assignments

        # Step L computation of characteristics of the randomization
        if (is.null(trial$nearestDiscord)) {
            # calculate nearest discord here
        }

        sd_distance <- stats::sd(trial$nearestDiscord)

        if (!is.null(input.parameters)) {

            calculateCRTpower <- calculateCRTpower(locations = CRT.design$locations, ICC = input.parameters$ICC,
                pC = input.parameters$pC, effect = input.parameters$effect,
                outcome.type = input.parameters$outcome.type, alpha = input.parameters$alpha,
                desiredPower = input.parameters$desiredPower, mean_h = mean_h,
                sd_h = sd_h)
            CRT.design <- calculateCRTpower$CRT.design.full
        }
    }
    CRT.design$mean_h <- mean_h
    CRT.design$sd_h <- sd_h
    CRT.design$sd_distance <- sd_distance
    return(CRT.design)
}


#' Summarise CRT
#'
#' \code{summary.CRT} generates a summary description of a CRT object
#' @param ... other arguments
#' @param object name of CRT
#' @export
summary.CRT <- function(object, ...) {
    cat("===============================CLUSTER RANDOMISED TRIAL ===========================\n")
    # create matrix
    output <- matrix("  ", nrow = 20, ncol = 2)
    rownames(output) <- c("Locations and Clusters\n-----------------------                               ",
        "Coordinate system            ", "Buffer width :               ",
        "Locations:                                            ", "Available clusters (across both arms)                 ",
        "  Per cluster mean number of locations                ", "  Per cluster s.d. number of locations                ",
        "S.D. of distance to nearest discordant location (km): ", "Cluster randomization:            ",
        "\nSpecification of Requirements\n-----------------------------",
        "Significance level:               ", "Type of Outcome                   ",
        "Expected outcome in control arm:  ", "Required effect size:             ",
        "\nPower calculations (ignoring contamination)\n------------------                    ",
        "Intra-cluster correlation: ", "Design effect:                         ",
        "Nominal power (%)                      ", paste0("Clusters required for power of ",
            object$input.parameters$desiredPower * 100, "%:     "), "Sufficient clusters for required power?")
    output[1, 1] <- "-"
    if (!is.null(object$x) & !is.null(object$y)) {
        output[2, 1] <- "(x, y)"
    } else if (!is.null(object$lat) & !is.null(object$long)) {
        output[2, 1] <- "Lat-Long"
    } else {
        output[2, 1] <- "No coordinates"
    }

    if (output[2, 1] == "(x, y)") {
        cat("\nSummary of coordinates\n----------------------\n")
        coordinate.summary <- with(object, summary(cbind(x, y)))
        rownames(coordinate.summary) <- substr(coordinate.summary[, 1], 1,
            8)
        coordinate.summary[, ] <- substr(coordinate.summary[, ], 9, 13)
        print(t(coordinate.summary))
        area <- pi * (as.numeric(coordinate.summary[6, 1]) - as.numeric(coordinate.summary[1,
            1])) * (as.numeric(coordinate.summary[6, 2]) - as.numeric(coordinate.summary[1,
            2]))
        cat("Approximate area (based on ellipse) : ", area, "sq.km\n\n")
    }
    if (is.null(object$cluster)) {
        output[5, 1] <- "Not assigned"
    } else {
        output[5, 1] <- "Assigned"
    }
    sd1 <- ifelse(is.null(object$CRT.design.full$sd_distance), NA, object$CRT.design.full$sd_distance)
    sd2 <- ifelse(is.null(object$CRT.design.core$sd_distance), NA, object$CRT.design.core$sd_distance)
    clustersAvailableFull <- with(object$CRT.design.full, floor(locations/mean_h))
    output[4, 1] <- object$CRT.design.full$locations
    output[5, 1] <- clustersAvailableFull
    output[6, 1] <- round(object$CRT.design.full$mean_h, digits = 1)
    output[7, 1] <- round(object$CRT.design.full$sd_h, digits = 1)
    if (!is.null(object$CRT.design.core)) {
      output[1, 1] <- "Full"
      output[1, 2] <- "Core"
      clustersAvailableCore <- with(object$CRT.design.core, floor(locations/mean_h))
        output[4, 2] <- object$CRT.design.core$locations
        output[5, 2] <- clustersAvailableCore
        output[6, 2] <- round(object$CRT.design.core$mean_h, digits = 1)
        output[7, 2] <- round(object$CRT.design.core$sd_h, digits = 1)
    }
    if (!is.null(object$arm)) {
        output[8, 1] <- ifelse(is.na(sd1), "", round(sd1, digits = 2))
        output[8, 2] <- ifelse(is.na(sd2), "", round(sd2, digits = 2))
        if (is.null(object$pair)) {
            output[9, 1] <- "Independently randomized"
        } else {
            output[9, 1] <- "Matched pairs randomized"
        }
    } else {
        # TODO: change wording for the case where there are no GIS data
        output[9, 1] <- "No randomization"
    }
    if (!is.null(object$input.parameters)) {
        output[10, 1] <- "-"
        output[11, 1] <- object$input.parameters$alpha
        output[12, 1] <- object$input.parameters$outcome.type
        output[13, 1] <- object$input.parameters$pC
        output[14, 1] <- object$input.parameters$effect
        output[16, 1] <- object$input.parameters$ICC
        if (!is.null(object$input.parameters$buffer.width)) {
            if (object$input.parameters$buffer.width > 0) {
                output[3, 1] <- paste0(object$input.parameters$buffer.width,
                  " km.")
            } else {
                output[3, 1] <- "No buffer"
            }
        }
    }
    output[15, 1] <- "-"
    if (is.null(object$input.parameters)) {
        rownames(output)[15] <- "No power calculations to report"
    } else {
        sufficient <- ifelse(clustersAvailableFull >= object$CRT.design.full$clustersRequired,
            "Yes", "No")
        output[17, 1] <- round(object$CRT.design.full$DE, digits = 1)
        output[18, 1] <- round(object$CRT.design.full$power * 100, digits = 1)
        output[19, 1] <- object$CRT.design.full$clustersRequired
        output[20, 1] <- sufficient

        if (is.null(object$CRT.design.core)) {
            output <- subset(output, select = -c(2))
        } else {
            output[15, 1] <- "Full"
            output[15, 2] <- "Core"
            clustersAvailableCore <- with(object$CRT.design.core, floor(locations/mean_h))
            sufficientCore <- ifelse(clustersAvailableCore >= object$CRT.design.core$clustersRequired,
                "Yes", "No")

            output[17, 2] <- round(object$CRT.design.core$DE, digits = 1)
            output[18, 2] <- round(object$CRT.design.core$power * 100, digits = 1)
            output[19, 2] <- object$CRT.design.core$clustersRequired
            output[20, 2] <- sufficientCore
        }
    }
    output <- output[(output[, 1] != "  "), ]
    # display and return table
    utils::write.table(output, quote = FALSE, col.names = FALSE, sep = "          ")
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
