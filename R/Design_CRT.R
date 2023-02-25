#' Design a CRT of a malaria intervention with contamination
#'
#' \code{Design_CRT} estimates the required number of clusters and the extent of contamination between arms for a CRT with the input set of locations.
#' Outputs are:
#' (i) Estimates of the required numbers of clusters.
#' (ii) A proposal for the cluster and arm assignments to the input coordinates. (A warning is output if the number of locations is too small to allow randomisation of sufficient clusters).
#' (iii) the proportion of households in the input geography falling within the core of the clusters (i.e. outside the contamination range of locations in the opposite arm)
#'
#' @param alpha confidence level
#' @param desiredPower desired power
#' @param effect required effect size
#' @param ICC Intra-Cluster Correlation obtained from other studies
#' @param pC baseline prevalence
#' @param buffer.width contamination range in km, obtained from other studies
#' @param coordinates dataframe containing coordinates of households. Columns 'x' and 'y' should contain Cartesian (x,y) coordinates. Units are expected to be km.
#' @param h  proposal for the number of coordinates in each cluster
#' @param algo algorithm for cluster boundaries, choose between
#' 'TSP': travelling salesman problem heuristic;
#' 'NN': nearest neighbor;
#' 'kmeans': kmeans algorithm
#' @param reuseTSP indicator of whether a pre-existing path should be used by the TSP algorithm
#' @return A list comprising a list the following elements:
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
#' \item \code{algo}: algorithm used for cluster boundaries. Options are
#' 'TSP': travelling salesman heuristic; 'NN': nearest neighbor; 'kmeans': kmeans algorithm
#' \item \code{n_ind}: algorithm used for cluster boundaries
#' \item \code{assignments}: data frame containing locations and cluster assignments
#' \item \code{min_c}: minimum number of clusters required
#' }
#' @export
#' @examples
#'
#' exampleDesign = Design_CRT(coordinates=readdata('test_site.csv'),
#'                 ICC=0.10, effect=0.4, pC=0.35, buffer.width=0.25, h=100)
Design_CRT <- function(alpha = 0.05, desiredPower = 0.8, effect, ICC, pC, buffer.width = 0, coordinates, h, algo = "kmeans",
    reuseTSP = FALSE) {

    # Step A: confidence level Step B: power Step C: Required effect size Step D: ICC, obtained from other studies
    # Step E: baseline prevalence Step F: buffer width based on postulated contamination range in km, obtained
    # from other studies Step G\t\tcoordinates of households in study area Step H: proposal for the number of
    # households in each cluster

    # convert power and significance level to normal deviates
    Zsig <- -qnorm(alpha/2)
    Zpow <- qnorm(desiredPower)

    # Power and sample size calculations based on Hemming et al, 2011
    # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-102


    # Step I: Calculations for the required minimum number of clusters for both arms

    pI <- pC * (1 - effect)  # probability in intervened group
    d <- pC - pI  # difference between groups
    sigma2 <- 1/2 * (pI * (1 - pI) + pC * (1 - pC))
    n_ind <- 2 * sigma2 * ((Zsig + Zpow)/d)^2  #required individuals per arm in individually randomized trial

    # see below for calculations of design effect and minimum numbers of clusters required

    # Step J: specify or compute cluster boundaries

    trial <- DefineClusters(trial = coordinates, h = h, algo = algo, reuseTSP = reuseTSP)

    # Step K: Random assignment of clusters to arms
    trial <- Randomize_CRT(trial)

    # augment the trial data frame with distance to nearest discordant coordinate (specifyBuffer assigns a buffer
    # only if a buffer width is > 0 is input)
    trial <- Specify_CRTbuffer(trial = trial, bufferWidth = buffer.width)

    # set the class to data.frame
    class(trial) <- "data.frame"

    input.parameters <- list(pC = pC, alpha = alpha, n_ind = n_ind, desiredPower = desiredPower, inputClusterSize = h,
        algo = algo, buffer.width = buffer.width, effect = effect, ICC = ICC, h = h)


    CRT.design.full <- describeTrial(trial = trial, pC = pC, d = d, desiredPower = desiredPower, n_ind = n_ind, sigma2 = sigma2,
        Zsig = Zsig, ICC = ICC)
    if (buffer.width > 0) {

        CRT.design.core <- describeTrial(trial = trial[trial$buffer == FALSE, ], pC = pC, d = d, desiredPower = desiredPower,
            n_ind = n_ind, sigma2 = sigma2, Zsig = Zsig, ICC = ICC)
    }
    class(trial) <- "CRT"
    trial$CRT.design.full <- CRT.design.full
    trial$CRT.design.core <- CRT.design.core
    trial$input.parameters <- input.parameters
    return(trial)
}

# Characteristics of a trial design
describeTrial <- function(trial, pC, d, desiredPower, n_ind, sigma2, Zsig, ICC) {

    arm <- unique(cbind(trial$cluster, trial$arm))[, 2]  #assignments

    k <- length(arm)/2  # number of clusters assigned to each arm

    ############################################################################## Step L computation of
    ############################################################################## characteristics of the
    ############################################################################## randomization

    sd_distance <- stats::sd(trial$nearestDiscord)

    # mean number of locations randomized to each arm
    mean_h <- mean(table(trial$cluster))

    # standard deviation of locations randomized to each arm
    sd_h <- stats::sd(table(trial$cluster))

    # coefficient of variation of the cluster sizes
    cv_h <- sd_h/mean_h

    # design effect (Variance Inflation Factor) allowing for varying cluster sizes(Hemming eqn 6)
    DE <- 1 + (cv_h^2 + 1) * (mean_h - 1) * ICC

    # number of individuals required per arm in CRT with equal cluster sizes
    n_crt <- n_ind * DE

    # minimum numbers of clusters required assuming varying cluster sizes per arm (Hemming eqn 8)
    min_k <- ceiling(n_ind * (1 + ((cv_h + 1) * mean_h - 1) * ICC)/mean_h)

    # power with k clusters per arm power = stats::pnorm(sqrt(k*mean_h/(2*(1+(mean_h-1)*ICC)))* d/sqrt(sigma2) -
    # Zsig) #Hemming eqn 27

    power <- stats::pnorm(sqrt(k * mean_h/(2 * DE)) * d/sqrt(sigma2) - Zsig)  #unequal cluster sizes


    CRT.design <- list(locations = length(trial$cluster), clusters = length(arm), nominalDE = DE, sd_distance = sd_distance,
        mean_h = mean_h, sd_h = sd_h, min_k = min_k, clustersRequired = 2 * min_k, clustersAssigned = length(arm), DE = DE,
        power = power)
    return(CRT.design)
}


#' Summarise CRT
#'
#' \code{summary.CRT} generates a summary description of a CRT object
#'
#' @param filename name of CRT
#' @export
summary.CRT <- function(trial) {
    cat("=====================CLUSTER RANDOMISED TRIAL DESIGN =================\n")
    cat("Site with ", length(trial$x), " total locations\n")
    if (is.null(trial$cluster)) {
        cat("No clusters assigned\n")
    } else {
        cat(length(table(trial$cluster)), "clusters assigned\n")
        if (is.null(trial$arm)) {
            cat("No randomization\n")
        } else if (is.null(trial$matchedPair)) {
            cat("Clusters independently randomized (no stratification)\n")
        } else {
            cat("Matched pairs of clusters randomized\n")
        }
    }
    if (!is.null(trial$input.parameters)) {
        cat("Significance level: ", trial$input.parameters$alpha, "\n")
        cat("Required effect size: ", trial$input.parameters$effect, "\n")
        cat("Assumed prevalence in absence of intervention ", trial$input.parameters$pC, "\n")
        cat("Pre-specified intra-cluster correlation: ", trial$input.parameters$ICC, "\n")
        cat("Buffer of width ", trial$input.parameters$buffer.width, " km.\n")
    }
    if (is.null(trial$CRT.design.full)) {
        cat("No power calculations included\n")
    } else {
        # create matrix
        output <- matrix(NA, nrow = 9, ncol = 2)

        # specify row and column names of matrix
        rownames(output) <- c("Locations", "  Per cluster mean number of locations", "  Per cluster s.d. number of locations",
            "Available clusters (across both arms)", "Clusters required", "Nominal power (%)", "Design effect", "S.D. of distance to nearest discordant location (km)",
            "Sufficient clusters for required power?")
        colnames(output) <- c("Full area", "Core area")
        sufficient <- ifelse(trial$CRT.design.full$clustersAssigned >= trial$CRT.design.full$clustersRequired, "Yes",
            "No")
        output[, 1] <- unlist(with(trial$CRT.design.full, c(locations, round(mean_h, digits = 1), round(sd_h, digits = 1),
            clustersAssigned, clustersRequired, round(power * 100, digits = 1), round(nominalDE, digits = 1), round(sd_distance,
                digits = 2), sufficient)))

        if (is.null(trial$CRT.design.core)) {
            output[, 2] <- NULL
        } else {
            sufficient <- ifelse(trial$CRT.design.core$clustersAssigned >= trial$CRT.design.core$clustersRequired, "Yes",
                "No")
            output[, 2] <- unlist(with(trial$CRT.design.core, c(locations, round(mean_h, digits = 1), round(sd_h, digits = 1),
                clustersAssigned, clustersRequired, round(power * 100, digits = 1), round(nominalDE, digits = 1), round(sd_distance,
                  digits = 2), sufficient)))
        }
        output <- format(output, digits = 1)
        # convert matrix to table
        output <- as.table(output)
        # display table
        output
    }
}
