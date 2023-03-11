#' Aggregate data across records with duplicated locations
#'
#' \code{aggregateCRT} aggregates data from a \code{CRT} object or data frame
#' containing multiple records with the same location, and outputs a list of S3 class \code{CRT}
#' containing single values for each location, for both the coordinates and the auxiliary variables.
#' @param trial a data frame containing locations (x,y) and other variables to be summed
#' @param auxiliaries vector of the names of auxiliary variables to be summed across each location
#' @returns A list with class \code{CRT} containing the following components:
#'  \tabular{lll}{
#'  \code{CRT.design.full}   \tab list \tab summary statistics describing the site\cr
#'  \code{x} \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y} \tab numeric vector \tab y-coordinates of locations \cr
#'  \code{...} \tab numeric vectors \tab auxiliary variables containing the sum(s) of the input auxiliaries\cr
#'  }
#' @export
aggregateCRT <- function(trial, auxiliaries = NULL) {
    trial <- convertCRT.data.frame(trial)
    x <- y <- NULL
    df <- trial %>%
        distinct(x, y, .keep_all = TRUE)
    df <- df[order(df$x, df$y), ]
    if (!is.null(auxiliaries)) {
        for (i in 1:length(auxiliaries)) {
            varName <- auxiliaries[i]
            df1 <- trial %>%
                group_by(x, y) %>%
                summarize(sumVar = sum(get(varName), na.rm = TRUE), .groups = "drop")
            df1 <- df1[order(df1$x, df1$y), ]
            df[[varName]] <- df1$sumVar
        }
    }
    CRT <- convert.data.frame.CRT(trial = df, input.parameters = NULL)
    return(CRT)
}

#' Specification of buffer zone in a cluster randomized trial
#'
#' \code{specify.buffer} specifies a buffer zone in a cluster randomized
#' trial (CRT) by specifying which locations are within a defined distance of
#' those in the opposite arm.
#'
#' @param trial an object of class \code{CRT} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}).
#' @param buffer.width minimum distance between locations in the core areas of
#'   opposing arms (km)
#' @returns an object of class \code{CRT} comprising the following components:
#'  \tabular{llll}{
#'  \code{CRT.design.full}   \tab list \tab summary statistics describing the site,
#'  cluster assignments, and randomization.\cr
#'  \code{CRT.design.core}   \tab list \tab summary statistics describing the core area \cr
#'  \code{x} \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y} \tab numeric vector \tab y-coordinates of locations \cr
#'  \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \code{arm} \tab factor \tab assignments to \code{control} or \code{intervention} for each location \cr
#'  \code{nearestDiscord} \tab numeric vector \tab Euclidean distance to nearest discordant location (km) \cr
#'  \code{buffer} \tab logical \tab indicator of whether the point is within the buffer \cr
#'  \code{...};   \tab other objects included in the input \code{CRT} object or data frame \cr
#'  }
#' @export
#' @examples
#' #Specify a buffer of 200m
#' exampletrial <- specify.buffer(trial = readdata('test_Arms.csv'), buffer.width = 0.2)
#'
specify.buffer <- function(trial, buffer.width = 0) {
  trial <- convertCRT.data.frame(trial)
  # nearestDiscord: nearest coordinate in the discordant arm, for the
  # control coordinates return the minimal distance with a minus sign
  if (is.null(trial$nearestDiscord)) trial$nearestDiscord <- get_nearestDiscord(trial)
  if (buffer.width > 0) {
    trial$buffer <- (abs(trial$nearestDiscord) < buffer.width)
  }
  CRT <- convert.data.frame.CRT(trial, input.parameters = NULL)
  return(CRT)
}

#' Randomize a two-armed cluster randomized trial
#'
#' \code{randomizeCRT} carries out randomization of clusters for a CRT and
#' augments the trial dataframe with assignments to arms \cr
#'
#' @param trial an object of class \code{CRT} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}). Optionally: specification of a buffer zone (logical \code{buffer});
#'   any other variables required for subsequent analysis.
#' @param matchedPair logical: indicator of whether pair-matching on the
#'   baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (required for
#'   matched-pair randomization)
#' @param baselineDenominator name of denominator variable for baseline data (required for
#'   matched-pair randomization)
#' @returns A list of class \code{CRT} containing the following components:
#'  \tabular{llll}{
#'  \code{CRT.design.full}   \tab list \tab summary statistics of the site \cr
#'  \code{CRT.design.core}   \tab list \tab summary statistics of the core area (if a buffer is present)\cr
#'  \code{x} \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y} \tab numeric vector \tab y-coordinates of locations \cr
#'  \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \code{pair} \tab factor \tab assignments to matched pair of each location  \cr
#'  (if \code{matchedPair} randomisation was carried out) \cr
#'  \code{arm} \tab factor \tab assignments to \code{control} or \code{intervention} for each location \cr
#'  \code{...};   \tab other objects included in the input \code{CRT} object or data frame  \cr
#'  }
#' @export
#' @examples
#' #Randomize the clusters in an example trial
#' set.seed(352)
#' example.CRT <- randomizeCRT(trial = readdata('test_Clusters.csv'), matchedPair = TRUE)

randomizeCRT <- function(trial, matchedPair = FALSE, baselineNumerator = "base_num",
    baselineDenominator = "base_denom") {

    trial <- convertCRT.data.frame(trial)

    # remove any preexisting assignments and coerce matchedPair to FALSE if there are no baseline data
    if(is.null(trial[[baselineNumerator]]) & matchedPair) {
        cat("** Warning: no baseline data for matching. Unmatched randomisation **\n")
        matchedPair <- FALSE
    }
    trial$arm <- NULL
    trial$pair <-NULL

    pair <- cluster <- base_num <- base_denom <- NULL
    trial$cluster <- as.factor(trial$cluster)
    # Randomization, assignment to arms
    nclusters <- length(unique(trial$cluster))
    if ((nclusters%%2) == 1 & matchedPair) {
        cat("** Warning: odd number of clusters: assignments are not matched on baseline data **\n")
        matchedPair <- FALSE
    }
    # uniformly distributed numbers, take mean and boolean of that
    rand_numbers <- runif(nclusters, 0, 1)
    if (matchedPair) {
        trial$base_num <- trial[[baselineNumerator]]
        trial$base_denom <- trial[[baselineDenominator]]
        cdf <- data.frame(trial %>%
            group_by(cluster) %>%
            dplyr::summarize(positives = sum(base_num), total = sum(base_denom)))
        cdf$p <- cdf$positives/cdf$total
        cdf <- cdf[order(cdf$p), ]
        cdf$pair <- rep(seq(1, nclusters/2), 2)
        cdf$rand_numbers <- rand_numbers
        cdf <- cdf[with(cdf, order(pair, rand_numbers)), ]
        cdf$arm <- rep(c(1, 0), nclusters/2)
        arm <- cdf$arm[order(cdf$cluster)]
        pair <- cdf$pair[order(cdf$cluster)]
    } else {
        arm <- ifelse(rand_numbers > median(rand_numbers), 1, 0)
    }
    if (matchedPair) trial$pair <- factor(pair[trial$cluster[]])
    trial$arm <- factor(arm[trial$cluster[]], levels = c(0, 1), labels = c("control",
        "intervention"))
    trial$nearestDiscord <- get_nearestDiscord(trial)
    CRT <- convert.data.frame.CRT(trial, input.parameters = NULL)
    return(CRT)
}



convert.data.frame.CRT <- function(trial, input.parameters){
  trial <- convertCRT.data.frame(trial)
  CRT.design.full <- describeTrial(trial, input.parameters = input.parameters)
  if (!is.null(trial$buffer)) {
    CRT.design.core <- describeTrial(trial = trial[trial$buffer == FALSE, ],
                    input.parameters = input.parameters)
  } else {
    CRT.design.core <- NULL
  }
  class(trial) <- "CRT"
  CRT <- trial
  CRT$CRT.design.full <-  CRT.design.full
  CRT$CRT.design.core <-  CRT.design.core
  CRT$input.parameters <- input.parameters
return(CRT)}



#' Assign each location to a cluster
#'
#' \code{specify.clusters} algorithmically assigns locations to clusters by grouping them geographically
#'
#' @param trial A CRT object or data frame containing (x,y) coordinates of
#'   households
#' @param h integer: proposal for the number of locations in each cluster
#' @param nclusters integer: number of clusters
#' @param algorithm algorithm for cluster boundaries, options are:
#' \code{NN} (the default),  \code{kmeans},  \code{TSP}
#' @param reuseTSP logical: indicator of whether a pre-existing path should be used by
#'   the TSP algorithm
#' @returns A list of class \code{CRT} containing the following components:
#'  \tabular{llll}{
#'  \code{CRT.design.full}   \tab list \tab summary statistics describing the site and cluster assignments.\cr
#'  \code{x} \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y} \tab numeric vector \tab y-coordinates of locations \cr
#'  \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \code{path} \tab numeric \tab travelling salesman path (created only if the TSP algorithm is used) \cr
#'  \code{...};   \tab other objects included in the input \code{CRT} object or data frame  \cr
#'  }
#' @details Clustering is carried out using one of three algorithms:
#' \tabular{llll}{
#' \code{NN}\tab Nearest neighbour:\tab assigns equal numbers of locations to each cluster \cr
#' \code{kmeans}\tab kmeans clustering  :\tab aims to partition locations so that each
#' belongs to the cluster with the nearest centroid \cr
#' \code{TSP} \tab travelling salesman problem heuristic :\tab Assigns locations sequentially
#' along a travelling salesman path.The \code{reuseTSP} parameter is used to allow the path to be reused
#' for creating alternative allocations with different cluster sizes. \cr
#' }
#' Either \code{nclusters} or \code{h} must be specified. If both are specified
#' the input value of \code{nclusters} is ignored.\cr
#' @export
#'
#' @examples
#' #Assign clusters to the test trial dataset averaging h = 40 using the kmeans algorithm
#' exampletrial <- specify.clusters(trial = readdata('test_site.csv'),
#'                             h = 40, algorithm = 'kmeans', reuseTSP = FALSE)
specify.clusters <- function(trial = trial, h = NULL, nclusters = NULL, algorithm = "NN",
    reuseTSP = FALSE) {
     trial <- convertCRT.data.frame(trial)


    # Local data from study area (ground survey and/or satellite
    # images)
    coordinates <- data.frame(x = as.numeric(as.character(trial$x)), y = as.numeric(as.character(trial$y)))

    # the number of clusters and the target cluster size must be integers.
    # cluster size can only be exactly equal to the input value of h if this is a factor of
    # the number of locations
    if (is.null(nclusters)) {
        nclusters <- ceiling(nrow(coordinates)/h)
    }
    if (is.null(h)) {
        h <- ceiling(nrow(coordinates)/nclusters)
    }
    # derive cluster boundaries

    if (algorithm == "TSP") {
        TSPoutput <- TSP_ClusterDefinition(coordinates, h, nclusters, reuseTSP)
        trial$path <- TSPoutput$path
        trial$cluster <- TSPoutput$cluster
    } else if (algorithm == "NN") {
        trial$cluster <- NN_ClusterDefinition(coordinates, h, nclusters)$cluster
    } else if (algorithm == "kmeans") {
        trial$cluster <- kmeans_ClusterDefinition(coordinates, nclusters)$cluster
    } else {
        stop("unknown method")
    }

    CRT <- convert.data.frame.CRT(trial, input.parameters = NULL)
    return(CRT)
}


#' Convert lat long co-ordinates to x,y
#'
#' \code{convert.latlong.xy} converts co-ordinates expressed as decimal degrees
#' into x,y
#' @param df data frame containing latitudes and longitudes in decimal degrees
#' @param latvar name of column containing latitudes in decimal degrees
#' @param longvar name of column containing longitudes in decimal degrees
#' @details An object containing the input locations replaced with Cartesian
#'   coordinates in units of km, centred on (0,0). Other data are unchanged.
#'   The equirectangular projection (valid for small areas) is used.
#' @returns A list of class \code{CRT} containing the following components:
#'  \tabular{llll}{
#'  \code{CRT.design.full}   \tab list \tab summary statistics describing the site \cr
#'  \code{x} \tab numeric vector \tab x-coordinates of locations \cr
#'  \code{y} \tab numeric vector \tab y-coordinates of locations \cr
#'  \code{...};   \tab other objects included in the input \code{CRT} object or data frame  \cr
#'  }
#' @export
convert.latlong.xy <- function(df, latvar = "lat", longvar = "long") {
    colnames(df)[identical(colnames(df), latvar)] <- "lat"
    colnames(df)[identical(colnames(df), longvar)] <- "long"
    R <- 6371  # radius of the earth
    latradians <- with(df, pi/180 * lat)
    longradians <- with(df, pi/180 * long)
    meanlat <- mean(latradians)
    meanlong <- mean(longradians)
    df$y <- R * (latradians - meanlat) * cos(longradians)
    df$x <- R * (longradians - meanlong)
    drops <- c("lat", "long")
    df <- df[, !(names(df) %in% drops)]
    CRT <- convert.data.frame.CRT(df, input.parameters = NULL)
    return(CRT)
}


#' Anonymize locations of a trial site
#'
#' \code{anonymize.site} carries out rotation of x,y coordinates a random angle
#' about a random origin.
#' @param trial CRT object or trial data frame with Cartesian co-ordinates of
#'   households (columns x and y)
#' @return object with re-centred co-ordinates the origin (other data are
#'   unchanged)
#' @format CRT object:
#' \itemize{
#' \item \code{CRT.design.full} list of characteristics of the full area
#' \item \code{x} x-coordinates of location
#' \item \code{y} y-coordinates of location
#' }
#' @export
#' @examples
#' #Rotate and reflect test site locations
#' transformedTestlocations <- anonymize.site(trial = readdata('test_site.csv'))

anonymize.site <- function(trial) {
    # Local data from study area (ground survey and/or satellite
    # images) random rotation angle
    trial <- convertCRT.data.frame(trial)
    theta <- 2 * pi * runif(n = 1)
    x <- trial$x
    y <- trial$y
    rangex <- max(x) - min(x)
    rangey <- max(y) - min(y)
    translation <- c(rangex * rnorm(n = 1), rangey * rnorm(n = 1))

    xy <- t(matrix(c(x, y), ncol = 2, nrow = length(x)))
    xytranslated <- xy + translation

    rotation <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
        nrow = 2, ncol = 2)

    # Rotate
    xytrans <- rotation %*% xytranslated

    # Recentre on origin
    recentred <- xytrans - c(mean(xytrans[1, ]), mean(xytrans[2, ]))
    trial$x <- recentred[1, ]
    trial$y <- recentred[2, ]

    CRT <- convert.data.frame.CRT(trial, input.parameters = NULL)
    return(CRT)
}


#' Read test dataset
#'
#' \code{readdata} reads a file from the package library of test datasets
#'
#' @param filename name of text file stored within the package
#' @return R object corresponding to the text file
#' @details The input file name should include the extension (either .csv or .txt).
#' The resulting object is a data frame if the extension is .csv.
#' @export
readdata <- function(filename) {
    fname <- eval(filename)
    extdata <- system.file("extdata", package = "CRTspillover")
    if (unlist(gregexpr(".txt", fname)) > 0)
        robject <- dget(file = paste0(extdata, "/", fname))
    if (unlist(gregexpr(".csv", fname)) > 0)
        robject <- read.csv(file = paste0(extdata, "/", fname), row.names = NULL)
    # remove variable 'X' if it is present
    robject$X <- NULL
    return(robject)
}

#' Convert object of S3 class CRT to a data frame
#'
#' \code{convertCRT.data.frame} removes lists of descriptors and summary statistics
#'
#' @param CRT name of CRT object
#' @return data frame with one row for each location
#' @details \code{convertCRT.data.frame} removes lists of descriptors and summary statistics
#' @export
convertCRT.data.frame <- function(CRT) {
    CRT$CRT.design.full <- NULL
    CRT$CRT.design.core <- NULL
    CRT$input.parameters <- NULL
    trial <- CRT
    class(trial) <- "data.frame"
    return(trial)
}

# Characteristics of a trial design. The input is a data frame. The output list
# conforms to the requirements for a CRT object
describeTrial <- function(trial, input.parameters = NULL) {

  # set the class to data.frame (removing the descriptors) so that
  # the function works both for data.frame and CRT input
  trial <- convertCRT.data.frame(trial)

  sd_distance <- mean_h <- sd_h <- min_k <- clustersRequired <- DE <- power <- NULL

  coordinates <- data.frame(cbind(x=trial$x,y=trial$y))
  CRT.design <- list(records = nrow(trial),
                     locations = nrow(dplyr::distinct(coordinates)))

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
#' @param object name of CRT
#' @param ... other arguments
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
    rownames(coordinate.summary) <- substr(coordinate.summary[, 1], 1, 8)
    coordinate.summary[, ] <- substr(coordinate.summary[, ], 9, 13)
    print(t(coordinate.summary))
    area <- pi * (as.numeric(coordinate.summary[6, 1]) -
                  as.numeric(coordinate.summary[1, 1])) *
      (as.numeric(coordinate.summary[6, 2]) - as.numeric(coordinate.summary[1, 2]))
    cat("Approximate area (based on ellipse) : ", area, "sq.km\n\n")
    # TODO replace with area calculation using sf::st_area(x, ...)
  }
  if (is.null(object$cluster)) {
    output[5, 1] <- "Not assigned"
  } else {
    clustersAvailableFull <- with(object$CRT.design.full, floor(locations/mean_h))
    output[5, 1] <- clustersAvailableFull
    output[6, 1] <- round(object$CRT.design.full$mean_h, digits = 1)
    output[7, 1] <- round(object$CRT.design.full$sd_h, digits = 1)
  }
  if(identical(object$CRT.design.full$locations,object$CRT.design.full$records)){
    output[4, 1] <- object$CRT.design.full$locations
  } else {
    rownames(output)[4] <- paste0("Aggregation required. Total records: ",object$CRT.design.full$records,
                                  ". Total locations:")
    output[4, 1] <- object$CRT.design.full$locations
  }
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
    sd1 <- ifelse(is.null(object$CRT.design.full$sd_distance), NA, object$CRT.design.full$sd_distance)
    sd2 <- ifelse(is.null(object$CRT.design.core$sd_distance), NA, object$CRT.design.core$sd_distance)
    output[8, 1] <- ifelse(is.na(sd1), "", round(sd1, digits = 2))
    output[8, 2] <- ifelse(is.na(sd2), "", round(sd2, digits = 2))
    if (is.null(object$pair)) {
      output[9, 1] <- "Independently randomized"
    } else {
      output[9, 1] <- "Matched pairs randomized"
    }
  } else {
    if (is.null(object$x)) {
      output[9, 1] <- "No locations to randomize"
    } else {
      output[9, 1] <- "No randomization"
    }
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



# compute vector of distances to nearest discordant location
get_nearestDiscord <- function(trial){
  dist_trial <- as.matrix(dist(cbind(trial$x, trial$y), method = "euclidean"))
  discord <- outer(trial$arm, trial$arm, "!=")  #true & false.
  discord_dist_trial <- ifelse(discord, dist_trial, Inf)
  nearestDiscord <- ifelse(trial$arm == "control", -apply(discord_dist_trial,
                              MARGIN = 2, min), apply(discord_dist_trial, MARGIN = 2, min))
  return(nearestDiscord)
}

TSP_ClusterDefinition <- function(coordinates, h, nclusters, reuseTSP) {

  if (!"path" %in% colnames(coordinates) | !reuseTSP) {
    # Code originally from Silkey and Smith, SolarMal

    # Order the coordinates along an optimised travelling
    # salesman path
    dist_coordinates <- dist(coordinates, method = "euclidean")
    tsp_coordinates <- TSP::TSP(dist_coordinates)  # object of class TSP
    tsp_coordinates <- TSP::insert_dummy(tsp_coordinates)
    tour <- TSP::solve_TSP(tsp_coordinates, "repetitive_nn")  #solves TSP, expensive
    path <- TSP::cut_tour(x = tour, cut = "dummy")
    coordinates$path <- path

  }
  # order coordinates
  coordinates$order <- seq(1:nrow(coordinates))
  coordinates <- coordinates[order(coordinates$path), ]

  n1 <- (nclusters - 1) * h
  nclusters_1 <- nclusters - 1
  # The last cluster may be a different size (if h is not a
  # factor of the population size) )
  coordinates$cluster <- NA
  coordinates$cluster[1:n1] <- c(rep(1:nclusters_1, each = h))  #add cluster assignment
  coordinates$cluster[which(is.na(coordinates$cluster))] <- nclusters
  coordinates <- coordinates[order(coordinates$order), ]
  return(coordinates)
}

NN_ClusterDefinition <- function(coordinates, h, nclusters) {

  # algorithm is inspired by this website: ??? (comment from Lea)

  # initialize cluster, calculate euclidean distance
  dist_coordinates <- as.matrix(dist(coordinates, method = "euclidean"))
  coordinates$cluster <- NA

  nclusters_1 <- nclusters - 1
  for (i in 1:nclusters_1) {

    # find unassigned coordinates
    cluster_unassigned <- which(is.na(coordinates$cluster))
    dist_coordinates_unassigned <- dist_coordinates[cluster_unassigned,
                                                    cluster_unassigned]
    cluster_na <- rep(NA, length(cluster_unassigned))

    # find the coordinate furthest away from all the others
    index <- which.max(rowSums(dist_coordinates_unassigned))

    # find the n nearest neighbors of index
    cluster_na[head(order(dist_coordinates_unassigned[index, ]),
                    h)] <- i
    coordinates$cluster[cluster_unassigned] <- cluster_na
  }
  # The last cluster may be a different size (if h is not a
  # factor of the population size) )
  coordinates$cluster[which(is.na(coordinates$cluster))] <- nclusters

  return(coordinates)
}

kmeans_ClusterDefinition <- function(coordinates, nclusters) {

  # kmeans as implemented in R base
  km <- kmeans(x = coordinates, centers = nclusters)
  coordinates$cluster <- km$cluster

  return(coordinates)
}
