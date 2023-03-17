#' Aggregate data across records with duplicated locations
#'
#' \code{aggregateCRT} aggregates data from a \code{"CRT"} object or data frame
#' containing multiple records with the same location, and outputs a list of S3 class \code{"CRT"}.
#' containing single values for each location, for both the coordinates and the auxiliary variables.
#' @param trial An object of S3 class \code{"CRT"} or a data frame containing locations (x,y) and variables to be summed
#' @param auxiliaries vector of names of auxiliary variables to be summed across each location
#' @returns A list of S3 class \code{"CRT"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site\tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  \tab \code{...} \tab numeric vectors: \tab auxiliary variables containing the sum(s) of the input auxiliaries\cr
#'  }
#' @export
aggregateCRT <- function(trial, auxiliaries = NULL) {
    if (identical(class(trial),"data.frame")){
      CRT <- list(trial = trial, design = NULL)
      class(CRT) <- "CRT"
    } else {
      CRT <- trial
      trial <- CRT$trial
    }
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
    CRT <- updateCRT(CRT = CRT, trial = trial)
    return(CRT)
}

#' Specification of buffer zone in a cluster randomized trial
#'
#' \code{specify.buffer} specifies a buffer zone in a cluster randomized
#' trial (CRT) by flagging those locations that are within a defined distance of
#' those in the opposite arm.
#'
#' @param trial an object of class \code{"CRT"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}).
#' @param buffer.width minimum distance between locations in
#'   opposing arms for them to qualify to be included in the core area (km)
#' @returns A list of S3 class \code{"CRT"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site,
#'  cluster assignments, and randomization.\tab\cr
#'  \code{geom.core}   \tab list: \tab summary statistics describing the core area \tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \tab \code{arm} \tab factor: \tab assignments to \code{"control"} or \code{"intervention"} for each location \cr
#'  \tab \code{nearestDiscord} \tab numeric vector: \tab Euclidean distance to nearest discordant location (km) \cr
#'  \tab \code{buffer} \tab logical: \tab indicator of whether the point is within the buffer \cr
#'  \tab \code{...} \tab \tab other objects included in the input \code{"CRT"} object or data frame \cr
#'  }
#' @export
#' @examples
#' #Specify a buffer of 200m
#' exampletrial <- specify.buffer(trial = readdata('test_Arms.csv'), buffer.width = 0.2)
specify.buffer <- function(trial, buffer.width = 0) {
  if (identical(class(trial),"data.frame")){
    CRT <- list(trial = trial, design = NULL)
    class(CRT) <- "CRT"
  } else {
    CRT <- trial
    trial <- CRT$trial
  }
  # nearestDiscord: nearest coordinate in the discordant arm, for the
  # control coordinates return the minimal distance with a minus sign
  if (is.null(trial$nearestDiscord)) trial$nearestDiscord <- get_nearestDiscord(trial)
  if (buffer.width > 0) {
    trial$buffer <- (abs(trial$nearestDiscord) < buffer.width)
  }
  CRT <- updateCRT(CRT = CRT, trial = trial)
  return(CRT)
}

#' Randomize a two-armed cluster randomized trial
#'
#' \code{randomizeCRT} carries out randomization of clusters for a CRT and
#' augments the trial dataframe with assignments to arms \cr
#'
#' @param trial an object of class \code{"CRT"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}). Optionally: specification of a buffer zone (logical \code{buffer});
#'   any other variables required for subsequent analysis.
#' @param matchedPair logical: indicator of whether pair-matching on the
#'   baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (required for
#'   matched-pair randomization)
#' @param baselineDenominator name of denominator variable for baseline data (required for
#'   matched-pair randomization)
#' @returns A list of S3 class \code{"CRT"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site,
#'  cluster assignments, and randomization.\tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \tab \code{pair} \tab factor \tab assigned matched pair of each location
#'  (if \code{matchedPair} randomisation was carried out) \cr
#'  \tab \code{arm} \tab factor: \tab assignments to \code{"control"} or \code{"intervention"} for each location \cr
#'  \tab \code{...} \tab \tab other objects included in the input \code{"CRT"} object or data frame \cr
#'  }
#' @export
#' @examples
#' #Randomize the clusters in an example trial
#' example.CRT <- randomizeCRT(trial = readdata('test_Clusters.csv'), matchedPair = TRUE)
randomizeCRT <- function(trial, matchedPair = FALSE, baselineNumerator = "base_num",
    baselineDenominator = "base_denom") {

    if (identical(class(trial),"data.frame")){
      CRT <- list(trial = trial, design = NULL)
      class(CRT) <- "CRT"
    } else {
      CRT <- trial
      trial <- CRT$trial
    }

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

    CRT <- updateCRT(CRT = CRT, trial = trial)
    return(CRT)
}


# Update CRT object with trial data frame and add geometry descriptions
updateCRT <- function(CRT, trial){
  CRT$trial <- trial
  CRT$geom.full <- get_geom(trial, design = CRT$design)
  if (!is.null(trial$buffer)) {
    CRT$geom.core <- get_geom(trial = trial[trial$buffer == FALSE, ],
                               design = CRT$design)
  } else {
    CRT$geom.core <- NULL
  }
  return(CRT)
}



as_CRT <- function(trial, design = NULL){
  geom.full <- get_geom(trial, design = design)
  if (!is.null(trial$buffer)) {
    geom.core <- get_geom(trial = trial[trial$buffer == FALSE, ],
                    design = design)
  } else {
    geom.core <- NULL
  }
  CRT <- list(trial = trial,
              geom.full = geom.full,
              geom.core = geom.core,
              design = design)
  class(CRT) <- "CRT"
return(CRT)}



#' algorithmically assigns locations to clusters in a CRT
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
#' @returns A list of S3 class \code{"CRT"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site,
#'  and cluster assignments.\tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor \tab assignments to cluster of each location  \cr
#'  \tab \code{...} \tab \tab other objects included in the input \code{"CRT"} object or data frame \cr
#'  }
#' @details Clustering is carried out using one of three algorithms:
#' \tabular{lll}{
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

   if (identical(class(trial),"data.frame")){
      CRT <- list(trial = trial, design = NULL)
      class(CRT) <- "CRT"
    } else {
      CRT <- trial
      trial <- CRT$trial
    }

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

    CRT <- updateCRT(CRT = CRT, trial = trial)
    return(CRT)
}


#' Convert lat long co-ordinates to x,y
#'
#' \code{latlong_as_xy} converts co-ordinates expressed as decimal degrees
#' into x,y
#' @param df data frame containing latitudes and longitudes in decimal degrees
#' @param latvar name of column containing latitudes in decimal degrees
#' @param longvar name of column containing longitudes in decimal degrees
#' @details An object containing the input locations replaced with Cartesian
#'   coordinates in units of km, centred on (0,0). Other data are unchanged.
#'   The equirectangular projection (valid for small areas) is used.
#' @returns A list of S3 class \code{"CRT"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site \tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  \tab \code{...} \tab \tab other objects included in the input \code{"CRT"} object or data frame \cr
#'  }
#' @export
latlong_as_xy <- function(df, latvar = "lat", longvar = "long") {
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
    CRT <- as_CRT(df, design = NULL)
    return(CRT)
}


#' Anonymize locations of a trial site
#'
#' \code{anonymize.site} carries out
#' @param trial \code{"CRT"} object or trial data frame with Cartesian co-ordinates of
#'   households (columns x and y)
#' @returns A list of S3 class \code{"CRT"} containing the following components:
#'  \tabular{llll}{
#'  \code{geom.full}   \tab list: \tab summary statistics describing the site \tab\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\tab\cr
#'  \tab \code{x} \tab numeric vector: \tab x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: \tab y-coordinates of locations \cr
#'  \tab \code{...} \tab \tab other objects included in the input \code{"CRT"} object or data frame \cr
#'  }
#' @export
#' @details The x,y coordinates are rotated by a random angle about a random origin. The returned object
#' has transformed co-ordinates re-centred at the origin. Other data are unchanged.
#' @examples
#' #Rotate and reflect test site locations
#' transformedTestlocations <- anonymize.site(trial =  readdata("test_Simulate_CRT.csv"))

anonymize.site <- function(trial) {
    # Local data from study area (ground survey and/or satellite
    # images) random rotation angle

    if (identical(class(trial),"data.frame")){
      CRT <- list(trial = trial, design = NULL)
      class(CRT) <- "CRT"
    } else {
      CRT <- trial
      trial <- CRT$trial
    }

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

    CRT <- updateCRT(CRT = CRT, trial = trial)
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
#' \code{CRT_as_data.frame} removes lists of descriptors and summary statistics
#' and returns the trial data frame
#'
#' @param CRT name of CRT object
#' @return data frame with one row for each location
#' @export
CRT_as_data.frame <- function(CRT) {
    CRT$geom.full <- NULL
    CRT$geom.core <- NULL
    CRT$design <- NULL
    trial <- CRT$trial
    class(trial) <- "data.frame"
    return(trial)
}

# Characteristics of a trial design. The input is a data frame. The output list
# conforms to the requirements for a CRT object
get_geom <- function(trial, design = NULL) {

  if (identical(class(trial),"data.frame")){
    CRT <- list(trial = trial, design = NULL)
    class(CRT) <- "CRT"
  } else {
    CRT <- trial
    trial <- CRT$trial
  }

  sd_distance <- mean_h <- sd_h <- min_k <- clustersRequired <- DE <- power <- NULL

  coordinates <- data.frame(cbind(x=trial$x,y=trial$y))
  geom <- list(records = nrow(trial),
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

    if (!is.null(design)) {

      calculateCRTpower <- calculateCRTpower(locations = geom$locations,
                                             alpha = design$alpha,
                                             desiredPower = design$desiredPower,
                                             effect = design$effect,
                                             yC = design$yC,
                                             outcome.type = design$outcome.type,
                                             ICC = design$ICC,
                                             sigma2 = design$sigma2,
                                             phi = design$phi,
                                             N = design$N,
                                             mean_h = mean_h,
                                             sd_h = sd_h)

      geom <- calculateCRTpower$geom.full
    }
  }
  geom$mean_h <- mean_h
  geom$sd_h <- sd_h
  geom$sd_distance <- sd_distance
  return(geom)
}



#' Summarise CRT
#'
#' \code{summary.CRT} generates a summary description of a CRT object
#' @param object name of CRT
#' @param maskbuffer buffer in km around locations to use for calculating total area
#' @param ... other arguments
#' @export
summary.CRT <- function(object, maskbuffer = 0.2, ...) {
  cat("===============================CLUSTER RANDOMISED TRIAL ===========================\n")
  output <- matrix("  ", nrow = 22, ncol = 2)
  rownames(output) <- paste0("row ", 1:nrow(output))
  rownames(output)[1] <- "Locations and Clusters\n----------------------                                "
  output[1, 1] <- "-"
  rownames(output)[2] <- "Coordinate system            "
  if (!is.null(object$trial$x) & !is.null(object$trial$y)) {
    output[2, 1] <- "(x, y)"
  } else if (!is.null(object$trial$lat) & !is.null(object$trial$long)) {
    output[2, 1] <- "Lat-Long"
  } else {
    output[2, 1] <- "No coordinates in dataset"
  }

  if (output[2, 1] == "(x, y)") {
    cat("\nSummary of coordinates\n----------------------\n")
    coordinate.summary <- with(object$trial, summary(cbind(x, y)))
    rownames(coordinate.summary) <- substr(coordinate.summary[, 1], 1, 8)
    coordinate.summary[, ] <- substr(coordinate.summary[, ], 9, 13)
    print(t(coordinate.summary))
    xycoords <- data.frame(cbind(x=object$trial$x,y=object$trial$y))
    tr <- sf::st_as_sf(xycoords, coords = c("x","y"))
    buf1 <- sf::st_buffer(tr, maskbuffer)
    buf2 <- sf::st_union(buf1)
    area <- sf::st_area(buf2)
    cat("Total area (within ", maskbuffer,"km of a location) : ", format(area, digits = 3), "sq.km\n\n")
  }
  rownames(output)[5] <- "Available clusters (across both arms)                 "
  if (is.null(object$geom.full$mean_h)) {
    output[5, 1] <- "Not assigned"
  } else {
    clustersAvailableFull <- with(object$geom.full, floor(locations/mean_h))
    output[5, 1] <- clustersAvailableFull
    rownames(output)[6] <- "  Per cluster mean number of points                   "
    output[6, 1] <- round(object$geom.full$mean_h, digits = 1)
    rownames(output)[7] <- "  Per cluster s.d. number of points                   "
    output[7, 1] <- round(object$geom.full$sd_h, digits = 1)
  }
  rownames(output)[4] <- "Locations:                                            "
  if(identical(object$geom.full$locations,object$geom.full$records)){
    output[4, 1] <- object$geom.full$locations
  } else {
    if(!is.null(object$geom.full$records)) {
      rownames(output)[4] <- paste0("Not aggregated. Total records: ",
                          object$geom.full$records,". Unique locations:")
    }
    output[4, 1] <- object$geom.full$locations
  }
  if (!is.null(object$geom.core)) {
    output[1, 1] <- "Full"
    output[1, 2] <- "Core"
    clustersAvailableCore <- with(object$geom.core, floor(locations/mean_h))
    output[4, 2] <- object$geom.core$locations
    output[5, 2] <- clustersAvailableCore
    output[6, 2] <- round(object$geom.core$mean_h, digits = 1)
    output[7, 2] <- round(object$geom.core$sd_h, digits = 1)
  }
  if (!is.null(object$trial$arm)) {
    sd1 <- ifelse(is.null(object$geom.full$sd_distance), NA, object$geom.full$sd_distance)
    sd2 <- ifelse(is.null(object$geom.core$sd_distance), NA, object$geom.core$sd_distance)
    rownames(output)[8] <- "S.D. of distance to nearest discordant location (km): "
    output[8, 1] <- ifelse(is.na(sd1), "", round(sd1, digits = 2))
    output[8, 2] <- ifelse(is.na(sd2), "", round(sd2, digits = 2))
    rownames(output)[9] <- "Cluster randomization:            "
    if (is.null(object$trial$pair)) {
      output[9, 1] <- "Independently randomized"
    } else {
      output[9, 1] <- "Matched pairs randomized"
    }
  } else {
    if (is.null(object$trial$x)) {
      rownames(output)[9] <- "No locations to randomize"
    } else {
      rownames(output)[9] <- "No randomization"
    }
    output[9, 1] <- "-"
  }
  if (!is.null(object$design)) {
    rownames(output)[10] <- "\nSpecification of Requirements\n-----------------------------"
    output[10, 1] <- "-"
    rownames(output)[11] <- "Significance level (2-sided):    "
    output[11, 1] <- object$design$alpha
    rownames(output)[12] <- "Type of Outcome                   "
    output[12, 1] <- switch(object$design$outcome.type,
                   'y' = "continuous",
                   "n" = "count",
                   "e" = "event rate",
                   'p' = "proportion",
                   'd' = "dichotomous")
    rownames(output)[13] <- "Expected outcome in control arm:  "
    output[13, 1] <- object$design$yC
    link <- switch(object$design$outcome.type,
                   'y' = "identity",
                   "n" = "log",
                   "e" = "log",
                   'p' = "logit",
                   'd' = "logit")
    rownames(output)[14] <- switch(link,
            "identity" = "Expected variance of outcome:     ",
                 "log" = "Expected overdispersion:          ",
               "logit" = "Mean denominator:                 ")
    output[14, 1] <- switch(link,
                    "identity" =  object$design$sigma2,
                    "log" = object$design$phi,
                    "logit" = object$design$N)
    rownames(output)[15] <- "Required effect size:             "
    output[15, 1] <- object$design$effect
    rownames(output)[16] <- "Intra-cluster correlation:        "
    output[16, 1] <- object$design$ICC
    if (!is.null(object$design$buffer.width)) {
      rownames(output)[3] <- "Buffer width :               "
      if (object$design$buffer.width > 0) {
        output[3, 1] <- paste0(object$design$buffer.width,
                               " km.")
      } else {
        output[3, 1] <- "No buffer"
      }
    }
  }

  output[17, 1] <- "-"
  if (is.null(object$design)) {
    rownames(output)[17] <- "No power calculations to report"
  } else {
    rownames(output)[17] <- "\nPower calculations (ignoring contamination)\n------------------                    "
    sufficient <- ifelse(clustersAvailableFull >= object$geom.full$clustersRequired,
                         "Yes", "No")
    rownames(output)[18] <- "Design effect:                         "
    output[18, 1] <- round(object$geom.full$DE, digits = 1)
    rownames(output)[19] <- "Nominal power (%)                      "
    output[19, 1] <- round(object$geom.full$power * 100, digits = 1)
    rownames(output)[20] <- paste0("Clusters required for power of ",
                                   object$design$desiredPower * 100, "%:     ")
    output[20, 1] <- object$geom.full$clustersRequired
    rownames(output)[21] <- "Sufficient clusters for required power?"
    output[21, 1] <- sufficient

    if (is.null(object$geom.core)) {
      output <- subset(output, select = -c(2))
    } else {
      output[17, 1] <- "Full"
      output[17, 2] <- "Core"
      clustersAvailableCore <- with(object$geom.core, floor(locations/mean_h))
      sufficientCore <- ifelse(clustersAvailableCore >= object$geom.core$clustersRequired,
                               "Yes", "No")

      output[18, 2] <- round(object$geom.core$DE, digits = 1)
      output[19, 2] <- round(object$geom.core$power * 100, digits = 1)
      output[20, 2] <- object$geom.core$clustersRequired
      output[21, 2] <- sufficientCore
    }
  }
  standard.names <- c("x", "y", "cluster", "arm", "buffer", "nearestDiscord",
                      "geom.full", "geom.core", "design")
  rownames(output)[22] <- "\nOther variables in dataset\n--------------------------"
  output[22, 1] <- paste(dplyr::setdiff(names(object$trial), standard.names), collapse = "  ")
  output <- output[trimws(output[, 1]) != "", ]
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

