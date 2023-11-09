#' Aggregate data across records with duplicated locations
#'
#' \code{aggregateCRT} aggregates data from a \code{"CRTsp"} object or trial data frame containing multiple records with the same location,
#' and outputs a list of class \code{"CRTsp"} containing single values for each location, for both the coordinates and the auxiliary variables.
#' @param trial An object of class \code{"CRTsp"} containing locations (x,y) and variables to be summed
#' @param auxiliaries vector of names of auxiliary variables to be summed across each location
#' @returns A list of class \code{"CRTsp"}
#' @details
#' Variables that in the trial dataframe that are not included in \code{auxiliaries} are retained in the output
#' algorithm \code{"CRTsp"} object, with the value corresponding to that of the first record for the location
#' in the input data frame
#' @examples {
#' trial <- readdata('example_site.csv')
#' trial$base_denom <- 1
#' aggregated <- aggregateCRT(trial, auxiliaries = c("RDT_test_result","base_denom"))
#' }
#' @export
#'
aggregateCRT <- function(trial, auxiliaries = NULL) {
    CRT <- CRTsp(trial)
    location <- NULL
    trial <- CRT$trial[order(CRT$trial$x, CRT$trial$y),]
    trial$location <- paste(trial$x,trial$y)
    trial1 <- trial
    if (length(auxiliaries) > 0) {
      auxvars <- names(trial) %in% auxiliaries
      trial1 <- dplyr::distinct(trial, location, .keep_all = TRUE)
      trial1 <- trial1[, !auxvars]
      # This code is a mess, but the smarter options seem to be work in progress in dplyr
      for(var in auxiliaries){
        if(var %in% names(trial)) {
          trial2 <- with(trial,
                         trial %>%
                         dplyr::group_by(location) %>%
                         dplyr::summarize(var = sum(get(var))))
          class(trial2) <- "data.frame"
          colnames(trial2) <- c('location',var)
          trial1 <- merge(trial1, trial2, by = 'location', all.x = FALSE, all.y = FALSE)
        } else {
          message('*** Variable', var,' not present in input data ***')
        }
      }
    }
    trial1$location <- NULL
    CRT$trial <- trial1
    return(CRTsp(CRT))
}

#' Specification of buffer zone in a cluster randomized trial
#'
#' \code{specify_buffer} specifies a buffer zone in a cluster randomized
#' trial (CRT) by flagging those locations that are within a defined distance of
#' those in the opposite arm.
#'
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}).
#' @param buffer_width minimum distance between locations in
#'   opposing arms for them to qualify to be included in the core area (km)
#' @returns A list of class \code{"CRTsp"} containing the following components:
#'  \tabular{lll}{
#'  \code{geom_full}   \tab list: \tab summary statistics describing the site,
#'  cluster assignments, and randomization.\cr
#'  \code{geom_core}   \tab list: \tab summary statistics describing the core area \cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\cr
#'  \tab \code{x} \tab numeric vector: x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor: assignments to cluster of each location  \cr
#'  \tab \code{arm} \tab factor: assignments to \code{"control"} or \code{"intervention"} for each location \cr
#'  \tab \code{nearestDiscord} \tab numeric vector: signed Euclidean distance to nearest discordant location (km) \cr
#'  \tab \code{buffer} \tab logical: indicator of whether the point is within the buffer \cr
#'  \tab \code{...} \tab other objects included in the input \code{"CRTsp"} object or data frame \cr
#'  }
#' @export
#' @examples
#' #Specify a buffer of 200m
#' exampletrial <- specify_buffer(trial = readdata('exampleCRT.txt'), buffer_width = 0.2)
specify_buffer <- function(trial, buffer_width = 0) {
  CRT <- CRTsp(trial)
  trial <- CRT$trial
  if (is.null(trial$arm)) return('*** Randomization is required before buffer specification ***')
  if (is.null(trial$nearestDiscord)) trial <- compute_distance(trial, distance = "nearestDiscord")
  if (buffer_width > 0) {
    trial$buffer <- (abs(trial$nearestDiscord) < buffer_width)
  }
  CRT$trial <- trial
  return(CRTsp(CRT))
}

#' Randomize a two-armed cluster randomized trial
#'
#' \code{randomizeCRT} carries out randomization of clusters for a CRT and
#' augments the trial dataframe with assignments to arms \cr
#'
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}). Optionally: specification of a buffer zone (logical \code{buffer});
#'   any other variables required for subsequent analysis.
#' @param matchedPair logical: indicator of whether pair-matching on the
#'   baseline data should be used in randomization
#' @param baselineNumerator name of numerator variable for baseline data (required for
#'   matched-pair randomization)
#' @param baselineDenominator name of denominator variable for baseline data (required for
#'   matched-pair randomization)
#' @returns A list of class \code{"CRTsp"} containing the following components:
#'  \tabular{lll}{
#'  \code{design}   \tab list: \tab parameters required for power calculations\cr
#'  \code{geom_full}   \tab list: \tab summary statistics describing the site\cr
#'  \code{geom_core}   \tab list: \tab summary statistics describing the core area
#'  (when a buffer is specified)\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\cr
#'  \tab \code{x} \tab numeric vector: x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor: assignments to cluster of each location  \cr
#'  \tab \code{pair} \tab factor: assigned matched pair of each location
#'  (for \code{matchedPair} randomisations) \cr
#'  \tab \code{arm} \tab factor: assignments to \code{"control"} or \code{"intervention"} for each location \cr
#'  \tab \code{...} \tab other objects included in the input \code{"CRTsp"} object or data frame \cr
#'  }
#' @export
#' @examples
#' # Randomize the clusters in an example trial
#' exampleCRT <- randomizeCRT(trial = readdata('exampleCRT.txt'), matchedPair = TRUE)
randomizeCRT <- function(trial, matchedPair = FALSE, baselineNumerator = "base_num",
    baselineDenominator = "base_denom") {

    CRT <- CRTsp(trial)
    CRT$design <- NULL
    trial <- CRT$trial

    # remove any preexisting assignments and coerce matchedPair to FALSE if there are no baseline data
    if(is.null(trial[[baselineNumerator]]) & matchedPair) {
        warning("*** No baseline data for matching. Unmatched randomisation ***")
        matchedPair <- FALSE
    }
    trial$arm <- trial$pair <- trial$nearestDiscord <- trial$hdep <- trial$sdep <- trial$disc <- trial$kern <- NULL
    pair <- cluster <- base_num <- base_denom <- NULL

    trial$cluster <- as.factor(trial$cluster)
    # Randomization, assignment to arms
    nclusters <- length(unique(trial$cluster))
    if ((nclusters%%2) == 1 & matchedPair) {
        warning("*** odd number of clusters: assignments are not matched on baseline data ***")
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
    CRT$trial <- trial
    CRT <- compute_distance(CRT, distance = "nearestDiscord")
    return(CRTsp(CRT))
}


plt <- function(object) {
  UseMethod("plt")
}



simulate_site <- function(geoscale, locations, kappa, mu) {
  scaling = geoscale * 10
  # Poisson point pattern with Thomas algorithm
  p <- spatstat.random::rThomas(kappa, geoscale, mu, win = spatstat.geom::owin(c(0, scaling), c(0, scaling)))
  # expected number of points: kappa*mu*scaling^2

  # create locations and specify co-ordinates
  hhID <- c(1:locations)
  x <- p$x[seq(1:locations)]
  y <- p$y[seq(1:locations)]
  coordinates <- data.frame(x = x - mean(x), y = y - mean(y))
  trial <- coordinates
  return(trial)
}


#' Algorithmically assign locations to clusters in a CRT
#'
#' \code{specify_clusters} algorithmically assigns locations to clusters by grouping them geographically
#'
#' @param trial A CRT object or data frame containing (x,y) coordinates of
#'   households
#' @param c integer: number of clusters in each arm
#' @param h integer: number of locations per cluster
#' @param algorithm algorithm for cluster boundaries, with options:
#' \tabular{ll}{
#' \code{NN}\tab Nearest neighbour: assigns equal numbers of locations to each cluster \cr
#' \code{kmeans}\tab kmeans clustering: aims to partition locations so that each
#' belongs to the cluster with the nearest centroid.\cr
#' \code{TSP}\tab travelling salesman problem heuristic: Assigns locations sequentially
#' along a travelling salesman path.\cr
#' }
#' @param reuseTSP logical: indicator of whether a pre-existing path should be used by
#'   the TSP algorithm
#' @returns A list of class \code{"CRTsp"} containing the following components:
#'  \tabular{lll}{
#'  \code{geom_full}   \tab list: \tab summary statistics describing the site,
#'  and cluster assignments.\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\cr
#'  \tab \code{x} \tab numeric vector: x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor: assignments to cluster of each location  \cr
#'  \tab \code{...} \tab other objects included in the input \code{"CRTsp"} object or data frame \cr
#'  }
#' @details
#' The \code{reuseTSP} parameter is used to allow the path to be reused
#' for creating alternative allocations with different cluster sizes.\cr\cr
#' Either \code{c} or \code{h} must be specified. If both are specified
#' the input value of \code{c} is ignored.\cr
#' @export
#'
#' @examples
#' #Assign clusters of average size h = 40 to a test set of co-ordinates, using the kmeans algorithm
#' exampletrial <- specify_clusters(trial = readdata('exampleCRT.txt'),
#'                             h = 40, algorithm = 'kmeans', reuseTSP = FALSE)
specify_clusters <- function(trial = trial, c = NULL, h = NULL, algorithm = "NN",
    reuseTSP = FALSE) {

    CRT <- CRTsp(trial)
    trial <- CRT$trial

    # Local data from study area (ground survey and/or satellite
    # images)
    coordinates <- data.frame(x = as.numeric(as.character(trial$x)), y = as.numeric(as.character(trial$y)))

    # the number of clusters and the target cluster size must be integers.
    # cluster size can only be exactly equal to the input value of h if this is a factor of
    # the number of locations
    if (is.null(c)) {
        c <- ceiling(nrow(coordinates)/(2 * h))
    }
    if (is.null(h)) {
        h <- ceiling(nrow(coordinates)/(2 * c))
    }
    nclusters <- 2 * c
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

    # remove any pre-existing arm assignments
    trial$arm <- NULL
    CRT$trial <- trial
    return(CRTsp(CRT))
}


#' Convert lat long co-ordinates to x,y
#'
#' \code{latlong_as_xy} converts co-ordinates expressed as decimal degrees into x,y
#' @param trial A trial dataframe or list of class \code{"CRTsp"} containing latitudes and longitudes in decimal degrees
#' @param latvar name of column containing latitudes in decimal degrees
#' @param longvar name of column containing longitudes in decimal degrees
#' @details The output object contains the input locations replaced with Cartesian
#'   coordinates in units of km, centred on (0,0), corresponding to using the equirectangular projection
#'   (valid for small areas). Other data are unchanged.
#' @returns A list of class \code{"CRTsp"} containing the following components:
#'  \tabular{lll}{
#'  \code{geom_full}   \tab list: \tab summary statistics describing the site \cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\cr
#'  \tab \code{x} \tab numeric vector: x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: y-coordinates of locations \cr
#'  \tab \code{...} \tab other objects included in the input \code{"CRTsp"} object or data frame \cr
#'  }
#' @export
#' @examples
#' examplexy <- latlong_as_xy(readdata("example_latlong.csv"))
#'
latlong_as_xy <- function(trial, latvar = "lat", longvar = "long") {
  CRT <- CRTsp(trial)
  trial <- CRT$trial
  colnames(trial)[colnames(trial) == latvar] <- "lat"
  colnames(trial)[colnames(trial) == longvar] <- "long"
  # scalef is the number of degrees per kilometer
  scalef <- 180/(6371*pi)
  centroid <- list(lat = mean(trial$lat),
                   long = mean(trial$long))
  trial$y <- (trial$lat - centroid$lat)/scalef
  trial$x <- (trial$long - centroid$long) * cos(trial$lat * pi/180)/scalef
  drops <- c("lat", "long")
  trial <- trial[, !(names(trial) %in% drops)]
  CRT <- CRTsp(trial, design = NULL)
  CRT$geom_full$centroid <- centroid
  return(CRT)
}


#' Anonymize locations of a trial site
#'
#' \code{anonymize_site} transforms coordinates to remove potential identification information.
#' @param trial \code{"CRTsp"} object or trial data frame with co-ordinates of households
#' @param ID name of column used as an identifier for the points
#' @param latvar name of column containing latitudes in decimal degrees
#' @param longvar name of column containing longitudes in decimal degrees
#' @returns A list of class \code{"CRTsp"}.
#' @export
#' @details
#' The coordinates are transformed to support confidentiality of
#' information linked to households by replacing precise geo-locations with transformed co-ordinates which preserve distances
#' but not positions. The input may have either \code{lat long} or \code{x,y} coordinates.
#' The function first searches for any \code{lat long} co-ordinates and converts these to \code{x,y}
#' Cartesian coordinates. These are then are rotated by a random angle about a random origin. The returned object
#' has transformed co-ordinates re-centred at the origin. Centroids stored in the \code{"CRTsp"} object are removed.
#' Other data are unchanged.
#' @examples
#' #Rotate and reflect test site locations
#' transformedTestlocations <- anonymize_site(trial =  readdata("exampleCRT.txt"))
anonymize_site <- function(trial, ID = NULL, latvar = "lat", longvar = "long") {
    # Local data from study area (ground survey and/or satellite
    # images) random rotation angle
    CRT <- CRTsp(trial)
    trial <- CRT$trial
    if (latvar %in% colnames(trial)) {
      CRT <- latlong_as_xy(trial, latvar = latvar, longvar = longvar)
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

    # Remove ID variable
    if (!is.null(ID)) {
      colnames(trial)[colnames(trial) == ID] <- "ID"
      trial$ID <- NULL
    }

    # Remove centroid information
    CRT$geom_full$centroid <- NULL
    CRT$trial <- trial
    return(CRTsp(CRT))
}


#' Read example dataset
#'
#' \code{readdata} reads a file from the package library of example datasets
#'
#' @param filename name of text file stored within the package
#' @return R object corresponding to the text file
#' @details The input file name should include the extension (either .csv or .txt).
#' The resulting object is a data frame if the extension is .csv.
#' @export
#' @examples
#' exampleCRT <- readdata('exampleCRT.txt')
#'
readdata <- function(filename) {
    fname <- eval(filename)
    extdata <- system.file("extdata", package = "CRTspat")
    if (unlist(gregexpr("mesh", fname)) > 0) {
      # The mesh was stored using saveRDS e.g.
      # library(Matrix)
      # saveRDS(inla_mesh,file = "inst/extdata/examplemesh100.rds")
      robject <- readRDS(file = paste0(extdata, "/", fname))
    } else if (unlist(gregexpr("analysis", fname)) > 0) {
      # Analysis objects should be stored using 'dump' but are easy to reproduce
      sourced <- load(file = paste0(extdata, "/", fname))
      robject <- sourced$value
    } else if (unlist(gregexpr(".csv", fname)) > 0) {
      robject <- read.csv(file = paste0(extdata, "/", fname), row.names = NULL)
      # remove variable 'X' if it is present
      robject$X <- NULL
    } else if (unlist(gregexpr(".txt", fname)) > 0) {
      sourced <- source(file = paste0(extdata, "/", fname))
      robject <- sourced$value
    }
    if (unlist(gregexpr("CRT", fname)) > 0) robject <- CRTsp(robject)
    return(robject)
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

map_scale_to_link <- function(scale) {
  scales <- c("proportion", "count", "continuous")
  links <- c("logit", "log", "identity")
  link <-  links[which(scale == scales)]
return(link)}
