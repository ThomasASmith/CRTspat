#' Create or update a \code{"CRTsp"} object
#'
#' \code{CRTsp} coerces data frames containing co-ordinates and location attributes
#' into objects of class \code{"CRTsp"} or creates a new \code{"CRTsp"} object by simulating a set of Cartesian co-ordinates for use as the locations in a simulated trial site
#' @param x an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}). Optionally specification of a buffer zone (logical \code{buffer});
#'   any other variables required for subsequent analysis.
#' @param design list: an optional list containing the requirements for the power of the trial
#' @param geoscale numeric: standard deviation of random displacement from each settlement cluster center (for new objects)
#' @param locations integer: number of locations in population (for new objects)
#' @param kappa numeric: intensity of Poisson process of settlement cluster centers (for new objects)
#' @param mu numeric: mean number of points per settlement cluster (for new objects)
#' @param geometry  with valid values \code{'point'} (the default, corresponding to point locations), \code{'triangle'},
#' \code{'square'} and \code{'hexagon'} corresponding to grids constructed from pixels of regular polygons.
#' @export
#' @returns A list of class \code{"CRTsp"} containing the following components:
#'  \tabular{lll}{
#'  \code{design}   \tab list: \tab parameters required for power calculations\cr
#'  \code{geom_full}   \tab list: \tab summary statistics describing the site \cr
#'  \code{geom_core}   \tab list: \tab summary statistics describing the core area
#'  (when a buffer is specified)\cr
#'  \code{trial} \tab data frame: \tab rows correspond to geolocated points, as follows:\cr
#'  \tab \code{x} \tab numeric vector: x-coordinates of locations \cr
#'  \tab \code{y} \tab numeric vector: y-coordinates of locations \cr
#'  \tab \code{cluster} \tab factor: assignments to cluster of each location  \cr
#'  \tab \code{arm} \tab factor: assignments to \code{"control"} or \code{"intervention"} for each location \cr
#'  \tab \code{nearestDiscord} \tab numeric vector: Euclidean distance to nearest discordant location (km) \cr
#'  \tab \code{buffer} \tab logical: indicator of whether the point is within the buffer \cr
#'  \tab \code{...} \tab other objects included in the input \code{"CRTsp"} object or data frame \cr
#'  }
#' @details
#' If a data frame or \code{"CRTsp"} object is input then the output \code{"CRTsp"} object is validated,
#' a description of the geography is computed and power calculations are carried out.\cr\cr
#' If \code{geoscale, locations, kappa} and \code{mu} are specified then a new trial dataframe is constructed
#' corresponding to a novel simulated human settlement pattern. This is generated using the
#' Thomas algorithm (\code{rThomas}) in [\code{spatstat.random}](https://CRAN.R-project.org/package=spatstat.random)
#' allowing the user to defined the density of locations and degree of spatial clustering.
#' The resulting trial data frame comprises a set of Cartesian coordinates centred at the origin.\cr\cr
#' @export
#' @examples
#' {# Generate a simulated area with 10,000 locations
#' example_area = CRTsp(geoscale = 1, locations=10000, kappa=3, mu=40)
#' summary(example_area)
#' }
CRTsp <- function(x = NULL, design = NULL, geoscale = NULL, locations = NULL,
                  kappa = NULL, mu = NULL, geometry = 'point') {
  centroid <- list(lat = NULL, long = NULL)
  if(identical(class(x),"CRTsp")) {
    CRT <- x
    if(!is.null(design)) CRT$design <- design
    centroid <- if(!is.null(CRT$geom_full$centroid$lat)) CRT$geom_full$centroid
  } else if("data.frame" %in% class(x)) {
    CRT <- list(trial = data.frame(x), design = design)
  } else if(is.null(x)) {
    if (!is.null(geoscale) & !is.null(locations)
        & !is.null(kappa) & !is.null(mu)){
      trial <- simulate_site(geoscale = geoscale, locations=locations, kappa=kappa, mu=mu)
      CRT <- list(trial = trial, design = design)
    } else {
      warning("*** All of geoscale, locations, kappa, mu needed to simulate a new site ***")
      CRT <- list(trial = data.frame(x=numeric(0),y=numeric(0)), design = design)
    }
  }
  if(is.null(CRT$design)) CRT$design <- list(locations = NULL, geometry = geometry,
                                             alpha = NULL, desiredPower = NULL, effect = NULL, yC = NULL,
                                             outcome_type = NULL, sigma2 = NULL, denominator = NULL, N = NULL,
                                             ICC = NULL, k = NULL, d_h = NULL, spillover_interval = NULL,
                                             contaminate_pop_pr = 0, distance_distribution = NULL)
  if(is.null(CRT$trial)) CRT$trial <- data.frame(x=numeric(0),y=numeric(0))
  CRT$geom_full <- get_geom(trial = CRT$trial, design = CRT$design)
  CRT$geom_full$centroid <- centroid
  empty_geom <- list(locations = 0, sd_h = NULL, k= NULL, uniquelocations = 0, mean_h = NULL,
    DE = NULL, power = NULL, clustersRequired = NULL)
  if (is.null(CRT$trial$buffer)) {
    CRT$geom_core <- empty_geom
  } else {
    if (all(CRT$trial$buffer == TRUE)) {
      CRT$geom_core <- empty_geom
    } else {
      CRT$geom_core <- get_geom(trial = CRT$trial[CRT$trial$buffer == FALSE, ],
                                design = CRT$design)
    }
  }
  return(validate_CRTsp(new_CRTsp(CRT)))
}

is_CRTsp <- function(x) {
  return(inherits(x, "CRTsp"))
}

new_CRTsp <- function(x = list()) {
  stopifnot(is.data.frame(x$trial))
  stopifnot(is.list(x$design))
  stopifnot(is.list(x$geom_full))
  stopifnot(is.list(x$geom_core))
  return(structure(x, class = "CRTsp"))
}

validate_CRTsp <- function(x) {
  stopifnot(inherits(x, "CRTsp"))
  values <- unclass(x)
  if (is.null(values$trial) & is.null(values$design)) {
    stop("There must be either a design or a trial data frame in `x`")
  }
  if (!is.null(values$trial)){
    if (!identical(class(values$trial),"data.frame")){
      stop("The trial object in `x` must be a data frame")
    }
    if (nrow(values$trial) > 0) {
      if (nrow(values$trial) != values$geom_full$locations){
        stop("The geom_full object in `x` is invalid")
      }
    }
  }
  return(x)
}

#' Summary description of a \code{"CRTsp"} object
#'
#' \code{summary.CRTsp} provides a description of a \code{"CRTsp"} object
#' @param object an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), arm assignments (factor \code{arm}) and buffer zones (logical \code{buffer}),
#'   together with any other variables required for subsequent analysis.
#' @param maskbuffer radius of area around a location to include in calculation of areas
#' @param ... other arguments used by summary
#' @method summary CRTsp
#' @export
#' @return No return value, write text to the console.
#' @examples
#' summary(CRTsp(readdata('exampleCRT.txt')))
summary.CRTsp <- function(object, maskbuffer = 0.2, ...) {
  defaultdigits <- getOption("digits")
  on.exit(options(digits = defaultdigits))
  options(digits = 3)
  cat("===============================CLUSTER RANDOMISED TRIAL ===========================\n")
  output <- matrix("  ", nrow = 23, ncol = 3)
  required_cols <- c(1, 2, 3)
  # The third column is to avoid loss of row names when the object is reduced to a vector
  rownames(output) <- paste0("row ", 1:nrow(output))
  rownames(output)[1] <- "Locations and Clusters\n----------------------                                "
  output[1, 1] <- "-"
  rownames(output)[2] <- "Coordinate system            "
  if (length(object$trial$x) > 0 & length(object$trial$y) > 0) {
    output[2, 1] <- "(x, y)"
  } else if(!is.null(object$trial$lat)) {
    output[2, 1] <- "Lat-Long"
  } else {
    output[2, 1] <- "No coordinates in dataset"
  }

  if (identical(unname(output[2, 1]),"(x, y)")) {
    cat("\nSummary of coordinates\n----------------------\n")
    if(is.null(object$trial$nearestDiscord)) {
      coordinate.summary <- with(object$trial, summary(cbind(x, y)))
    } else {
      coordinate.summary <- with(object$trial, summary(cbind(x, y, nearestDiscord)))
    }
    rownames(coordinate.summary) <- substr(coordinate.summary[, 1], 1, 8)
    coordinate.summary[, ] <- substr(coordinate.summary[, ], 9, 13)
    print(t(coordinate.summary))
    xycoords <- data.frame(cbind(x=object$trial$x,y=object$trial$y))
    tr <- sf::st_as_sf(xycoords, coords = c("x","y"))
    buf1 <- sf::st_buffer(tr, maskbuffer)
    buf2 <- sf::st_union(buf1)
    area <- sf::st_area(buf2)
    cat("\nTotal area (within ", maskbuffer,"km of a location) : ", area, "sq.km\n")
    convex_hull <- sf::st_area(sf::st_convex_hull(sf::st_union(tr)))
    cat("Total area (convex hull) : ", convex_hull, "sq.km\n\n")
    if (!is.null(object$geom_full$centroid$lat)) {
      cat("Geolocation of centroid (radians):
          latitude: ", object$geom_full$centroid$lat,
          "longitude: ", object$geom_full$centroid$long,"\n\n")
    }
  }
  rownames(output)[5] <- "Available clusters (across both arms)                 "
  if (is.na(object$geom_full$c)) {
    output[5, 1] <- "Not assigned"
  } else {
    clustersAvailableFull <- with(object$geom_full, floor(locations/mean_h))
    output[5, 1] <- clustersAvailableFull
    rownames(output)[6] <- "  Per cluster mean number of points                   "
    output[6, 1] <- round(object$geom_full$mean_h, digits = 1)
    rownames(output)[7] <- "  Per cluster s.d. number of points                   "
    if (!is.null(object$geom_full$sd_h))
      output[7, 1] <- round(object$geom_full$sd_h, digits = 1)
  }
  rownames(output)[4] <- "Locations:                                            "

  if(identical(object$geom_full$locations,object$geom_full$uniquelocations)){
    output[4, 1] <- object$geom_full$locations
  } else {
    if(!is.null(object$geom_full$locations)) {
      rownames(output)[4] <- paste0("Not aggregated. Total records: ",
                                    object$geom_full$locations,". Unique locations:")
    }
    if(!is.null(object$geom_full$uniquelocations))
      output[4, 1] <- object$geom_full$uniquelocations
  }
  if (object$geom_core$locations > 0 &
      object$geom_core$locations < object$geom_full$locations) {
    output[1, 1] <- "Full"
    output[1, 2] <- "Core"
    output[4, 2] <- object$geom_core$uniquelocations
    if (!is.na(object$geom_core$mean_h)) {
      clustersAvailableCore <- with(object$geom_core, floor(locations/mean_h))
      output[5, 2] <- clustersAvailableCore
      output[6, 2] <- round(object$geom_core$mean_h, digits = 1)
    }
    if (!is.null(object$geom_core$sd_h))
      output[7, 2] <- round(object$geom_core$sd_h, digits = 1)
  }
  if (!is.null(object$trial$arm) & !identical(object$trial$arm,character(0))) {
    sd1 <- ifelse(is.null(object$geom_full$sigma_x), NA, object$geom_full$sigma_x)
    sd2 <- ifelse(is.null(object$geom_core$sigma_x), NA, object$geom_core$sigma_x)
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
      rownames(output)[8] <- "    "
      rownames(output)[9] <- "No locations to randomize"
    } else {
      rownames(output)[9] <- "No randomization"
    }
    output[9, 1] <- "-"
  }
  if (!is.null(object$design$alpha)) {
    rownames(output)[10] <- "\nSpecification of Requirements\n-----------------------------"
    output[10, 1] <- "-"
    rownames(output)[11] <- "Significance level (2-sided):    "
    output[11, 1] <- object$design$alpha
    rownames(output)[12] <- "Type of Outcome                   "
    output[12, 1] <- switch(object$design$outcome_type,
                            'y' = "continuous",
                            "n" = "count",
                            "e" = "event rate",
                            'p' = "proportion",
                            'd' = "dichotomous")
    rownames(output)[13] <- "Expected outcome in control arm:  "
    output[13, 1] <- object$design$yC
    link <- switch(object$design$outcome_type,
                   'y' = "identity",
                   "n" = "log",
                   "e" = "log",
                   'p' = "logit",
                   'd' = "logit")
    rownames(output)[14] <- switch(link,
                                   "identity" = "Expected variance of outcome:     ",
                                   "log" = "Mean rate multiplier:             ",
                                   "cloglog" = "Mean rate multiplier:             ",
                                   "logit" = "Mean denominator:                 ")
    output[14, 1] <- switch(link,
                            "identity" =  object$design$sigma2,
                            "log" = object$design$denominator,
                            "cloglog" = object$design$denominator,
                            "logit" = object$design$N)
    if (identical(object$design$outcome_type, 'd')) output[14, 1] <- ""
    rownames(output)[15] <- "Required effect size:             "
    output[15, 1] <- object$design$effect
    if (is.na(object$design$ICC)) {
      rownames(output)[16] <- "Coefficient of variation (%):     "
      output[16, 1] <- object$design$cv_percent
    } else {
      rownames(output)[16] <- "Intra-cluster correlation:        "
      output[16, 1] <- object$design$ICC
    }
    if (!is.null(object$design$buffer_width)) {
      rownames(output)[3] <- "Buffer width :               "
      if (object$design$buffer_width > 0) {
        output[3, 1] <- paste0(object$design$buffer_width,
                               " km.")
      } else {
        output[3, 1] <- "No buffer"
      }
    } else {
      rownames(output)[3] <- " "
    }
  }

  output[17, 1] <- "-"
  if (is.null(object$design$effect)) {
    rownames(output)[17] <- "No power calculations to report"
  } else {
    rownames(output)[17] <- "\nPower calculations\n------------------                    "

    sufficient <- ifelse(clustersAvailableFull >= object$geom_full$clustersRequired,
                         "Yes", "No")
    rownames(output)[18] <- "Design effect:                         "
    output[18, 1] <- round(object$geom_full$DE, digits = 1)
    if (object$geom_full$contaminate_pop_pr > 0) {
      rownames(output)[19] <- paste0("Spillover affecting ",
        round(object$geom_full$contaminate_pop_pr*100, digits = 1),
      "% of data,\n  ", object$design$distance_distribution," model gives bias estimate:")
      output[19, 1] <- round(object$geom_full$delta, digits = 3)
    } else {
      rownames(output)[19] <- "Calculations ignoring spillover"
    }
    rownames(output)[20] <- "Nominal power (%)                      "
    output[20, 1] <- round(object$geom_full$power * 100, digits = 1)
    rownames(output)[21] <- paste0("Total clusters required (power of ",
                                   object$design$desiredPower * 100, "%):")
    output[21, 1] <- object$geom_full$clustersRequired
    rownames(output)[22] <- "Sufficient clusters for required power?"
    output[22, 1] <- sufficient

    if (is.null(object$geom_core$power)) {
      required_cols <- c(1, 3)
    } else {
      required_cols <- c(1, 2)
      output[17, 1] <- "Full"
      output[17, 2] <- "Core"
      clustersAvailableCore <- with(object$geom_core, floor(locations/mean_h))
      sufficientCore <- ifelse(clustersAvailableCore >= object$geom_core$clustersRequired,
                               "Yes", "No")
      output[18, 2] <- round(object$geom_core$DE, digits = 1)
      if (object$geom_core$contaminate_pop_pr > 0) {
        output[19, 2] <- round(object$geom_core$delta, digits = 3)
      }
      output[20, 2] <- round(object$geom_core$power * 100, digits = 1)
      output[21, 2] <- object$geom_core$clustersRequired
      output[22, 2] <- sufficientCore
    }
  }
  standard_names <- c("x", "y", "cluster", "arm", "buffer", "nearestDiscord",
                      "geom_full", "geom_core", "design")
  rownames(output)[23] <- "\nOther variables in dataset\n--------------------------"
  output[23, 1] <- paste(dplyr::setdiff(names(object$trial), standard_names), collapse = "  ")
  required_rows <- seq(1:nrow(output))[nchar(output[, 1]) > 0 & substr(rownames(output), 1, 3) != 'row']
  output1 <- output[required_rows, required_cols]
  # display and return table
  utils::write.table(output1, quote = FALSE, col.names = FALSE, sep = "          ")
  options(digits = defaultdigits)
  invisible(object)
}

