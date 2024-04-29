#' Compute distance or surround values for a cluster randomized trial
#'
#' \code{compute_distance} computes distance or surround values for a cluster randomized trial (CRT)
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}).
#' @param distance the quantity(s) to be computed. Options are:
#'  \tabular{ll}{
#' \code{"nearestDiscord"} \tab distance to nearest discordant location (km)\cr
#' \code{"distanceAssigned"} \tab distance to the nearest pixel in the assigned cluster (km)\cr
#' \code{"disc"} \tab disc \cr
#' \code{"kern"} \tab kernel-based measure \cr
#' \code{"hdep"} \tab Tukey half space depth\cr
#' \code{"sdep"} \tab simplicial depth\cr
#' }
#' @param scale_par scale parameter equal to the disc radius in km if \code{distance = "disc"}
#' or to the standard deviance of the kernels if \code{distance = "kern"}
#' @param auxiliary \code{"CRTsp"} object containing external cluster and or arm assignments.
#' @returns The input \code{"CRTsp"} object with additional column(s) added to the \code{trial} data frame
#' with variable name corresponding to the input value of \code{distance}.
#' @details
#' For each selected distance measure, the function first checks whether the variable is already present, and carries out
#' the calculations only if the corresponding field is absent from the \code{trial} data frame.\cr\cr
#' If \code{distance = "nearestDiscord"} is selected the computed values are Euclidean distances
#' assigned a positive sign for the intervention arm of the trial, and a negative sign for the control arm.\cr\cr
#' If \code{distance = "distanceAssigned"} is selected the computed values are Euclidean distances
#' to the nearest pixel in the \code{auxiliary} \code{"CRTsp"} object.\cr\cr
#' If \code{distance = "disc"} is specified, the disc statistic is computed for each location as the number of locations
#' within the specified radius that are in the intervention arm
#' ([Anaya-Izquierdo & Alexander(2020)](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.13316)). The input
#' value of \code{scale_par} is stored in the \code{design} list
#' of the output \code{"CRTsp"} object. Recalculation is carried out if the input value of
#' \code{scale_par} differs from the one in the input \code{design} list. The value of the the surround calculated
#' based on intervened locations is divided by the value of the surround calculated on the basis of all locations, so the
#' value returned is a proportion.\cr\cr
#' If \code{distance = "kern"} is specified, the Normal curve with standard deviation
#' \code{scale_par} is used to simulate diffusion of the intervention effect by Euclidean
#' distance. For each location in the trial, the contributions of all intervened locations are
#' summed. As with \code{distance = "disc"}, when \code{distance = "kern"} the surround calculated
#' based on intervened locations is divided by the value of the surround calculated on the basis of all locations, so the
#' value returned is a proportion.\cr\cr
#' If either \code{distance = "hdep"} or \code{distance = "sdep"} is specified then both the simplicial depth and
#' Tukey half space depth are calculated using the algorithm of
#' [Rousseeuw & Ruts(1996)](https://www.jstor.org/stable/2986073). The half-depth probability within the intervention cloud (di) is computed
#' with respect to other locations in the intervention arm ([Anaya-Izquierdo & Alexander(2020)](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.13316)). The half-depth within
#' the half-depth within the control cloud (dc) is also computed. \code{CRTspat} returns the proportion di/(dc + di). \cr\cr
#' If an auxiliary \code{auxiliary} \code{"CRTsp"} object is specified then either \code{distanceAssigned} or \code{nearestDiscord} (the default)
#' is computed with respect to the assignments in the auxiliary. If the auxiliary is a grid with \code{design$geometry} set to \code{'triangle'},
#' \code{'square'} or \code{'hexagon'} then the distance is computed to the edge of the nearest grid pixel in the discordant arm
#' (using a circular approximation for the perimeter) rather than to the point location itself.\cr
#' @export
#' @examples{
#' # Calculate the disc with a radius of 0.5 km
#' exampletrial <- compute_distance(trial = readdata('exampleCRT.txt'),
#' distance = 'disc', scale_par = 0.5)
#' }
compute_distance <- function(trial, distance = "nearestDiscord", scale_par = NULL, auxiliary = NULL) {
  CRT <- CRTsp(trial)
  trial <- CRT$trial
  require_nearestDiscord <- is.null(trial$nearestDiscord) & identical(distance, "nearestDiscord")
  require_hdep <- is.null(trial$hdep) & identical(distance, "hdep")
  require_sdep <- is.null(trial$sdep) & identical(distance, "sdep")
  require_disc <- identical(distance, "disc") &
                (is.null(trial$disc) | !identical(CRT$design$disc$scale_par, scale_par))
  require_kern <- identical(distance, "kern") &
                (is.null(trial$kern) | !identical(CRT$design$kern$scale_par, scale_par))
  kern <- NULL
  if (!is.null(auxiliary)) {
    if (!(distance %in% c("nearestDiscord", "distanceAssigned"))) {
      stop('*** Distance measure is incompatible with use of auxiliary data object ***')
    }
    # function calls for computing distance to auxiliary
    dist_vec <- interrogate_auxiliary(CRT = CRT, auxiliary = auxiliary, distance = distance)
    if (all(is.na(dist_vec$dist_corrected))) {
      warning("*** arm or cluster assignments not available for calculation of distance ***")
    } else {
      if (identical(distance, "nearestDiscord")) {
        trial$nearestDiscord <- ifelse(trial$arm == "control",
                                       -dist_vec$dist_corrected, dist_vec$dist_corrected)
        message("*** nearestDiscord is distance to nearest pixel in discordant arm ***")
      } else {
        trial$distanceAssigned <- dist_vec$dist_corrected
        message("*** distanceAssigned is distance to nearest pixel in the assigned cluster ***")
      }
    }
  } else {
    require_nearestDiscord <- is.null(trial$nearestDiscord) & identical(distance, "nearestDiscord")
    if (is.null(trial[[distance]] & is.null(trial$arm)))
      stop('*** Randomization is required for computation of distances or surrounds ***')
    if (require_hdep | require_sdep){
      depthilist <- apply(trial, MARGIN = 1, FUN = depths, trial = trial, cloud = 'intervention')
      depthi_df <- as.data.frame(do.call(rbind, lapply(depthilist, as.data.frame)))
      depthclist <- apply(trial, MARGIN = 1, FUN = depths, trial = trial, cloud = 'control' )
      depthc_df <- as.data.frame(do.call(rbind, lapply(depthclist, as.data.frame)))
      trial$hdep <- depthi_df$hdep/(depthc_df$hdep + depthi_df$hdep)
      trial$sdep <- depthi_df$sdep/(depthc_df$sdep + depthi_df$sdep)

      # replace NA with limiting value, depending which arm it is in (these points are on the outside of the cloud)
      trial$hdep[is.na(trial$hdep)] <- ifelse(trial$arm[is.na(trial$hdep)] == 'intervention', 1, 0)
      trial$sdep[is.na(trial$sdep)] <- ifelse(trial$arm[is.na(trial$sdep)] == 'intervention', 1, 0)

      CRT$design$hdep <- distance_stats(trial, distance = "hdep")
      CRT$design$sdep <- distance_stats(trial, distance = "sdep")
    }

    if ((require_nearestDiscord | require_disc | require_kern)){
        dist_trial <- as.matrix(dist(cbind(trial$x, trial$y), method = "euclidean"))
        if (require_nearestDiscord){
            discord <- outer(trial$arm, trial$arm, "!=")  #true & false.
            discord_dist_trial <- ifelse(discord, dist_trial, Inf)
            trial$nearestDiscord <- ifelse(trial$arm == "control", -apply(discord_dist_trial,
                       MARGIN = 2, min), apply(discord_dist_trial, MARGIN = 2, min))
            CRT$design$nearestDiscord <- distance_stats(trial, distance = "nearestDiscord")
            message("*** computed distance to nearest measurements in discordant arm ***")
        }
        if (require_disc){
            if (is.null(scale_par)) {
              stop("*** radius (scale_par) must be specified for computation of disc ***")
            }
            neighbours <- colSums(dist_trial <= scale_par)
            intervened_neighbours <- colSums(trial$arm =='intervention' & (dist_trial <= scale_par))
            trial$disc <- intervened_neighbours/neighbours
            CRT$design$disc <- distance_stats(trial, distance = "disc")
            CRT$design$disc$scale_par <- scale_par
        }
        if (require_kern){
            if (is.null(scale_par)) {
              stop("*** s.d. (scale_par) must be specified for computation of kern ***")
            }
            weighted_neighbours <- colSums(dnorm(dist_trial, mean = 0, sd = scale_par))
            weighted_intervened <- colSums(dnorm(dist_trial, mean = 0, sd = scale_par) * matrix(data = (trial$arm == 'intervention'),
                         nrow = nrow(trial), ncol = nrow(trial)))
            trial$kern <- weighted_intervened/weighted_neighbours
            CRT$design$kern <- distance_stats(trial, distance = "kern")
            CRT$design$kern$scale_par <- scale_par
        }
    }
  }
  CRT$trial <- trial
  return(CRT)
}


depths <- function(X, trial, cloud) {
  # this is an R translation of the fortran code in
  # Rousseeuw & Ruts https://www.jstor.org/stable/2986073
  # algorithm as 307.1 Appl.Statist. (1996), vol.45, no.4
  # calculation of the simplicial depth and
  # the half space depth
  # u and v are the coordinates of the arbitrary point
  u <- as.numeric(X[["x"]])
  v <- as.numeric(X[["y"]])

  # for the CRT application, depth is computed with respect to the set of intervention individuals
  # excluding the point itself (if it is in the intervention arm)
  trial <- trial[trial$arm == cloud & (trial$x != u | trial$y != v),]

  n <- nrow(trial)
  x <- trial$x
  y <- trial$y
  nums <- 0
  numh <- 0
  sdep <- 0  # simplicial depth
  hdep <- 0  # half-space depth
  eps <- 1e-06
  nt <- 0
  # construct the vector alpha
  alpha <- fval <- rep(NA, nrow(trial))
  for (i in 1:n) {
    d <- sqrt((x[i] - u) * (x[i] - u) + (y[i] - v) * (y[i] - v))
    if (d <= eps) {
      nt <- nt + 1
    } else {
      xu <- (x[i] - u)/d
      yu <- (y[i] - v)/d
      if (abs(xu) > abs(yu)) {
        if (x[i] >= u) {
          alpha[i - nt] <- asin(yu)
          if (alpha[i - nt] < 0.0) {
            alpha[i - nt] <- 2 * pi + alpha[i - nt]
          }
        } else {
          alpha[i - nt] <- pi - asin(yu)
        }
      } else {
        if (y[i] >= v) {
          alpha[i - nt] <- acos(xu)
        } else {
          alpha[i - nt] <- 2 * pi - acos(xu)
        }
      }
      if (alpha[i - nt] >= (2 * pi - eps)) alpha[i - nt] <- 0.0
    }
  }
  nn <- n - nt
  if (nn > 1) {
    # nn is the number of elements of alpha that have been assigned a value
    # the missing elements should be removed
    #call sort (alpha, nn)
    alpha <- alpha[!is.na(alpha)]
    alpha <- alpha[order(alpha)]
    # check whether theta=(u,v) lies outside the data cloud
    angle <- alpha[1] - alpha[nn] + 2 * pi
    for (i in 2:nn) {
      angle <- max(angle, (alpha[i] - alpha[i - 1]))
    }
    if (angle <= (pi + eps)) {
      # make smallest alpha equal to zero, and compute nu = number of alpha < pi
      angle <- alpha[1]
      nu <- 0
      for (i in 1:nn) {
        alpha[i] <- alpha[i] - angle
        if (alpha[i] < (pi - eps)) nu <- nu + 1
      }
      if (nu < nn) {
        # merge sort the alpha with their antipodal angles beta, and at the same time
        # update i,fval[i], and nbad
        ja <- 1
        jb <- 1
        alphk <- alpha[1]
        betak <- alpha[nu + 1] - pi
        nn2 <- nn * 2
        nbad <- 0
        i <- nu
        nf <- nn
        for (j in 1:nn2) {
          if ((alphk + eps) < betak) {
            nf <- nf + 1
            if (ja < nn) {
              ja <- ja + 1
              alphk <- alpha[ja]
            } else {
              alphk <- 2 * pi + 1
            }
          } else {
            i <- i + 1
            if (identical(i,(nn + 1))) {
              i <- 1
              nf <- nf - nn
            }
            fval[i] <- nf
            nbad <- nbad + k((nf - i), 2)
            if (jb < nn) {
              jb <- jb + 1
              if ((jb + nu) <= nn) {
                betak <- alpha[jb + nu] - pi
              } else {
                betak <- alpha[jb + nu - nn] + pi
              }
            } else {
              betak <- 2 * pi + 1.0
            }
          }
        }
        nums <- k(nn, 3) - nbad
        # computation of numh for half space depth
        gi <- 0
        ja <- 1
        angle <- alpha[1]
        numh <- min(fval[1], (nn - fval[1]))
        for (i in 2:nn) {
          if (alpha[i] <= (angle + eps)) {
            ja <- ja + 1
          } else {
            gi <- gi + ja
            ja <- 1
            angle <- alpha[i]
          }
          ki <- fval[i] - gi
          numh <- min(numh, min(ki, (nn - ki)))
        }
        # adjust for the number nt of datapoints equal to theta
      }
    }
  }
  nums <- nums + k(nt, 1) * k(nn, 2) + k(nt, 2) * k(nn, 1) + k(nt, 3)
  if (n >= 3) sdep <- nums/k(n, 3)
  numh <- numh + nt
  hdep <- numh/n
  depths <- list(numh = numh, hdep = hdep, sdep = sdep)
  return(depths)
}


k <- function(m, j) {
  # algorithm as 307.2 appl.statist. (1996),vol.45, no.4
  # returns the value zero if m <j; otherwise
  # computes the number of combinations of j out of m

  if (m < j) {
    k <- 0
  } else {
    if (j == 1)
      k <- m
    if (j == 2)
      k <- (m * (m - 1))/2
    if (j == 3)
      k <- (m * (m - 1) * (m - 2))/6
  }
  return(k)
}

# This could be incorporated into calculate_distance
distance_stats <- function(trial, distance){
  trial$distance <- trial[[distance]]
  formula <- stats::as.formula("distance ~ cluster")
  aov <- summary(aov(data = trial, formula = formula))
  within_cluster_sd <- sqrt(aov[[1]]$`Mean Sq`[2])
  rSq <- aov[[1]]$`Sum Sq`[1]/(aov[[1]]$`Sum Sq`[1] + aov[[1]]$`Sum Sq`[2])
  distance_stats <- c(as.list(summary(trial[[distance]])),
                      list(sd = sd(trial[[distance]]),
                           within_cluster_sd = within_cluster_sd, rSq = rSq))
  return(distance_stats)
}

# Functions for trials with cluster assignments in an auxiliary list

interrogate_auxiliary <- function(CRT, auxiliary, distance){
  # first compute the offset for converting the coordinates onto the same (x,y) scale
  scalef <- 180/(6371*pi)
  offset <- list(y = (CRT$geom_full$centroid$lat - auxiliary$geom_full$centroid$lat)/scalef,
                 x = ((CRT$geom_full$centroid$long - auxiliary$geom_full$centroid$long) *
                   cos(CRT$geom_full$centroid$lat * pi/180))/scalef)

  # dist_df is a data frame containing the relevant distance measure
  dist_list <- apply(CRT$trial, FUN = dist_auxiliary, MARGIN = 1,
                            auxiliary = auxiliary, offset = offset, distance = distance)
  dist_df <- as.data.frame(do.call(rbind, lapply(dist_list, as.data.frame)))

  # subtract the average pixel radius using a circular approximation (unless the auxiliary contains points)
  if(identical(auxiliary$design$geometry, 'point')){
    dist_df$dist_corrected <- dist_df$dist_uncorrected
  } else {
    # first compute the minimum distance between grid points
    d <- min(dist(cbind(auxiliary$trial$x, auxiliary$trial$y), method = "euclidean"))
    dist_df$dist_corrected <- dist_df$dist_uncorrected - switch(auxiliary$design$geometry,
                                                'triangle' = d * sqrt(3 * sqrt(3)/(4 * pi)),
                                                'square' = d * sqrt(1/ pi),
                                                'hexagon' = d * sqrt(sqrt(3)/(2 * pi)))
  }
  return(dist_df)
}

dist_auxiliary <- function(X, auxiliary, offset, distance){
# minimum distance to discordant element in auxiliary CRTsp object
# the data coordinates are recentered to match the centroid to that of the auxiliary
  x <- as.numeric(X[['x']]) + offset[['x']]
  y <- as.numeric(X[['y']]) + offset[['y']]
  cluster <- arm <- NA
  if (identical(distance, "nearestDiscord")) {
    if(!is.null(auxiliary$trial$arm) & 'arm' %in% names(X)) {
      d2 <- (x-auxiliary$trial$x[auxiliary$trial$arm != X[['arm']]])^2 +
            (y-auxiliary$trial$y[auxiliary$trial$arm != X[['arm']]])^2
    } else {
      d2 <- NA
    }
  } else if (identical(distance, "distanceAssigned")) {
    if(!is.null(auxiliary$trial$cluster) & 'cluster' %in% names(X)) {
        d2 <- (x-auxiliary$trial$x[auxiliary$trial$cluster == X[['cluster']]])^2 +
          (y-auxiliary$trial$y[auxiliary$trial$cluster == X[['cluster']]])^2
    } else {
        d2 <- NA
    }
  } else if (identical(distance, "nearestPixel")) {
    d2 <- (x-auxiliary$trial$x)^2 +
          (y-auxiliary$trial$y)^2
    equal_to_min <- logical(length(d2))
    equal_to_min[match(min(d2), d2)] <- TRUE
    recd <- which(equal_to_min)[1]
    if(!is.null(auxiliary$trial$cluster)) cluster <- as.character(auxiliary$trial$cluster[recd])
    if(!is.null(auxiliary$trial$arm)) arm <- as.character(auxiliary$trial$arm[recd])
  }
return(list(dist_uncorrected = sqrt(min(d2)),
            cluster = cluster, arm = arm))}



