#' Map of CRT trial area
#'
#' \code{plotCRTmap} returns a graphics object created using ggplot2.
#' Cartesian (x,y) coordinates are used. Units are expected to be km.
#'
#' @param trial either:\cr
#' an object of class \code{"CRT"}; \cr
#' a data frame containing locations in (x,y) coordinates,
#' (optionally) cluster assignments (factor \code{cluster}), and arm assignments (factor \code{arm});\cr
#' an object of class \code{CRTanalysis} produced by an \code{"INLA"} analysis.
#' @param fill fill layer of map
#' \tabular{ll}{
#' \code{"cluster"} \tab cluster assignment \cr
#' \code{"arms"}   \tab arm assignment \cr
#' \code{"distance"} \tab distance to the nearest discordant location\cr
#' \code{"prediction"}\tab model prediction of the outcome \cr
#' \code{"none"}\tab No fill \cr
#' }
#' @param showLocations logical: determining whether locations are shown
#' @param showClusterBoundaries logical: determining whether cluster boundaries are shown
#' @param showClusterLabels logical: determining whether the cluster numbers are shown
#' @param showContamination logical: determining whether the estimated contamination range should be shown
#' @param cpalette colour palette (to use different colours for clusters must be at least as long as the number of clusters, defaults to rainbow())
#' @param maskbuffer radius of buffer drawn around inhabited areas (km)
#' @param labelsize size of cluster number labels
#' @param legend.position (using ggplot2::themes syntax)
#' @return graphics object produced by the ggplot2 package
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct group_by summarize
#' @importFrom ggplot2 geom_polygon
#' @export
#'
#' @examples
#' #Plot locations only
#' plotCRTmap(trial = readdata('test_CRT2.csv'), fill = "none", showLocations = TRUE,
#'            showClusterBoundaries=FALSE, maskbuffer=0.2)
#'
#' #Show cluster boundaries and number clusters
#' plotCRTmap(trial = readdata('test_CRT2.csv'), fill ="none", showClusterBoundaries=TRUE,
#'            showClusterLabels=TRUE, maskbuffer=0.2)
#'
#' #Plot clusters in colour
#' plotCRTmap(trial=readdata('test_CRT2.csv'), fill = "clusters", showClusterLabels = TRUE,
#'           labelsize=2, maskbuffer=0.2)
#'
#' #Plot arms
#' plotCRTmap(trial=readdata('test_CRT2.csv'), maskbuffer=0.2, legend.position=c(0.2,0.8))
#'
plotCRTmap <- function(trial = trial, fill = 'arms', showLocations = FALSE, showClusterBoundaries=TRUE,
                       showClusterLabels = FALSE,  showContamination = FALSE,
                       cpalette = NULL, maskbuffer = 0.2, labelsize = 4, legend.position = NULL){

  pixel <- arm <- cluster <- x <- y <- prediction <- nearestDiscord <- NULL

  colourClusters <- identical(fill,"clusters")
  showArms <- identical(fill,"arms")
  showDistance <- identical(fill,"distance")
  showPrediction <- identical(fill,"prediction")

  if(identical(class(trial), "CRTanalysis")) {
    pixel <- trial$inla.mesh$pixel
    raster <- trial$inla.mesh$prediction
    contamination <- trial$contamination
    trial <- trial$trial
  }
  trial <- CRT_as_data.frame(trial)

  # The plotting routines require unique locations
  trial <- aggregateCRT(trial = trial)

  # The plotting routines use (x,y) coordinates
  if (is.null(trial$x)) {
    trial <- latlong_as_xy(trial)
  }

  # remove any buffer zones
  class(trial) <- "data.frame"
  if (!is.null(trial$buffer)) {
    trial <- trial[!trial$buffer, ]
  }

  # Adjust the required plots to exclude those for which there is no data or combinations that are too
  # cluttered or overprinted

  if (is.null(pixel)) {
    showdistance <- FALSE
    showPrediction <- FALSE
  } else {
     if (showPrediction) {
       showdistance <- FALSE
       colourClusters <- FALSE
       showArms <- FALSE
     } else if (showDistance) {
        colourClusters <- FALSE
        showArms <- FALSE
     }
  }
  if (is.null(trial$cluster)) {
    trial$cluster <- 1
    showClusterBoundaries <- FALSE
    showClusterLabels <- FALSE
  }
  if (is.null(trial$arm)) {
    trial$arm <- 0
    showArms <- FALSE
  }
  if (!showClusterBoundaries) {
    showClusterLabels <- FALSE
  }
  if (showClusterLabels) {
    showLocations <- FALSE
  }

  # create pts
  pts <- tidyr::tibble(y = trial$y, x = trial$x) %>%
    sf::st_as_sf(coords = c('y', 'x')) %>%
    sf::st_set_crs("Euclidean")

  tr <- sf::st_as_sf(trial, coords = c("x","y"))
  # voronoi of pts
  vor <- sf::st_voronoi(sf::st_combine(tr))

  vor <- sf::st_collection_extract(vor, "POLYGON")
  vor <- sf::st_as_sf(vor)

  clusters <- vor %>%
    sf::st_join(tr, sf::st_intersects) %>%
    group_by(cluster) %>%
    dplyr::summarize()

  totalClusters <- length(unique(trial$cluster))

# the plot limits are determined by the min and max of the coordinates
  xlim <- c(min(trial$x - maskbuffer),max(trial$x + maskbuffer))
  ylim <- c(min(trial$y - maskbuffer),max(trial$y + maskbuffer))

# the mask needs to extend outside the plot area
  x0 <- xlim[1] - 0.5
  x1 <- xlim[2] + 0.5
  y0 <- ylim[1] - 0.5
  y1 <- ylim[2] + 0.5
  bbox = sf::st_polygon(
    list(
      cbind(
        x = c(x0,x1,x1,x0,x0),
        y = c(y0,y0,y1,y1,y0))
    )
  )
  bbox <- sf::st_sfc(bbox)

# mask for uninhabited areas
  buf1 <- sf::st_buffer(tr, maskbuffer)
  buf2 <- sf::st_union(buf1)
  mask <- sf::st_difference(bbox, buf2)

  # Positions of centroids of clusters for locating the labels
  cc <- data.frame(trial %>%
                     group_by(cluster) %>%
                     dplyr::summarize(x = mean(x), y = mean(y), .groups = "drop"))


  if (is.null(cpalette)) cpalette <- sample(rainbow(totalClusters))
  if (totalClusters == 1) cpalette <- c("white")

  g <- ggplot2::ggplot()
  if (showClusterBoundaries) {
    g <- g + ggplot2::geom_sf(data = clusters, color = "black", fill = "white")
  }
  if (colourClusters) {
    g <- g + ggplot2::geom_sf(data = clusters, aes(fill = cluster), fill = cpalette, alpha=0.8)
  }
  if (showArms) {
    arms <- vor %>%
      sf::st_join(tr, sf::st_intersects) %>%
      group_by(arm) %>%
      dplyr::summarize()
    g <- g + ggplot2::geom_sf(data = arms, aes(fill = arm))
      # use standard colour-blind compatible palette
    g <- g + ggplot2::scale_fill_manual(name = "Arms", values = c("#b2df8a", "#1f78b4"),
                  labels = c("Control", "Intervention"))

  }

  # raster images derived from inla analysis
  if (showPrediction) {
    g <- g + ggplot2::geom_tile(data = raster,
                      aes(x = x, y = y, fill = prediction),
                      alpha = 0.5, width = pixel, height = pixel)
    g <- g + ggplot2::scale_fill_gradient(name = "Prediction", low = "blue", high = "orange")

  }
  if (showDistance) {
    g <- g + ggplot2::geom_tile(data = raster,
                                aes(x = x, y = y, fill = nearestDiscord),
                                alpha = 0.5, width = pixel, height = pixel)
    g <- g + ggplot2::scale_fill_gradient(name = "Distance", low = "blue", high = "orange")

  }
  if (showContamination) {
    # augment the prediction grid with a classifier of whether the point is within the contamination interval
    range <- contamination$contamination.limits
    raster$contaminated <- ifelse(dplyr::between(raster$nearestDiscord,
                                                                       range[1], range[2]), TRUE, FALSE)
    g <- g + ggplot2::geom_raster(data = raster[raster$contaminated,
    ], aes(x = x, y = y), fill = "black", alpha = 0.2)

  }
  ####################################################################

  if (showLocations) {
    g <- g + ggplot2::geom_point(data = trial, aes(x = x, y = y), size = 0.5)
  }
  if (showClusterLabels) {
    g <- g + ggplot2::geom_text(data = cc, aes(x = x, y = y, label = cluster), hjust = 0.5, vjust = 0.5, size = labelsize)
  }
  # mask for remote areas
  g <- g + ggplot2::geom_sf(data = mask, fill = "grey")
  g <- g + ggplot2::theme(legend.position = legend.position)
  g <- g + ggplot2::theme(panel.border = ggplot2::element_blank())
  g <- g + ggplot2::theme(axis.title = ggplot2::element_blank())
  g <- g + ggplot2::coord_sf(expand = FALSE, xlim = xlim, ylim = ylim)
  return(g)
}


#' \code{plotContamination} returns a plot of an estimated contamination function
#' @param analysis list object of class CRT produced by \code{analyseCRT()}
#' @return graphics object produced by the ggplot2 package
#' @importFrom ggplot2 aes
#' @details The fitted contamination function is plotted as a continous blue line against the distance from the nearest discordant
#' location.Using the same axes, data summaries are plotted for ten categories of distance from the boundary. Both the
#' average of the outcome and confidence intervals are plotted. \cr
#' For analyses with logit link function the outcome is plotted as a proportion. \cr
#' For analyses with log link function the outcome is plotted on a scale of the Williams' mean
#' (mean of exp(log(x + 1))) - 1) \cr
#' @export
#'
#' @examples
#' plotContamination(analysis = readdata('test_Analyse_CRT.txt'))

plotContamination <- function(analysis) {
  d <- average <- upper <- lower <- contaminationFunction <- NULL
  interval <- analysis$contamination$contamination.limits
  g <- ggplot2::ggplot(data = analysis$contamination$data,aes(x = d, y = average))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::geom_point(size = 2)
  g <- g + ggplot2::geom_errorbar(mapping = aes(x = d, ymin = upper, ymax = lower),linewidth = 0.5, width = 0.1)
  g <- g + ggplot2::geom_line(data = analysis$contamination$FittedCurve,
                              aes(x = d, y = contaminationFunction), linewidth = 2, colour = "#0072A7")
  g <- g + ggplot2::geom_vline(xintercept = interval, linewidth = 1)
  g <- g + ggplot2::geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed")
  g <- g + ggplot2::geom_rect(aes(xmin = interval[1], xmax = interval[2], ymin = -Inf, ymax = Inf), fill = alpha("#2C77BF", 0.02))
  g <- g + ggplot2::xlab("Distance from boundary (km)")
  g <- g + ggplot2::ylab("Outcome")
  return(g)
}


#' \code{plotDataByDistance} returns a stacked bar chart of the grouped data
#' the outcome grouped by distance from the arm boundary
#' @param trial a dataframe containing locations (x,y), cluster assignments, and arm assignments
#' @param num name of numerator variable
#' @param denom name of denominator variable
#' @param cpalette colour palette (to use different colours for clusters must be at least as long as the number of clusters, defaults to rainbow())
#' @return graphics object produced by the ggplot2 package
#' @importFrom ggplot2 aes alpha
#' @export
#'
#' @examples
#' plotDataByDistance(trial=readdata('test_Simulate_CRT.csv'))

plotDataByDistance <- function(trial = trial, num = num, denom = denom,
                                cpalette = c("#D55E00", "#0072A7")) {
  outcome <- positives <- negatives <- frequency <- dcat <- NULL
  analysis <- analyseCRT(trial = trial, method = "EMP")
  analysis$contamination$data$negatives <- with(analysis$contamination$data, total - positives)
  analysis$contamination$data$dcat <- with(analysis$contamination,
                                           min(FittedCurve$d) + (data$cats - 0.5) * (max(FittedCurve$d) -
                                                                                       min(FittedCurve$d))/10)
  data <- tidyr::gather(analysis$contamination$data[, c("dcat", "negatives", "positives")], outcome, frequency, positives:negatives,
                        factor_key = TRUE)
  plot <- ggplot2::ggplot(data = data,
                          aes(x = dcat, y = frequency, fill = outcome)) +
    ggplot2::theme_bw() + ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = cpalette, labels = c("Positive", "Negative")) +
    ggplot2::geom_vline(xintercept = 0,                                                                                                                 linewidth = 1, linetype = "dashed") + ggplot2::xlab("Distance from boundary (km)") + ggplot2::ylab("Frequency") +
    ggplot2::theme(legend.position = "bottom")
  return(plot)
}
