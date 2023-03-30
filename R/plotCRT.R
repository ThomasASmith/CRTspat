

#' Graphical displays of the geography of a CRT
#'
#' \code{plotCRT} returns graphical displays of the geography of a CRT
#' \code{plotCRT} creates graphics displaying the geographical showing the results of statistical analysis of a CRT
#' @param object object of class \code{'CRTanalysis'} produced by \code{CRTanalysis()}
#' @param map logical: indicator of whether a map is required
#' @param fill fill layer of map
#' \tabular{ll}{
#' \code{'cluster'} \tab cluster assignment \cr
#' \code{'arms'}   \tab arm assignment \cr
#' \code{'distance'} \tab distance to the nearest discordant location\cr
#' \code{'prediction'}\tab model prediction of the outcome \cr
#' \code{'none'}\tab No fill \cr
#' }
#' @param showLocations logical: determining whether locations are shown
#' @param showClusterLabels logical: determining whether the cluster numbers are shown
#' @param showClusterBoundaries logical: determining whether cluster boundaries are shown
#' @param showContamination logical: determining whether the estimated contamination range should be mapped
#' @param cpalette colour palette (to use different colours for clusters must be at
#' least as long as the number of clusters, defaults to \code{rainbow())}
#' @param maskbuffer radius of buffer drawn around inhabited areas (km)
#' @param labelsize size of cluster number labels
#' @param legend.position (using \code{ggplot2::themes} syntax)
#' @return graphics object produced by the \code{ggplot2} package
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct group_by summarize
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 aes
#' @details
#' If \code{map = TRUE} a plot corresponding to the value of \code{fill} is generated.
#'  \itemize{
#' \item \code{fill = 'clusters'} or leads to thematic map showing the locations of the clusters
#' \item \code{fill = 'arms'} leads to a thematic map showing the geography of the randomization
#' \item \code{fill = 'distance'} leads to a raster plot of the distance to the nearest discordant location.
#' \item \code{fill = 'prediction'} leads to a raster plot of predictions from an \code{'INLA'} model.
#' }
#' If \code{showContamination = TRUE} the map is overlaid with a grey transparent layer showing showing which
#' areas are within the contamination zone estimated by an \code{'INLA'} model.\cr
#' If \code{map = FALSE} and the input is a trial data frame or a \code{CRTsp} object, containing a randomisation
#' to arms, a stacked bar chart of the the outcome grouped by distance from the boundary is produced.\cr
#' If \code{map = FALSE} and the input is a \code{CRTanalysis} object plot of the
#' estimated contamination function is generated.
#' The fitted contamination function is plotted as a continuous blue line against the distance
#' from the nearest discordant location.Using the same axes, data summaries are plotted for
#' ten categories of distance from the boundary. Both the
#' average of the outcome and confidence intervals are plotted. \cr
#' For analyses with logit link function the outcome is plotted as a proportion. \cr
#' For analyses with log link function the outcome is plotted on a scale of the Williams' mean
#' (mean of exp(log(x + 1))) - 1) \cr
#' @export
#' @importFrom ggplot2 aes alpha
#' @examples
#' #Plot of data by distance
#' plotCRT(readdata('exampleCRT.csv'))
#' #Map of locations only
#' plotCRT(readdata('exampleCRT.csv'), map = TRUE, fill = 'none', showLocations = TRUE,
#'            showClusterBoundaries=FALSE, maskbuffer=0.2)
#' #show cluster boundaries and number clusters
#' plotCRT(readdata('exampleCRT.csv'), map = TRUE, fill ='none', showClusterBoundaries=TRUE,
#'            showClusterLabels=TRUE, maskbuffer=0.2, labelsize = 2)
#' #show clusters in colour
#' plotCRT(readdata('exampleCRT.csv'), map = TRUE, fill = 'clusters', showClusterLabels = TRUE,
#'           labelsize=2, maskbuffer=0.2)
#' #show arms
#' plotCRT(readdata('exampleCRT.csv'), map = TRUE,
#'         fill = 'arms', maskbuffer=0.2, legend.position=c(0.2,0.8))
#' #contamination plot
#' analysis <- CRTanalysis(readdata('exampleCRT.csv')); plotCRT(analysis, map = FALSE)
#'
#' @export
plotCRT <- function(object, map = FALSE, fill = "arms", showLocations = FALSE,
    showClusterBoundaries = TRUE, showClusterLabels = FALSE, showContamination = FALSE,
    cpalette = NULL, maskbuffer = 0.2, labelsize = 4, legend.position = NULL) {
    g <- NULL
    if (!isa(object, what = 'CRTanalysis')) object <- CRTsp(object)
    trial <- object$trial
    if (is.null(trial)) {
        stop("*** No data points for plotting ***")
    }
    if (!map) {
        if (isa(object, what = 'CRTanalysis')) {
            # if the object is the output from analysisCRT
            analysis <- object
            if (is.null(analysis$contamination$contamination_limits))
                stop("*** No analysis of contamination available ***")
            d <- average <- upper <- lower <- contaminationFunction <- NULL
            interval <- analysis$contamination$contamination_limits
            g <- ggplot2::ggplot(data = analysis$contamination$data, ggplot2::aes(x = d,
                                                                           y = average))
            g <- g + ggplot2::theme_bw()
            g <- g + ggplot2::geom_point(size = 2)
            g <- g + ggplot2::geom_errorbar(mapping = ggplot2::aes(x = d, ymin = upper,
                                        ymax = lower), linewidth = 0.5, width = 0.1)
            g <- g + ggplot2::geom_line(data = analysis$contamination$FittedCurve,
                                        ggplot2::aes(x = d, y = contaminationFunction), linewidth = 2,
                                        colour = "#0072A7")
            g <- g + ggplot2::geom_vline(xintercept = interval, linewidth = 1)
            g <- g + ggplot2::geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed")
            g <- g + ggplot2::geom_rect(ggplot2::aes(xmin = interval[1], xmax = interval[2],
                                         ymin = -Inf, ymax = Inf), fill = alpha("#2C77BF", 0.02))
            g <- g + ggplot2::xlab("Distance from boundary (km)")
            g <- g + ggplot2::ylab("Outcome")
        } else {
            # Plot of frequency by distance to boundary
            if (is.null(cpalette))
                cpalette <- c("#D55E00", "#0072A7")
            outcome <- positives <- negatives <- frequency <- dcat <- NULL
            if (is.null(object$trial$num)) {
                return(plot(object$trial))
            }
            an <- CRTanalysis(trial = object$trial, method = "EMP")
            an$contamination$data$negatives <- with(an$contamination$data, total -
                                                        positives)
            an$contamination$data$dcat <- with(an$contamination, min(FittedCurve$d) +
                                 (data$cats - 0.5) * (max(FittedCurve$d) - min(FittedCurve$d))/10)
            data <- tidyr::gather(an$contamination$data[, c("dcat", "negatives",
                                "positives")], outcome, frequency, positives:negatives, factor_key = TRUE)
            g <- ggplot2::ggplot(data = data, aes(x = dcat, y = frequency, fill = outcome)) +
                ggplot2::theme_bw() + ggplot2::geom_bar(stat = "identity") +
                ggplot2::scale_fill_manual(values = cpalette, labels = c("Positive",
                                "Negative")) + ggplot2::geom_vline(xintercept = 0, linewidth = 1,
                                linetype = "dashed") + ggplot2::xlab("Distance from boundary (km)") +
                ggplot2::ylab("Frequency") + ggplot2::theme(legend.position = "bottom")
        }
    } else {
        if (isa(object, what = 'CRTanalysis')) {
            # raster map
            analysis <- object
            contamination <- analysis$contamination
            showDistance <- identical(fill, "distance")
            showPrediction <- identical(fill, "prediction")
            if (showPrediction)
                showdistance <- FALSE
            if(showDistance | showPrediction | showContamination) {
                # raster images derived from inla analysis
                x <- y <- prediction <- nearestDiscord <- NULL

                g <- ggplot2::ggplot() + ggplot2::theme(aspect.ratio = 1)
                if (!identical(analysis$options$method, "INLA")) {
                    stop("*** Raster plots only available for outputs from INLA analysis ***")
                } else {
                    pixel <- analysis$inla_mesh$pixel
                    raster <- analysis$inla_mesh$prediction
                    if (showPrediction) {
                        g <- g + ggplot2::geom_tile(data = raster, aes(x = x, y = y,
                                                                       fill = prediction), alpha = 0.5, width = pixel, height = pixel)
                        g <- g + ggplot2::scale_fill_gradient(name = "Prediction",
                                                              low = "blue", high = "orange")

                    }
                    if (showDistance) {
                        g <- g + ggplot2::geom_tile(data = raster, aes(x = x, y = y,
                                                                       fill = nearestDiscord), alpha = 0.5, width = pixel, height = pixel)
                        g <- g + ggplot2::scale_fill_gradient(name = "Distance",
                                                              low = "blue", high = "orange")

                    }
                    if (showContamination) {
                        # augment the prediction grid with a classifier of
                        # whether the point is within the contamination
                        # interval
                        range <- contamination$contamination_limits
                        raster$contaminated <- ifelse(dplyr::between(raster$nearestDiscord,
                                                         range[1], range[2]), TRUE, FALSE)
                        g <- g + ggplot2::geom_tile(data = raster[raster$contaminated, ],
                            aes(x = x, y = y), fill = "black", alpha = 0.2)
                    }
                }
            }
        }
        g <- vectorPlot(trial = trial, fill = fill, showLocations = showLocations,
            showClusterBoundaries = showClusterBoundaries, cpalette = cpalette,
            showClusterLabels = showClusterLabels, maskbuffer = maskbuffer,
            labelsize = labelsize, legend.position = legend.position, g = g)
    }
return(g)
}


vectorPlot <- function(trial, fill, showLocations, showClusterBoundaries,
    showClusterLabels, cpalette = NULL, maskbuffer, labelsize, legend.position,
    g = NULL) {
    arm <- cluster <- x <- y <- NULL
    colourClusters <- identical(fill, "clusters")
    showArms <- identical(fill, "arms")

    # The plotting routines require unique locations
    CRT <- aggregateCRT(trial)

    # The plotting routines use (x,y) coordinates
    if (is.null(CRT$trial$x)) {
        CRT <- latlong_as_xy(CRT)
    }

    # remove any buffer zones
    if (!is.null(trial$buffer)) {
        trial <- trial[!trial$buffer, ]
    }

    # Adjust the required plots to exclude those for which there is no
    # data or combinations that are too cluttered or overprinted

    if (is.null(trial$cluster)) {
        trial$cluster <- rep(1, nrow(trial))
        showClusterBoundaries <- FALSE
        showClusterLabels <- FALSE
        colourClusters <- FALSE
    }
    if (is.null(trial$arm)) {
        trial$arm <- 0
        showArms <- FALSE
    }
    if (!showClusterBoundaries) {
        showClusterLabels <- FALSE
    }

    # create pts
    pts <- tidyr::tibble(y = trial$y, x = trial$x) %>%
        sf::st_as_sf(coords = c("y", "x")) %>%
        sf::st_set_crs("Euclidean")

    tr <- sf::st_as_sf(trial, coords = c("x", "y"))
    # voronoi of pts
    vor <- sf::st_voronoi(sf::st_combine(tr))

    vor <- sf::st_collection_extract(vor, "POLYGON")
    vor <- sf::st_as_sf(vor)

    if (!is.null(trial$cluster)) {
        clusters <- vor %>%
            sf::st_join(tr, sf::st_intersects) %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarize()
    }

    totalClusters <- length(unique(trial$cluster))

    if (is.null(cpalette))
        cpalette <- sample(rainbow(totalClusters))
    if (totalClusters == 1)
        cpalette <- c("white")

    if (is.null(g))
        g <- ggplot2::ggplot() + ggplot2::theme(aspect.ratio = 1)

    if (colourClusters) {
        g <- g + ggplot2::geom_sf(data = clusters, aes(fill = cluster), fill = cpalette,
            alpha = 0.8)
    }
    if (showArms) {
        arms <- vor %>%
            sf::st_join(tr, sf::st_intersects) %>%
            group_by(arm) %>%
            dplyr::summarize()
        g <- g + ggplot2::geom_sf(data = arms, aes(fill = arm))
        # use standard colour-blind compatible palette
        g <- g + ggplot2::scale_fill_manual(name = "Arms", values = c("#b2df8a",
            "#1f78b4"), labels = c("Control", "Intervention"))
    }
    if (showClusterBoundaries) {
        g <- g + ggplot2::geom_sf(data = clusters, color = "black", fill = NA)
    }
    g <- add_annotations(trial = trial, showLocations = showLocations, showClusterLabels = showClusterLabels,
        maskbuffer = maskbuffer, labelsize = labelsize, legend.position = legend.position,
        g = g)
    return(g)
}



add_annotations <- function(trial, showLocations, showClusterLabels, maskbuffer,
    labelsize, legend.position, g) {
    cluster <- x <- y <- NULL
    # mask for remote areas the plot limits are determined by the min
    # and max of the coordinates
    xlim <- c(min(trial$x - maskbuffer), max(trial$x + maskbuffer))
    ylim <- c(min(trial$y - maskbuffer), max(trial$y + maskbuffer))

    # mask for excluded areas the mask needs to extend outside the plot
    # area
    x0 <- xlim[1] - 0.5
    x1 <- xlim[2] + 0.5
    y0 <- ylim[1] - 0.5
    y1 <- ylim[2] + 0.5
    bbox <- sf::st_polygon(list(cbind(x = c(x0, x1, x1, x0, x0), y = c(y0,
        y0, y1, y1, y0))))
    bbox <- sf::st_sfc(bbox)
    tr <- sf::st_as_sf(trial, coords = c("x", "y"))
    buf1 <- sf::st_buffer(tr, maskbuffer)
    buf2 <- sf::st_union(buf1)
    mask <- sf::st_difference(bbox, buf2)

    g <- g + ggplot2::geom_sf(data = mask, fill = "grey")

    # Labels
    if (showClusterLabels) {
        showLocations <- FALSE
        # Positions of centroids of clusters for locating the labels
        cc <- data.frame(trial %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarize(x = mean(x), y = mean(y), .groups = "drop"))
        g <- g + ggplot2::geom_text(data = cc, aes(x = x, y = y, label = cluster),
            hjust = 0.5, vjust = 0.5, size = labelsize)
    }
    if (showLocations) {
        g <- g + ggplot2::geom_point(data = trial, aes(x = x, y = y), size = 0.5)
    }


    g <- g + ggplot2::theme(legend.position = legend.position)
    g <- g + ggplot2::theme(panel.border = ggplot2::element_blank())
    g <- g + ggplot2::theme(axis.title = ggplot2::element_blank())
    g <- g + ggplot2::coord_sf(expand = FALSE, xlim = xlim, ylim = ylim)

    return(g)
}


