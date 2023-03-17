#' Graphical display of the geography of a CRT
#'
#' \code{plot.CRT} returns a stacked bar chart of the
#' the outcome grouped by distance from the arm boundary
#' @param CRT an object of class \code{'CRT'}; \cr
#' @param ... other arguments of \code{base::plot}
#' @param map logical: indicator of whether a map is required
#' @param fill fill layer of map
#' \tabular{ll}{
#' \code{'cluster'} \tab cluster assignment \cr
#' \code{'arms'}   \tab arm assignment \cr
#' \code{'none'}\tab No fill \cr
#' }
#' @param num name of numerator variable
#' @param denom name of denominator variable
#' @param showLocations logical: determining whether locations are shown
#' @param showClusterBoundaries logical: determining whether cluster boundaries are shown
#' @param showClusterLabels logical: determining whether the cluster numbers are shown
#' @param cpalette colour palette (to use different colours for clusters must be at least as
#' long as the number of clusters, defaults to rainbow())
#' @param maskbuffer radius of buffer drawn around inhabited areas (km)
#' @param labelsize size of cluster number labels
#' @param legend.position (using ggplot2::themes syntax)
#' @return graphics object produced by the ggplot2 package
#' @details If \code{map = TRUE} a map is produced
#' If \code{map = FALSE} a stacked bar plot of the data is produced
#' @importFrom ggplot2 aes alpha
#' @export
#' @examples
#' #Plot locations only
#' plot(CRT = readdata('test.CRT.csv'), map = TRUE, fill = 'none', showLocations = TRUE,
#'            showClusterBoundaries=FALSE, maskbuffer=0.2)
#'
#' #Show cluster boundaries and number clusters
#' plot(CRT = readdata('test.CRT.csv'), map = TRUE, fill ='none', showClusterBoundaries=TRUE,
#'            showClusterLabels=TRUE, maskbuffer=0.2)
#'
#' #Plot clusters in colour
#' plot(CRT=readdata('test.CRT.csv'), map = TRUE, fill = 'clusters', showClusterLabels = TRUE,
#'           labelsize=4, maskbuffer=0.2)
#'
#' #Plot arms
#' plot(CRT=readdata('test.CRT.csv'), maskbuffer=0.2, legend.position=c(0.2,0.8))
#' @examples
#' CRT=readdata('test_Simulate_CRT.csv'); CRT <- as_CRT(trial); plot(CRT)
plot.CRT <- function(CRT, ..., map = FALSE, fill = "arms", num = num, denom = denom,
    showLocations = FALSE, showClusterBoundaries = TRUE, showClusterLabels = FALSE,
    cpalette = NULL, maskbuffer = 0.2, labelsize = 4,
    legend.position = NULL) {
    if (!map) {
        if (is.null(cpalette)) cpalette <- c("#D55E00", "#0072A7")
        outcome <- positives <- negatives <- frequency <- dcat <- NULL
        an <- analyseCRT(trial = CRT$trial, method = "EMP")
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
    } else {

        trial <- CRT$trial

        arm <- cluster <- x <- y <- NULL
        colourClusters <- identical(fill, "clusters")
        showArms <- identical(fill, "arms")

        # The plotting routines require unique locations
        CRT <- as_CRT(trial)
        CRT <- aggregate(CRT)

        # The plotting routines use (x,y) coordinates
        if (is.null(CRT$trial$x)) {
          CRT <- latlong_as_xy(CRT)
        }

        trial <- CRT$trial
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

        clusters <- vor %>%
          sf::st_join(tr, sf::st_intersects) %>%
          group_by(cluster) %>%
          dplyr::summarize()

        totalClusters <- length(unique(trial$cluster))

        if (is.null(cpalette))
          cpalette <- sample(rainbow(totalClusters))
        if (totalClusters == 1)
          cpalette <- c("white")

        g <- ggplot2::ggplot()
        if (showClusterBoundaries) {
          g <- g + ggplot2::geom_sf(data = clusters, color = "black", fill = "white")
        }
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
        g <- add_annotations(trial = trial, showLocations = showLocations,
            showClusterLabels = showClusterLabels, maskbuffer = maskbuffer,
            labelsize = labelsize, legend.position = legend.position, g = g)
    }
    return(g)
}



add_annotations <- function(trial, showLocations, showClusterLabels, maskbuffer,
    labelsize, legend.position, g) {
    if (showClusterLabels) {
        showLocations <- FALSE
    }
    if (showLocations) {
        g <- g + ggplot2::geom_point(data = trial, aes(x = x, y = y), size = 0.5)
    }

    # Positions of centroids of clusters for locating the labels
    cc <- data.frame(trial %>%
        group_by(cluster) %>%
        dplyr::summarize(x = mean(x), y = mean(y), .groups = "drop"))
    if (showClusterLabels) {
        g <- g + ggplot2::geom_text(data = cc, aes(x = x, y = y, label = cluster),
            hjust = 0.5, vjust = 0.5, size = labelsize)
    }

    # mask for remote areas the plot limits are determined by the min
    # and max of the coordinates
    xlim <- c(min(trial$x - maskbuffer), max(trial$x + maskbuffer))
    ylim <- c(min(trial$y - maskbuffer), max(trial$y + maskbuffer))

    # the mask needs to extend outside the plot area
    x0 <- xlim[1] - 0.5
    x1 <- xlim[2] + 0.5
    y0 <- ylim[1] - 0.5
    y1 <- ylim[2] + 0.5
    bbox <- sf::st_polygon(list(cbind(x = c(x0, x1, x1, x0, x0), y = c(y0,
        y0, y1, y1, y0))))
    bbox <- sf::st_sfc(bbox)

    # mask for uninhabited areas
    tr <- sf::st_as_sf(trial, coords = c("x", "y"))
    buf1 <- sf::st_buffer(tr, maskbuffer)
    buf2 <- sf::st_union(buf1)
    mask <- sf::st_difference(bbox, buf2)

    g <- g + ggplot2::geom_sf(data = mask, fill = "grey")
    g <- g + ggplot2::theme(legend.position = legend.position)
    g <- g + ggplot2::theme(panel.border = ggplot2::element_blank())
    g <- g + ggplot2::theme(axis.title = ggplot2::element_blank())
    g <- g + ggplot2::coord_sf(expand = FALSE, xlim = xlim, ylim = ylim)

    return(g)
}

#' Plot results of statistical analysis of a CRT
#'
#' \code{plot.CRTanalysis} creates graphics showing the results of statistical analysis of a CRT
#' @param analysis object of S3 class \code{'CRTanalysis'} produced by \code{analyseCRT()}
#' @param ... other arguments of \code{base::plot}
#' @param map logical: indicator of whether a map is required
#' @param fill fill layer of map
#' \tabular{ll}{
#' \code{'distance'} \tab distance to the nearest discordant location\cr
#' \code{'prediction'}\tab model prediction of the outcome \cr
#' \code{'none'}\tab No fill \cr
#' }
#' @param showLocations logical: determining whether locations are shown
#' @param showClusterLabels logical: determining whether the cluster numbers are shown
#' @param showContamination logical: determining whether the estimated contamination range should be shown
#' @param cpalette colour palette (to use different colours for clusters must be at
#' least as long as the number of clusters, defaults to rainbow())
#' @param maskbuffer radius of buffer drawn around inhabited areas (km)
#' @param labelsize size of cluster number labels
#' @param legend.position (using ggplot2::themes syntax)
#' @return graphics object produced by the ggplot2 package
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct group_by summarize
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 aes
#' @details
#' If \code{map = TRUE}
#' If \code{map = FALSE} a plot of the estimated contamination function is generated.\cr
#' The fitted contamination function is plotted as a continuous blue line against the distance
#' from the nearest discordant
#' location.Using the same axes, data summaries are plotted for ten categories of distance from the boundary. Both the
#' average of the outcome and confidence intervals are plotted. \cr
#' For analyses with logit link function the outcome is plotted as a proportion. \cr
#' For analyses with log link function the outcome is plotted on a scale of the Williams' mean
#' (mean of exp(log(x + 1))) - 1) \cr
#' @export
#' @examples
#' analysis <- readdata('test_Analyse_CRT.txt'); class(analysis) <- 'CRTanalysis'; plot(analysis)
plot.CRTanalysis <- function(analysis, ..., map = FALSE, fill = "arms", showLocations = FALSE,
    showClusterLabels = FALSE, showContamination = FALSE,
    cpalette = NULL, maskbuffer = 0.2, labelsize = 4, legend.position = NULL) {
    if (!map) {
        d <- average <- upper <- lower <- contaminationFunction <- NULL
        interval <- analysis$contamination$contamination.limits
        g <- ggplot2::ggplot(data = analysis$contamination$data, aes(x = d,
            y = average))
        g <- g + ggplot2::theme_bw()
        g <- g + ggplot2::geom_point(size = 2)
        g <- g + ggplot2::geom_errorbar(mapping = aes(x = d, ymin = upper,
            ymax = lower), linewidth = 0.5, width = 0.1)
        g <- g + ggplot2::geom_line(data = analysis$contamination$FittedCurve,
            aes(x = d, y = contaminationFunction), linewidth = 2, colour = "#0072A7")
        g <- g + ggplot2::geom_vline(xintercept = interval, linewidth = 1)
        g <- g + ggplot2::geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed")
        g <- g + ggplot2::geom_rect(aes(xmin = interval[1], xmax = interval[2],
            ymin = -Inf, ymax = Inf), fill = alpha("#2C77BF", 0.02))
        g <- g + ggplot2::xlab("Distance from boundary (km)")
        g <- g + ggplot2::ylab("Outcome")
    } else {
        pixel <- analysis$inla.mesh$pixel
        raster <- analysis$inla.mesh$prediction
        contamination <- analysis$contamination
        trial <- analysis$trial
        showDistance <- identical(fill, "distance")
        showPrediction <- identical(fill, "prediction")
        if (showPrediction) showdistance <- FALSE

        # raster images derived from inla analysis
        x <- y <- prediction <- nearestDiscord <- NULL

        g <- ggplot2::ggplot()
        if (showPrediction) {
            g <- g + ggplot2::geom_tile(data = raster, aes(x = x, y = y,
                fill = prediction), alpha = 0.5, width = pixel, height = pixel)
            g <- g + ggplot2::scale_fill_gradient(name = "Prediction", low = "blue",
                high = "orange")

        }
        if (showDistance) {
            g <- g + ggplot2::geom_tile(data = raster, aes(x = x, y = y,
                fill = nearestDiscord), alpha = 0.5, width = pixel, height = pixel)
            g <- g + ggplot2::scale_fill_gradient(name = "Distance", low = "blue",
                high = "orange")

        }
        if (showContamination) {
            # augment the prediction grid with a classifier of whether
            # the point is within the contamination interval
            range <- contamination$contamination.limits
            raster$contaminated <- ifelse(dplyr::between(raster$nearestDiscord,
                range[1], range[2]), TRUE, FALSE)
            g <- g + ggplot2::geom_raster(data = raster[raster$contaminated,
                ], aes(x = x, y = y), fill = "black", alpha = 0.2)
        }
        g <- add_annotations(trial = trial, showLocations = showLocations,
            showClusterLabels = showClusterLabels, maskbuffer = maskbuffer,
            labelsize = labelsize, legend.position = legend.position, g = g)
    }
    return(g)
}



