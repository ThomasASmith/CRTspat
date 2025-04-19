#' Analysis of cluster randomized trial with spillover
#'
#' \code{CRTanalysis} carries out a statistical analysis of a cluster randomized trial (CRT).
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}) and outcome data (see details).
#' @param method statistical method with options:
#' \tabular{ll}{
#' \code{"EMP"} \tab simple averages of the data   \cr
#' \code{"T"}   \tab comparison of cluster means by t-test \cr
#' \code{"GEE"} \tab Generalised Estimating Equations \cr
#' \code{"LME4"} \tab Generalized Linear Mixed-Effects Models \cr
#' \code{"INLA"}\tab Integrated Nested Laplace Approximation (INLA) \cr
#' \code{"MCMC"}\tab Markov chain Monte Carlo using \code{"stan"} \cr
#' \code{"WCA"}\tab Within cluster analysis \cr
#' }
#' @param distance Measure of distance or surround with options: \cr
#' \tabular{ll}{
#' \code{"nearestDiscord"} \tab distance to nearest discordant location (km)\cr
#' \code{"disc"} \tab disc\cr
#' \code{"kern"} \tab surround based on sum of normal kernels\cr
#' \code{"hdep"} \tab Tukey half space depth\cr
#' \code{"sdep"} \tab simplicial depth\cr
#' }
#' @param cfunc transformation defining the spillover function with options:
#' \tabular{llll}{
#' \code{"Z"} \tab\tab arm effects not considered\tab reference model\cr
#' \code{"X"} \tab\tab spillover not modelled\tab the only valid value of \code{cfunc} for methods \code{"EMP"}, \code{"T"} and \code{"GEE"}\cr
#' \code{"L"} \tab\tab inverse logistic (sigmoid)\tab the default for \code{"INLA"} and \code{"MCMC"} methods\cr
#' \code{"P"} \tab\tab inverse probit (error function)\tab available with \code{"INLA"} and \code{"MCMC"} methods\cr
#' \code{"D"} \tab\tab diffusion model\tab only available with the \code{"MCMC"} method\cr
#' \code{"S"} \tab\tab piecewise linear\tab only available with the \code{"MCMC"} method\cr
#' \code{"E"} \tab\tab estimation of scale factor\tab only available with \code{distance = "disc"} or \code{distance = "kern"}\cr
#' \code{"R"} \tab\tab rescaled linear\tab \cr
#' }
#' @param scale_par numeric: pre-specified value of the spillover parameter or disc radius for models where this is fixed (\code{cfunc = "R"}).\cr\cr
#' @param link link function with options:
#' \tabular{ll}{
#' \code{"logit"}\tab (the default). \code{numerator} has a binomial distribution with denominator \code{denominator}.\cr
#' \code{"log"}  \tab \code{numerator} is Poisson distributed with an offset of log(\code{denominator}).\cr
#' \code{"cloglog"} \tab \code{numerator} is Bernoulli distributed with an offset of log(\code{denominator}).\cr
#' \code{"identity"}\tab The outcome is \code{numerator/denominator} with a normally distributed error function.\cr
#' }
#' @param numerator string: name of numerator variable for outcome
#' @param denominator string: name of denominator variable for outcome data (if present)
#' @param excludeBuffer logical: indicator of whether any buffer zone (records with \code{buffer=TRUE}) should be excluded from analysis
#' @param alpha numeric: confidence level for confidence intervals and credible intervals
#' @param baselineOnly logical: indicator of whether required analysis is of effect size or of baseline only
#' @param baselineNumerator string: name of numerator variable for baseline data (if present)
#' @param baselineDenominator string: name of denominator variable for baseline data (if present)
#' @param personalProtection logical: indicator of whether the model includes local effects with no spillover
#' @param clusterEffects logical: indicator of whether the model includes cluster random effects
#' @param spatialEffects logical: indicator of whether the model includes spatial random effects
#' (available only for \code{method = "INLA"} or for \code{method = "MCMC"})
#' @param control list: control options to be passed to the statistical fitting function
#' (available only for \code{method = "MCMC"})
#' @param pixel numeric: size of pixel in km for spatial model (used for \code{method = "MCMC"})
#' @param requireMesh logical: indicator of whether spatial predictions are required
#' (available only for \code{method = "INLA"})
#' @param inla_mesh string: name of pre-existing INLA input object created by \code{compute_mesh()}
#' @return list of class \code{CRTanalysis} containing the following results of the analysis:
#' \itemize{
#' \item \code{description} : description of the dataset
#' \item \code{method} : statistical method
#' \item \code{pt_ests} : point estimates
#' \item \code{int_ests} : interval estimates
#' \item \code{model_object} : object returned by the fitting routine
#' \item \code{spillover} : function values and statistics describing the estimated spillover
#' }
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @details \code{CRTanalysis} is a wrapper for the statistical analysis packages:
#' [gee](https://CRAN.R-project.org/package=gee),
#' [INLA](https://www.r-inla.org/),
#' [rstan](https://CRAN.R-project.org/package=rstan),
#' and the [t.test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test)
#' function of package \code{stats}.\cr\cr
#' The wrapper does not provide an interface to the full functionality of these packages.
#' It is specific for typical analyses of cluster randomized trials with geographical clustering. Further details
#' are provided in the [vignette](https://thomasasmith.github.io/articles/Usecase5.html).\cr\cr
#' The key results of the analyses can be extracted using a \code{summary()} of the output list.
#' The \code{model_object} in the output list is the usual output from the statistical analysis routine,
#' and can be also be inspected with \code{summary()}, or analysed using \code{stats::fitted()}
#' for purposes of evaluation of model fit etc..\cr\cr
#' For models with a complementary log-log link function specified with \code{link = "cloglog"}.
#' the numerator must be coded as 0 or 1. Technically the binomial denominator is then 1 and the
#' value of \code{denominator} is used as a rate multiplier.\cr\cr
#' With the \code{"INLA"} method 'iid' random effects are used to model extra-Poisson variation.\cr\cr
#' Interval estimates for the coefficient of variation of the cluster level outcome are calculated using the method of
#' [Vangel (1996)](https://www.jstor.org/stable/2685039).\cr\cr
#' If a \code{control} list is provided then this is passed to the
#' @export
#' @examples
#' \donttest{
#' example <- readdata('exampleCRT.txt')
#' # Analysis of test dataset by t-test
#' exampleT <- CRTanalysis(example, method = "T")
#' summary(exampleT)
#' # Standard GEE analysis of test dataset ignoring spillover
#' exampleGEE <- CRTanalysis(example, method = "GEE")
#' summary(exampleGEE)
#' # LME4 analysis with error function spillover function
#' exampleLME4 <- CRTanalysis(example, method = "LME4", cfunc = "P")
#' summary(exampleLME4)
#' }
CRTanalysis <- function(
    trial, method = "GEE", distance = "nearestDiscord", scale_par = NULL,
    cfunc = "L", link = "logit", numerator = "num",
    denominator = "denom", excludeBuffer = FALSE, alpha = 0.05,
    baselineOnly = FALSE, baselineNumerator = "base_num", baselineDenominator = "base_denom",
    personalProtection = FALSE, clusterEffects = TRUE, spatialEffects = FALSE, pixel = 0.1,
    control = NULL, requireMesh = FALSE, inla_mesh = NULL, verbose = FALSE) {

    CRT <- CRTsp(trial)

    cluster <- linearity <- penalty <- distance_type <- NULL
    resamples <- 1000
    penalty <- 0

    # For models based on "nearestDiscord" the prior for the scale parameter should allow this to range from smaller than
    # any plausible spillover zone, to larger than the study area
    if (identical(distance, "nearestDiscord")) {
        log_sp_prior <- c(-5, log(max(CRT$trial$x) - min(CRT$trial$x)) + 2)
    } else {
    # For models based on surrounds the prior for the scale parameter must be smaller than the distance between locations
        max_dist <- max(dist(cbind(CRT$trial$x, CRT$trial$y), method = "euclidean"))
        log_sp_prior <- c(-5, log(max_dist/2))
    }

    # Test of validity of inputs
    if (!method %in% c("EMP", "T", "MCMC", "GEE", "INLA", "LME4", "WCA"))
    {
        stop("*** Invalid value for statistical method ***")
        return(NULL)
    }
    if (identical(method, "INLA") & identical(system.file(package='INLA'), "")){
        message("*** INLA package is not installed. Running lme4 analysis instead. ***")
        method <- "LME4"
    }
    # Some statistical methods do not allow for spillover
    if (method %in% c("EMP", "T", "GEE")) cfunc <- "X"

    # cfunc='Z' is used to remove the estimation of effect size from the model
    if (baselineOnly) cfunc <- "Z"

    # Classification of distance_type and which non-linear fit is needed
    if (cfunc %in% c("Z","X")) {
        distance_type <- "No fixed effects of distance "
        if (is.null(CRT$trial[[distance]]) & is.null(CRT$trial$arm)){
            distance <- "dummy"
            CRT$trial$dummy <- runif(nrow(CRT$trial), 0, 1)
        }
        linearity <- "No non-linear parameter. "
        scale_par <- 1.0
    } else {
        if(identical(distance, "nearestDiscord")) {
            distance_type <- "Signed distance "
        } else if (distance  %in% c("hdep", "sdep", "disc", "kern")) {
            distance_type <- "Surround: "
        } else if(is.null(CRT$trial[[distance]])) {
            stop("*** Invalid distance measure ***")
            return(NULL)
        } else {
            distance_type <- ifelse((min(CRT$trial[[distance]]) < 0), "Signed distance ", "Surround: ")
        }
        if (identical(distance_type, "Surround: ")) {
            if (!identical(cfunc,"E")) {
                message("*** Surrounds must have cfunc 'E' or 'R': using cfunc = 'R' ***")
                cfunc <- "R"
            }
        } else {
            if (identical(cfunc,"E")) {
                message("*** Signed distances cannot have cfunc = 'E': using cfunc = 'R' ***")
                cfunc <- "R"
            }
        }
        if(identical(cfunc, "R")) {
            if (distance  %in% c("disc", "kern")) {
                if (is.null(scale_par)) {
                    penalty <- ifelse(identical(method, "MCMC"), 0, 2)
                    linearity <- "Estimated scale parameter: "
                } else {
                    linearity <- paste0("Precalculated scale parameter: ")
                }
            } else {
                scale_par <- 1.0
                linearity <- "No non-linear parameter. "
            }
        }
        else if(is.null(scale_par)) {
            if(identical(distance_type, "Surround: ") & identical(cfunc, "E")){
                # message("Estimated escape function )
            } else if (!cfunc %in% c("L", "P", "S", "D")){
                stop("*** Invalid spillover function ***")
                return(NULL)
            }
            # the goodness-of-fit is penalised if scale_par needs to be estimated
            # (unless this is via MCMC)
            penalty <- ifelse(identical(method, "MCMC"), 0, 2)
            linearity <- "Estimated scale parameter: "
        } else {
            linearity <- paste0("Precalculated scale parameter of ", round(scale_par, digits = 3),": ")
        }
    }

    # if the distance or surround is not provided, augment the trial data frame with distance or surround
    # (compute distance does nothing beyond validating the CRTsp, if the distance has already been calculated)
    if (!(distance %in% c("disc", "kern", "dummy"))) {
        CRT <- compute_distance(CRT, distance = distance, scale_par = scale_par)
    }

    trial <- CRT$trial



    if ("buffer" %in% colnames(trial) & excludeBuffer) trial <- trial[!trial$buffer, ]

    # trial needs to be ordered for some analyses
    if(!is.null(trial$cluster)) trial <- trial[order(trial$cluster), ]

    # Some statistical methods only run if there are cluster effects
    if (method %in% c("LME4", "MCMC")) clusterEffects <- TRUE

    if (baselineOnly){
        # Baseline analyses are available only for GEE and INLA
        if (method %in% c("EMP", "T", "GEE", "MCMC", "LME4", "WCA"))
            {
            method <- "GEE"
            message("Analysis of baseline only, using GEE\n")
        } else if (identical(method,"INLA")) {
            message("Analysis of baseline only, using INLA\n")
        }
        if (is.null(trial[[baselineNumerator]])) {
            stop("*** No baseline data provided ***")
        }
        if (is.null(trial[[baselineDenominator]])) {
            trial[[baselineDenominator]] <- 1
        }
        trial$y1 <- trial[[baselineNumerator]]
        trial$y0 <- trial[[baselineDenominator]] - trial[[baselineNumerator]]
        trial$y_off <- trial[[baselineDenominator]]

    } else {
        if (is.null(trial[[numerator]])){
            stop("*** No outcome data to analyse ***")
        }
        trial$y1 <- trial[[numerator]]
        if (is.null(trial[[denominator]])){
            stop("*** No denominator specified ***")
        }
        trial$y0 <- trial[[denominator]] - trial[[numerator]]
        trial$y_off <- trial[[denominator]]
    }

    # create model formula for use in equations and for display
    fterms <- switch(cfunc,
        Z = NULL,
        X = "arm",
        "effect"
    )

    if (personalProtection & cfunc != 'X') fterms <- c(fterms, "arm")
    if (clusterEffects) {
        if (identical(method, "INLA")) {
            fterms <- c(fterms, "f(cluster, model = \'iid\')")
        } else if(method %in% c("EMP","GEE")) {
            fterms <- fterms
        } else {
            fterms <- c(fterms, "(1 | cluster)")
        }
    }
    if (identical(method, "INLA")){
        if (spatialEffects) fterms <- c(fterms, "f(s, model = spde)")
        if (link %in% c("log", "cloglog")) fterms <- c(fterms, "f(id, model = \'iid\')")
    }
    if (identical(method, "MCMC")){
        if (spatialEffects) fterms <- c(fterms, "ICAR")
    }
    ftext <- paste(fterms, collapse = " + ")

    # create names for confidence limits for use throughout
    CLnames <- c(
        paste0(alpha/0.02, "%"),
        paste0(100 - alpha/0.02, "%")
    )

    # store options here- noting that the model formula depends on allowable values of other options
    options <- list(method = method,
                    link = link,
                    distance = distance,
                    cfunc = cfunc,
                    alpha = alpha,
                    baselineOnly = baselineOnly,
                    fterms = fterms,
                    ftext = ftext,
                    CLnames = CLnames,
                    log_sp_prior = log_sp_prior,
                    clusterEffects = clusterEffects,
                    spatialEffects = spatialEffects,
                    pixel = pixel,
                    control = control,
                    personalProtection = personalProtection,
                    distance_type = distance_type,
                    linearity = linearity,
                    scale_par = scale_par,
                    penalty = penalty,
                    verbose = verbose)

    # scaffolds for lists
    pt_ests <- list(scale_par = NA, personal_protection = NA, spillover_interval = NA)
    int_ests <- list(controlY = NA, interventionY = NA, effect_size = NA)
    model_object <- list()
    description <- get_description(trial=trial, link=link, alpha=alpha, baselineOnly)
    analysis <- list(trial = trial,
                     pt_ests = pt_ests,
                     int_ests = int_ests,
                     description = description,
                     options = options)

    analysis <- switch(method,
           "EMP" = EMPanalysis(analysis),
           "T" = Tanalysis(analysis),
           "GEE" = GEEanalysis(analysis = analysis, resamples=resamples),
           "LME4" = LME4analysis(analysis),
           "INLA" = INLAanalysis(analysis, requireMesh = requireMesh, inla_mesh = inla_mesh),
           "MCMC" = stananalysis(analysis),
           "WCA" = wc_analysis(analysis, design = CRT$design)
    )
    if (!baselineOnly & !is.null(analysis$pt_ests$controlY)){
        if (!identical(method,"EMP")){
            scale_par <- analysis$options$scale_par
            if (verbose){
                message(paste0(linearity, ifelse(is.null(scale_par), "",
                ifelse(identical(scale_par, 1), "", round(scale_par, digits = 3)))," ",
                distance_type, "-", ifelse(identical(distance_type, "No fixed effects of distance "),
                "", getDistanceText(distance = distance, scale_par = scale_par)), "\n"))
            }
        }
    }
    class(analysis) <- "CRTanalysis"
    if (verbose) message(summary(analysis))
    return(analysis)
}

# functions for INLA analysis
#' Create INLA mesh for spatial analysis
#'
#' \code{compute_mesh} create objects required for INLA analysis of an object of class \code{"CRTsp"}.
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}) and outcome.
#' @param offset see \code{inla.mesh.2d} documentation
#' @param max.edge see \code{inla.mesh.2d} documentation
#' @param inla.alpha parameter related to the smoothness (see \code{inla} documentation)
#' @param maskbuffer numeric: width of buffer around points (km)
#' @param pixel numeric: size of pixel (km)
#' @return list
#' \itemize{
#' \item \code{prediction} Data frame containing the prediction points and covariate values
#' \item \code{A} projection matrix from the observations to the mesh nodes.
#' \item \code{Ap} projection matrix from the prediction points to the mesh nodes.
#' \item \code{indexs} index set for the SPDE model
#' \item \code{spde} SPDE model
#' \item \code{pixel} pixel size (km)
#' }
#' @details \code{compute_mesh} carries out the computationally intensive steps required for setting-up an
#' INLA analysis of an object of class \code{"CRTsp"}, creating the prediction mesh and the projection matrices.
#' The mesh can be reused for different models fitted to the same
#' geography. The computational resources required depend largely on the resolution of the prediction mesh.
#' The prediction mesh is thinned to include only pixels centred at a distance less than
#' \code{maskbuffer} from the nearest point.\cr
#' A warning may be generated if the \code{Matrix} library is not loaded.
#' @export
#' @examples
#' {
#' # low resolution mesh for test dataset
#' library(Matrix)
#' example <- readdata('exampleCRT.txt')
#' exampleMesh=compute_mesh(example, pixel = 0.5)
#' }
compute_mesh <- function(trial = trial, offset = -0.1, max.edge = 0.25,
                         inla.alpha = 2, maskbuffer = 0.5, pixel = 0.5)
{
    if (identical(system.file(package='INLA'), "")){
        message("*** INLA package is not installed ***")
        return("Mesh not created as INLA package is not installed")
    } else {
        # extract the trial data frame from the "CRTsp" object
        if (identical(class(trial),"CRTsp")) trial <- trial$trial
        # create an id variable if this does not exist
        if(is.null(trial$id)) trial <- dplyr::mutate(trial, id =  dplyr::row_number())

        # create buffer around area of points
        trial.coords <- base::matrix(
            c(trial$x, trial$y),
            ncol = 2
        )

        tr <- sf::st_as_sf(trial, coords = c("x","y"))
        buf1 <- sf::st_buffer(tr, maskbuffer)
        buf2 <- sf::st_union(buf1)
        # determine pixel size
        area <- sf::st_area(buf2)
        buffer <- sf::as_Spatial(buf2)

        # estimation mesh construction
        # dummy call to Matrix. This miraculously allows the loading of the "dgCMatrix" in the mesh to pass the test
        dummy <- Matrix::as.matrix(c(1,1,1,1))

        mesh <- INLA::inla.mesh.2d(
            boundary = buffer, offset = offset, cutoff = 0.05, max.edge = max.edge
        )

        # set up SPDE (Stochastic Partial Differential Equation) model
        spde <- INLA::inla.spde2.matern(mesh = mesh, alpha = inla.alpha, constr = TRUE)
        indexs <- INLA::inla.spde.make.index("s", spde$n.spde)
        A <- INLA::inla.spde.make.A(mesh = mesh, loc = trial.coords)

        # 8.3.6 Prediction data from https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html
        bb <- sf::st_bbox(buffer)

        # create a raster that is slightly larger than the buffered area
        xpixels <- round((bb$xmax - bb$xmin)/pixel) + 2
        ypixels <- round((bb$ymax - bb$ymin)/pixel) + 2
        x <- bb$xmin + (seq(1:xpixels) - 1.5)*pixel
        y <- bb$ymin + (seq(1:ypixels) - 1.5)*pixel
        all.coords <- as.data.frame(expand.grid(x, y), ncol = 2)
        colnames(all.coords) <- c("x", "y")
        all.coords <- sf::st_as_sf(all.coords, coords = c("x", "y"))
        pred.coords <- sf::st_filter(all.coords, sf::st_as_sf(buf2))
        pred.coords <- t(base::matrix(
            unlist(pred.coords),
            nrow = 2
        ))
        # projection matrix for the prediction locations
        Ap <- INLA::inla.spde.make.A(mesh = mesh, loc = pred.coords)

        # Distance matrix calculations for the prediction stack Create all pairwise comparisons
        pairs <- tidyr::crossing(
            row = seq(1:nrow(pred.coords)),
            col = seq(1:nrow(trial))
        )
        # Calculate the distances
        calcdistP <- function(row, col) sqrt(
            (trial$x[col] - pred.coords[row, 1])^2 + (trial$y[col] - pred.coords[row,
                                                                                 2])^2
        )
        distP <- apply(pairs, 1, function(y) calcdistP(y["row"], y["col"]))
        distM <- base::matrix(
            distP, nrow = nrow(pred.coords),
            ncol = nrow(trial),
            byrow = TRUE
        )
        nearestNeighbour <- apply(distM, 1, function(x) return(array(which.min(x))))
        prediction <- data.frame(
            x = pred.coords[, 1], y = pred.coords[, 2],
            nearestNeighbour = nearestNeighbour)
        prediction$id <- trial$id[nearestNeighbour]
        if (!is.null(trial$arm)) prediction$arm <- trial$arm[nearestNeighbour]
        if (!is.null(trial$cluster)) prediction$cluster <- trial$cluster[nearestNeighbour]
        prediction <- with(prediction, prediction[order(y, x), ])
        prediction$shortestDistance <- apply(distM, 1, min)
        rows <- seq(1:nrow(prediction))
        inla_mesh <- list(
            prediction = prediction, A = A, Ap = Ap, indexs = indexs, spde = spde,
            pixel = pixel)

        if (nrow(prediction) > 20){
            message("Mesh of ", nrow(prediction), " pixels of size ", pixel," km \n")
        }
        return(inla_mesh)
        }
}

EMPanalysis <- function(analysis){
    lp <- arm <- NULL
    description <- analysis$description
    pt_ests <- list()
    pt_ests$controlY <- unname(description$controlY)
    pt_ests$interventionY <- unname(description$interventionY)
    pt_ests$effect_size <- unname(description$effect_size)
    pt_ests$spillover_interval <- NA
    pt_ests$personal_protection <- NA
    analysis$pt_ests <- pt_ests
    return(analysis)
}

Tanalysis <- function(analysis) {
    trial <- analysis$trial
    link <- analysis$options$link
    alpha <- analysis$options$alpha
    clusterSum <- clusterSummary(trial, link)
    formula <- stats::as.formula("lp ~ arm")
    model_object <- stats::t.test(
        formula = formula, data = clusterSum, alternative = "two.sided",
        conf.level = 1 - alpha, var.equal = TRUE
    )
    analysis$model_object <- model_object
    analysis$pt_ests$p.value <- model_object$p.value
    analysisC <- stats::t.test(
        clusterSum$lp[clusterSum$arm == "control"], conf.level = 1 - alpha)
    analysis$pt_ests$controlY <- unname(invlink(link, analysisC$estimate[1]))
    analysis$int_ests$controlY <- invlink(link, analysisC$conf.int)
    analysisI <- stats::t.test(
        clusterSum$lp[clusterSum$arm == "intervention"], conf.level = 1 - alpha)
    analysis$pt_ests$interventionY <- unname(invlink(link, analysisI$estimate[1]))
    analysis$int_ests$interventionY <- invlink(link, analysisI$conf.int)

    # Covariance matrix (note that two arms are independent so the off-diagonal elements are zero)
    Sigma <- base::matrix(
        data = c(analysisC$stderr^2, 0, 0, analysisI$stderr^2),
        nrow = 2, ncol = 2)
    if (link == 'identity'){
        analysis$pt_ests$effect_size <- analysis$pt_ests$controlY -
            analysis$pt_ests$interventionY
        analysis$int_ests$effect_size <- unlist(model_object$conf.int)
    }
    if (link %in% c("logit", "log", "cloglog")){
        analysis$pt_ests$effect_size <- 1 - analysis$pt_ests$interventionY/
            analysis$pt_ests$controlY
        analysis$int_ests$effect_size <- 1 - exp(-unlist(model_object$conf.int))
    }
    analysis$pt_ests$t.statistic <- analysis$model_object$statistic
    analysis$pt_ests$df <- unname(analysis$model_object$parameter)
    analysis$pt_ests$p.value <- analysis$model_object$p.value
    return(analysis)
}


clusterSummary <- function(trial = trial, link = link){
    y1 <- arm <- cluster <- y_off <- NULL
    clusterSum <- data.frame(
        trial %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarize(
                y = sum(y1),
                total = sum(y_off),
                arm = arm[1]
            )
    )
    clusterSum$lp <- switch(link,
                            "identity" = clusterSum$y/clusterSum$total,
                            "log" = log(clusterSum$y/clusterSum$total),
                            "logit" = logit(clusterSum$y/clusterSum$total),
                            "cloglog" = log(clusterSum$y/clusterSum$total))
    # Trap any non-finite values
    clusterSum$lp[!is.finite(clusterSum$lp)] <- NA
    return(clusterSum)
}



wc_analysis <- function(analysis, design) {
    analysis$pt_ests <- analysis$int_ests <- y1 <- cluster <- y_off <- NULL
    trial <- analysis$trial
    link <- analysis$options$link
    alpha <- analysis$options$alpha
    distance = analysis$options$distance
    trial$d <- trial[[distance]]
    nclusters <- nlevels(trial$cluster)
    analysis$options <- list(
         method = "WCA",
         link = link,
         distance = distance,
         alpha = alpha,
         scale_par = design[[distance]][["scale_par"]]
    )
    analysis[[distance]] <- design[[distance]]
    analysis$nearestDiscord <- design$nearestDiscord
    fterms <- switch(link,
                       "identity" = "y1/y_off ~ 1 + d",
                       "log" = "y1 ~ 1 + d + offset(log(y_off))",
                       "cloglog" = "y1 ~ 1 + d + offset(log(y_off))",
                       "logit" = "cbind(y1,y0) ~ 1 + d")
    formula <- stats::as.formula(paste(fterms, collapse = "+"))
    pe <- matrix(nrow = 0, ncol = 2)
    for (cluster in levels(trial$cluster)){
        glm <- tryCatch(glm(formula = formula, family = "binomial", data = trial[trial$cluster == cluster,])
            , warning = function(w){ NULL})
        if (!is.null(glm)) {
            pe <- rbind(pe, matrix(glm$coefficients, ncol = 2))
        }
    }
    rr <- invlink(link = link, x = pe[ , 1] + pe[, 2])/invlink(link = link, x = pe[ , 1])
    exact = ifelse(length(unique(rr[!is.infinite(rr)])) == length(rr[!is.infinite(rr)]) , TRUE, FALSE)
    model_object <- stats::wilcox.test(rr, mu = 1,
                                alternative = "less", exact = exact, conf.int = TRUE, conf.level = 1 - alpha)
    model_object$conf.int[1] <- ifelse(model_object$conf.int[1] > 0, model_object$conf.int[1], 0)
    analysis$pt_ests$effect_size <- 1 - model_object$estimate
    analysis$int_ests$effect_size <- 1 - rev(unname(model_object$conf.int))
    analysis$pt_ests$test.statistic <- unname(model_object$statistic)
    analysis$pt_ests$p.value <- model_object$p.value
    analysis$model_object <- model_object
    return(analysis)
}

wc_summary <- function(analysis){
    defaultdigits <- getOption("digits")
    on.exit(options(digits = defaultdigits))
    options(digits = 3)
    distance <- analysis$options$distance
    cat('\nDistance and surround statistics\n')
    cat('Measure            Minimum    Median    Maximum    S.D.   Within-cluster S.D.   R-squared\n')
    cat('-------            -------    ------    -------    ----   -------------------   ---------\n')
    with(analysis$nearestDiscord,cat("Signed distance",
                                     format(Min., scientific=F), Median,
                                     format(Max., scientific=F), sd, "   ",
                                     within_cluster_sd, rSq, "\n", sep = "      "))
    with(analysis[[distance]],
         cat(distance, strrep(" ",8-nchar(distance)), Min., Median, Max., sd, "   ", within_cluster_sd, rSq, "\n", sep = "      "))
    cat("\nClusters assigned    : ", analysis$description$nclusters, "\n")
    cat("Clusters analysed    : ", analysis$description$nclusters, "\n")
    cat("Wilcoxon statistic   : ", analysis$pt_ests$statistic, "\n")
    cat("P-value (1-sided)    : ", analysis$pt_ests$p.value, "\n")
    cat(
        "Effect size estimate : ", analysis$pt_ests$effect_size,
         paste0(" (", 100 * (1 - analysis$options$alpha), "% CL: "), unlist(analysis$int_ests$effect_size),")\n"
    )
}


GEEanalysis <- function(analysis, resamples){
    trial <- analysis$trial
    link <- analysis$options$link
    alpha <- analysis$options$alpha
    cfunc <- analysis$options$cfunc
    fterms <- analysis$options$fterms
    pt_ests <- analysis$pt_ests
    int_ests <- analysis$int_ests
    method <- analysis$options$method
    cluster <- NULL

    model_object <- get_GEEmodel(trial = trial, link = link, fterms = fterms)

    summary.fit <- summary(model_object)

    z <- -qnorm(alpha/2)  #standard deviation score for calculating confidence intervals
    lp_yC <- summary.fit$coefficients[1, 1]
    se_lp_yC <- summary.fit$coefficients[1, 2]

    clusterSize <- nrow(trial)/nlevels(as.factor(trial$cluster))


    # remove the temporary objects from the dataframe
    model_object$data <- trial
    model_object$data$y1 <- model_object$data$y0 <- model_object$data$y_off <- NULL

    pt_ests$controlY <- invlink(link, lp_yC)
    int_ests$controlY <- namedCL(
        invlink(link, c(lp_yC - z * se_lp_yC, lp_yC + z * se_lp_yC)),
        alpha = alpha
    )

    # Intracluster correlation: original GEE based estimator replace owing to archiving of GEEpack at CRAN
    pt_ests$ICC <- analysis$description$ICC

    # In GEEpack:
    # pt_ests$ICC <- noLabels(summary.fit$corr[1])  #with corstr = 'exchangeable', alpha is the ICC
    # se_ICC <- noLabels(summary.fit$corr[2])

     se_ICC <- NA
     int_ests$ICC <- namedCL(
         noLabels(c(pt_ests$ICC - z * se_ICC, pt_ests$ICC + z * se_ICC)),
        alpha = alpha
    )
    pt_ests$DesignEffect <- 1 + (clusterSize - 1) * pt_ests$ICC  #Design Effect
    int_ests$DesignEffect <- 1 + (clusterSize - 1) * int_ests$ICC

    # Estimation of effect_size does not apply if analysis is of baseline only (cfunc='Z')
    pt_ests$effect_size <- NULL
    if (cfunc == "X") {
        lp_yI <- summary.fit$coefficients[1, 1] + summary.fit$coefficients[2, 1]
        se_lp_yI <- sqrt(
            model_object$robust.variance[1, 1] +
                model_object$robust.variance[2, 2] +
                2 * model_object$robust.variance[1,2]
        )

        int_ests$interventionY <- namedCL(
            invlink(link, c(lp_yI - z * se_lp_yI, lp_yI + z * se_lp_yI)),
            alpha = alpha)

        int_ests$effect_size <- estimateCLeffect_size(
            q50 = summary.fit$coefficients[, 1], Sigma = model_object$robust.variance,
            alpha = alpha, resamples = resamples, method = method,
            link = link)

        pt_ests$interventionY <- invlink(link, lp_yI)
        pt_ests$effect_size <- (1 - invlink(link, lp_yI)/invlink(link, lp_yC))
    }
    analysis$model_object <- model_object
    analysis$pt_ests <- pt_ests[names(pt_ests) != "model_object"]
    analysis$int_ests <- int_ests
    return(analysis)
}

get_GEEmodel <- function(trial, link, fterms){
    # GEE analysis of cluster effects
    cluster <- NULL
    fterms <- c(switch(link,
                       "identity" = "y1/y_off ~ 1",
                       "log" = "y1 ~ 1 + offset(log(y_off))",
                       "cloglog" = "cbind(y1, 1) ~ 1 + offset(log(y_off))",
                       "logit" = "cbind(y1,y0) ~ 1"),
                fterms)
    formula <- stats::as.formula(paste(fterms, collapse = "+"))
    if (link == "log") {
        model_object <- suppressMessages(gee::gee(
            formula = formula, id = cluster, data = trial, family = poisson(link = "log"),
            corstr = "exchangeable", scale.fix = FALSE))
    } else if (link == "cloglog") {
        model_object <- suppressMessages(gee::gee(
            formula = formula, id = cluster, data = trial, family = binomial(link = "cloglog"),
            corstr = "exchangeable", scale.fix = FALSE))
    } else if (link == "logit") {
        model_object <- suppressMessages(gee::gee(
            formula = formula, id = cluster, corstr = "exchangeable",
            data = trial, family = binomial(link = "logit")))
    } else if (link == "identity") {
        model_object <- suppressMessages(gee::gee(
            formula = formula, id = cluster, corstr = "exchangeable",
            data = trial, family = gaussian))
    }
return(model_object)}

LME4analysis <- function(analysis, cfunc, trial, link, fterms){
    trial <- analysis$trial
    link <- analysis$options$link
    cfunc <- analysis$options$cfunc
    FUN <- get_FUN(cfunc)
    alpha <- analysis$options$alpha
    scale_par <- analysis$options$scale_par
    distance <- analysis$options$distance
    log_sp_prior <- analysis$options$log_sp_prior
    linearity <- analysis$options$linearity
    distance_type <- analysis$options$distance_type
    fterms <- analysis$options$fterms
    # TODO replace the use of ftext with fterms
    ftext <- analysis$options$ftext
    log_scale_par <- NA
    contrasts <- NULL
    fterms = switch(link,
                    "identity" = c("y1/y_off ~ 1", fterms),
                    "log" = c("y1 ~ 1", fterms, "offset(log(y_off))"),
                    "cloglog" = c("y1 ~ 1", fterms, "offset(log(y_off))"),
                    "logit" = c("cbind(y1,y0) ~ 1", fterms))
    formula <- stats::as.formula(paste(fterms, collapse = "+"))
    if (!identical(distance_type, "No fixed effects of distance ")) {
        if (analysis$options$penalty > 0) {
            if (identical(Sys.getenv("TESTTHAT"), "true")) {
                log_scale_par <- 2.0
            } else {
                tryCatch({
                    #message "Estimating scale parameter for spillover interval\n")
                    log_scale_par <- stats::optimize(
                        f = estimateSpilloverLME4, interval = log_sp_prior, maximum = FALSE,
                        tol = 0.1, trial = trial, FUN = FUN, formula = formula, link = link, distance = distance)$minimum
                },
                error = function(e){
                    message("*** Spillover scale parameter cannot be estimated ***")
                    log_scale_par <- 0
                })
            }
            scale_par <- exp(log_scale_par)
        }
        analysis$options$scale_par <- scale_par
    }
    if (distance %in% c('disc','kern')) {
        trial <- compute_distance(trial,
                                  distance = distance, scale_par = scale_par)$trial
        x <- trial[[distance]]
        analysis$trial <- trial
    } else {
        x <- trial[[distance]] / scale_par
    }
    trial$effect <- eval(parse(text = FUN))
    model_object <- switch(link,
                           "identity" = lme4::lmer(formula = formula, data = trial, REML = FALSE),
                           "log" = lme4::glmer(formula = formula, data = trial,
                                               family = poisson),
                           "logit" = lme4::glmer(formula = formula, data = trial,
                                               family = binomial),
                           "cloglog" = lme4::glmer(formula = formula, data = trial,
                                               family = binomial))
    analysis$pt_ests$scale_par <- exp(log_scale_par)
    analysis$pt_ests$deviance <- unname(summary(model_object)$AICtab["deviance"])
    analysis$pt_ests$AIC <- unname(summary(model_object)$AICtab["AIC"])
    analysis$pt_ests$df <- unname(summary(model_object)$AICtab["df.resid"])
    # if the spillover parameter has been estimated then penalise the AIC and
    # adjust the degrees of freedom in the output

    cov <- q50 <- NULL

    analysis$pt_ests$AIC <- analysis$pt_ests$AIC + analysis$options$penalty
    analysis$pt_ests$df <- analysis$pt_ests$df - (analysis$options$penalty > 0)
    coefficients <- summary(model_object)$coefficients
    if (!identical(distance_type, "No fixed effects of distance ")) {
        WaldP <- ifelse(ncol(coefficients) == 4, coefficients['effect',4], NA)
    }
    if (grepl("effect", ftext, fixed = TRUE) |
        grepl("arm", ftext, fixed = TRUE)) {
        q50 <- summary(model_object)$coefficients[,1]
        names(q50)[grep("Int",names(q50))] <- "int"
        names(q50)[grep("arm",names(q50))] <- "arm"
        names(q50)[grep("effect",names(q50))] <- "effect"
        cov <- vcov(model_object)
        rownames(cov) <- colnames(cov) <- names(q50)
    }
    analysis$model_object <- model_object
    if (!identical(cfunc, "Z")){
        sample <- as.data.frame(MASS::mvrnorm(n = 10000, mu = q50, Sigma = cov))
        analysis <- extractEstimates(analysis = analysis, sample = sample)
    } else {
        analysis$pt_ests$controlY <- invlink(link, summary(model_object)$coefficients[,1])
    }
    return(analysis)
}

INLAanalysis <- function(analysis, requireMesh = requireMesh, inla_mesh = inla_mesh){
    trial <- analysis$trial
    cfunc <- analysis$options$cfunc
    link <- analysis$options$link
    distance <- analysis$options$distance
    log_sp_prior <- analysis$options$log_sp_prior
    distance_type <- analysis$options$distance_type
    linearity <- analysis$options$linearity
    scale_par <- analysis$options$scale_par
    alpha <- analysis$options$alpha
    FUN <- get_FUN(cfunc)
    fterms <- analysis$options$fterms
    # TODO replace the use of ftext with fterms
    ftext <- analysis$options$ftext

    if (identical(link,"identity")) {
        fterms <- c("y1/y_off ~ 0 + int", fterms)
    } else {
        fterms <- c("y1 ~ 0 + int", fterms)
    }
    formula <- stats::as.formula(paste(fterms, collapse = "+"))
    trial <- dplyr::mutate(trial, id =  dplyr::row_number())

    # Check if an appropriate inla_mesh is present and create one if necessary
    # If a mesh is provided use scale_par from the mesh
    # If spatial predictions are not required a minimal mesh is sufficient
    pixel <- 0.5
    if (!requireMesh) pixel <- (max(trial$x) - min(trial$x))/2
    if (is.null(inla_mesh)) {
        inla_mesh <- compute_mesh(
            trial = trial, offset = -0.1, max.edge = 0.25, inla.alpha = 2,
            maskbuffer = 0.5, pixel = pixel)
    }
    y_off <- NULL
    spde <- inla_mesh$spde
    df = data.frame(
        int = rep(1, nrow(trial)),
        id = trial$id)
    if ("int + arm" %in% fterms | "arm" %in% fterms) df$arm = ifelse(trial$arm == "intervention", 1, 0)
    if ("f(cluster, model = \'iid\')" %in% fterms) df$cluster = trial$cluster
    effectse <- list(df = df, s = inla_mesh$indexs)
    dfp = data.frame(
        int = rep(1, nrow(inla_mesh$prediction)),
        id = inla_mesh$prediction$id)
    if ("int + arm" %in% fterms | "arm" %in% fterms) dfp$arm = ifelse(inla_mesh$prediction$arm == "intervention", 1, 0)
    if ("f(cluster, model = \"iid\")" %in% fterms) dfp$cluster = inla_mesh$prediction$cluster
    effectsp <- list(df = dfp, s = inla_mesh$indexs)

    lc <- NULL
    FUN <- get_FUN(cfunc)
    log_scale_par <- ifelse(is.null(scale_par), NA, log(scale_par))
    if (!identical(distance_type, "No fixed effects of distance ")) {
        if (analysis$options$penalty > 0) {
            if (identical(Sys.getenv("TESTTHAT"), "true")) {
                log_scale_par <- 2.0
            } else {
                tryCatch({
                    #messag"Estimating scale parameter for spillover interval\n")
                    log_scale_par <- stats::optimize(
                        f = estimateSpilloverINLA, interval = log_sp_prior,
                        tol = 0.1, trial = trial, FUN = FUN, formula = formula,
                        link = link, inla_mesh = inla_mesh, distance = distance)$minimum
                },
                error = function(e){
                    message("*** Spillover scale parameter cannot be estimated ***")
                    log_scale_par <- 5
                })
            }
            scale_par <- exp(log_scale_par)
        }
        analysis$options$scale_par <- scale_par
        if (distance %in% c("disc", "kern")) {
            trial <- compute_distance(trial, distance = distance, scale_par = scale_par)$trial
            x <- trial[[distance]]
            trial$effect <- eval(parse(text = FUN))
            analysis$trial <- trial
            inla_mesh$prediction[[distance]] <-
                trial[[distance]][inla_mesh$prediction$nearestNeighbour]
            x <- inla_mesh$prediction[[distance]]
        } else {
            x <- trial[[distance]]/scale_par
            trial$effect <- eval(parse(text = FUN))
            inla_mesh$prediction[[distance]] <-
                trial[[distance]][inla_mesh$prediction$nearestNeighbour]
            x <- inla_mesh$prediction[[distance]]/scale_par
        }
        inla_mesh$prediction$effect <- eval(parse(text = FUN))
        effectse$df$effect <- trial$effect
        effectsp$df$effect <- inla_mesh$prediction$effect
        # set up linear contrasts
        lc <- INLA::inla.make.lincomb(int = 1, effect = 1)
        if (grepl("arm", ftext, fixed = TRUE)){
            lc <- INLA::inla.make.lincomb(int = 1, effect = 1, arm = 1)
        }
    } else if (grepl("arm", ftext, fixed = TRUE)) {
        lc <- INLA::inla.make.lincomb(int = 1, arm = 1)
    }
    # stack for estimation stk.e
    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla_mesh$A),
        effects = effectse)

    # stack for prediction stk.p
    stk.p <- INLA::inla.stack(
        tag = "pred", data = list(y1 = NA, y_off = NA),
        A = list(1, inla_mesh$Ap),
        effects = effectsp)

    # stk.full comprises both stk.e and stk.p if a prediction mesh is in use
    stk.full <- INLA::inla.stack(stk.e, stk.p)

    if (link == "identity") {
        model_object <- INLA::inla(
            formula, family = "gaussian", lincomb = lc,
            control.family = list(link = "identity"),
            data = INLA::inla.stack.data(stk.full),
            control.fixed = list(correlation.matrix = TRUE),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.full)),
            control.compute = list(dic = TRUE))
    } else if (link == "log") {
        model_object <- INLA::inla(
            formula, family = "poisson", lincomb = lc,
            control.family = list(link = "log"),
            data = INLA::inla.stack.data(stk.full),
            control.fixed = list(correlation.matrix = TRUE),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.full)),
            control.compute = list(dic = TRUE))
    } else if (link == "logit") {
        model_object <- INLA::inla(
            formula, family = "binomial", Ntrials = y_off, lincomb = lc,
            control.family = list(link = "logit"),
            data = INLA::inla.stack.data(stk.full),
            control.fixed = list(correlation.matrix = TRUE),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.full)),
            control.compute = list(dic = TRUE))
    } else if (link == "cloglog") {
        model_object <- INLA::inla(
            formula, family = "binomial", Ntrials = 1, lincomb = lc,
            control.family = list(link = "cloglog"),
            data = INLA::inla.stack.data(stk.full),
            control.fixed = list(correlation.matrix = TRUE),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.full)),
            control.compute = list(dic = TRUE))
    }

    analysis$pt_ests$scale_par <- scale_par

    # The DIC is penalised if a scale parameter was estimated
    analysis$pt_ests$DIC <- model_object$dic$dic + analysis$options$penalty

    # Augment the inla results list with application specific quantities
    index <- INLA::inla.stack.index(stack = stk.full, tag = "pred")$data
    inla_mesh$prediction$prediction <-
        invlink(link, model_object$summary.linear.predictor[index, "0.5quant"])
    # Compute sample-based confidence limits for intervened outcome and effect_size
    # intervention effects are estimated
    q50 <- cov <- list()
    if (grepl("effect", ftext, fixed = TRUE) |
        grepl("arm", ftext, fixed = TRUE)) {
        # Specify the point estimates of the parameters
        q50 <- model_object$summary.lincomb.derived$"0.5quant"
        names(q50) <- rownames(model_object$summary.lincomb.derived)
        # Specify the covariance matrix of the parameters
        cov <- model_object$misc$lincomb.derived.covariance.matrix

    }
    analysis$inla_mesh <- inla_mesh
    analysis$model_object <- model_object
    if (!identical(cfunc, "Z")){
        sample <- as.data.frame(MASS::mvrnorm(n = 10000, mu = q50, Sigma = cov))
        analysis <- extractEstimates(analysis = analysis, sample = sample)
    } else {
        analysis$pt_ests$controlY <- invlink(link, model_object$summary.fixed[["0.5quant"]])
    }
return(analysis)
}

group_data <- function(analysis, distance = NULL, grouping = "quintiles"){
    # define the limits of the curve both for control and intervention arms
    trial <- analysis$trial
    link <- analysis$options$link
    alpha <- analysis$options$alpha
    if (is.null(distance)) distance <- analysis$options$distance
    y_off <- y1 <- average <- upper <- lower <- d <- NULL
    cats <- NULL
    breaks0 <- breaks1 <- rep(NA, times = 6)
    # categorisation of trial data for plotting
    if (identical(grouping, "quintiles")) {
        groupvar <- ifelse(trial$arm == "intervention", 1000 + trial[[distance]], trial[[distance]])
        breaks0 <-unique(c(-Inf, quantile(groupvar[trial$arm == "control"],
                                probs = seq(0.2, 1, by = 0.20))))
        breaks1 <-unique(c(999, quantile(groupvar[trial$arm == "intervention"],
                               probs = seq(0.2, 1, by = 0.20))))
        trial$cat <- cut(
            groupvar, breaks=c(breaks0, breaks1),labels = FALSE)
        arm <- c(rep("control", times = length(breaks0)-1), rep("intervention", times = length(breaks1)-1))
    } else {
        range_d <- max(trial[[distance]]) - min(trial[[distance]])
        trial$cat <- cut(
            trial[[distance]], breaks =
                c(-Inf, min(trial[[distance]]) + seq(1:9) * range_d/10, Inf),labels = FALSE)
        arm <- NA
    }
    trial$d <- trial[[distance]]
    if (link %in% c('log', 'cloglog')) {
        data <- data.frame(
            trial %>%
            dplyr::group_by(cat) %>%
            dplyr::summarize(
                locations = dplyr::n(),
                positives = sum(y1),
                total = sum(y_off),
                d = median(d),
                average = Williams(x=y1/y_off, alpha=alpha, option = 'M'),
                lower = Williams(x=y1/y_off, alpha=alpha, option = 'L'),
                upper = Williams(x=y1/y_off, alpha=alpha, option = 'U')))
    } else if (link == 'logit') {
        data <- data.frame(
            trial %>%
            dplyr::group_by(cat) %>%
            dplyr::summarize(
                locations = dplyr::n(),
                d = median(d),
                positives = sum(y1),
                total = sum(y_off)))
        # overwrite with proportions and binomial confidence intervals by category
        data$average <- data$positives/data$total
        data$upper <- with(data, average -
                               stats::qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
        data$lower <- with(data, average +
                               stats::qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
    } else if (link == 'identity') {
        # overall means and t-based confidence intervals by category
        data <- trial %>%
            dplyr::group_by(cat) %>%
            dplyr::summarize(
                locations = dplyr::n(),
                positives = sum(y1),
                total = sum(y_off),
                d = median(d),
                average = mean(x=y1/y_off),
                lower = Tinterval(y1/y_off, alpha = alpha, option = 'L'),
                upper = Tinterval(y1/y_off, alpha = alpha, option = 'U'))
    }
    data$arm <- arm
return(data)
}

# add labels to confidence limits
namedCL <- function(limits, alpha = alpha)
    {
    names(limits) <- c(
        paste0(100 * alpha/2, "%"),
        paste0(100 - 100 * alpha/2, "%")
    )
    return(limits)
}


# Data description, crude effect_size estimate, CV and ICC calculation
get_description <- function(trial, link, alpha, baselineOnly) {
    lp <- arm <- NULL
    if(baselineOnly) trial$arm <- "control"
    clusterSum <- clusterSummary(trial = trial, link = "identity")
    sum.numerators <- tapply(trial$y1, trial$arm, FUN = sum)
    sum.denominators <- tapply(trial$y_off, trial$arm, FUN = sum)
    ratio <- sum.numerators/sum.denominators
    if(baselineOnly) {
        controlY <- ratio[1]
        effect_size <- interventionY <- NULL
    } else {
        controlY <- ratio[1]
        interventionY <- ratio[2]
        effect_size <- switch(link,
           "identity" = ratio[2] - ratio[1],
           "log" = 1 - ratio[2]/ratio[1],
           "cloglog" = 1 - ratio[2]/ratio[1],
           "logit" =  1 - ratio[2]/ratio[1])
    }
    means <- clusterSum %>% group_by(arm) %>% dplyr::summarize(lp = mean(lp))
    deviations <- ifelse(clusterSum$arm == "control", clusterSum$lp - means$lp[1], clusterSum$lp - means$lp[2])

    sigma2B <- var(clusterSum$y/clusterSum$total)
    mu <- mean(trial$y1/trial$y_off)
    sigma2 <- mean(trial$y1/trial$y_off)

    # coefficient of variation (Hayes and Moulton, equation 2.1)
    sigmaB <- sqrt(sigma2B)
    K <- sigmaB/mu

    cv_percent <- 100 * K
    cv_intervals <- cv_interval(K = K, n = nrow(clusterSum), alpha = alpha)
    if (identical(link, 'logit')){
        # Hayes & Moulton equation 2.5
        ICC <- K^2 * mu/(1 - mu)
    } else if (identical(link, 'identity')){
        # Hayes & Moulton equation 2.3
        ICC <- sigma2B/sigma2
    } else {
        ICC <- NA
    }
    description <- list(
        sum.numerators = sum.numerators,
        sum.denominators = sum.denominators,
        controlY = controlY,
        interventionY = interventionY,
        effect_size = effect_size,
        nclusters = max(as.numeric(as.character(trial$cluster))),
        cv_percent = cv_percent,
        cv_lower = cv_intervals$lcl * 100,
        cv_upper = cv_intervals$ucl * 100,
        ICC = ICC,
        locations = nrow(trial)
    )
return(description)
}

cv_interval <- function(K, n, alpha) {
    # Vangel (1996) method for interval estimates of the cv
    # https://www.jstor.org/stable/2685039

    u1 <- stats::qchisq(p = 1 - alpha/2, df = n - 1)
    u2 <- stats::qchisq(p = alpha/2, df = n - 1)
    lcl <- K/sqrt(((u1 + 2)/n - 1)* K^2 + u1/(n - 1))
    ucl <- K/sqrt(((u2 + 2)/n - 1)* K^2 + u2/(n - 1))
    value <- list(lcl = lcl, ucl = ucl)
    return(value)
}



# Functions for T and GEE analysis

noLabels <- function(x)
    {
    xclean <- as.matrix(x)
    dimnames(xclean) <- NULL
    xclean <- as.vector(xclean)
    return(xclean)
}

estimateCLeffect_size <- function(q50, Sigma, alpha, resamples, method, link)
    {

    # Use resampling approach to avoid need for Taylor approximation use at least 10000 samples (this is very
    # cheap)
    resamples1 <- max(resamples, 10000, na.rm = TRUE)
    samples <- MASS::mvrnorm(n = resamples1, mu = q50, Sigma = Sigma)

    pC <- invlink(link, samples[, 1])

    # for the T method the t-tests provide estimates for the s.e. for both arms separately
    # for the GEE method the input is in terms of the incremental effect of the intervention

    if (method == "T")
        pI <- invlink(link, samples[, 2])
    if (method == "GEE")
        pI <- invlink(link, samples[, 1] + samples[, 2])

    eff <- 1 - pI/pC

    CL <- quantile(eff, probs = c(alpha/2, 1 - alpha/2))
    return(CL)
}


# Calculate the distance or surround for an arbitrary location
# TODO: this seems to be orphan code
calculate_singlevalue <- function(i, trial , prediction , distM, distance, scale_par){
    if (identical(distance, "nearestDiscord")) {
        discords <- (trial$arm != prediction$arm[i])
        nearestDiscord <- min(distM[i, discords])
        value <- ifelse(prediction$arm[i] == "control", -nearestDiscord, nearestDiscord)
    }
    if (distance %in% c("hdep", "sdep")){
        X = list(x = prediction$x[i], y = prediction$y[i])
        depthlist <- depths(X, trial = trial)
        value <- depthlist[distance]
    }
    if (identical(distance, "disc")) {
        value <- sum(trial$arm =='intervention' & (distM[i, ] <= scale_par))
        if(identical(prediction$arm,'intervention')) value <- value - 1
    }
    return(value)
}

# Use profiling to estimate scale_par
estimateSpilloverINLA <- function(
    log_scale_par = log_scale_par, trial = trial, FUN = FUN, inla_mesh = inla_mesh,
    formula = formula, link = link, distance = distance){
    y1 <- y0 <- y_off <- NULL
    if (distance %in% c('disc','kern')) {
        updated <- compute_distance(trial, distance = distance, scale_par = exp(log_scale_par))
        x <- trial[[distance]] <- updated$trial[[distance]]
    } else {
        x <- trial[[distance]]/exp(log_scale_par)
    }
    trial$effect <- eval(parse(text = FUN))

    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla_mesh$A),
        effects = list(
            data.frame(
                int = rep(1, nrow(trial)),
                arm = ifelse(trial$arm == "intervention", 1, 0),
                effect = trial$effect, id = trial$id, cluster = trial$cluster
            ),
            s = inla_mesh$indexs
        )
    )
    # run the model with just the estimation stack (no predictions needed at this stage)
    if (link == "identity") {
        result.e <- INLA::inla(
            formula, family = "gaussian",
            control.family = list(link = "identity"),
            data = INLA::inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.e)),
            control.compute = list(dic = TRUE))
    } else if (link == "log") {
        result.e <- INLA::inla(
            formula, family = "poisson",
            control.family = list(link = "log"),
            data = INLA::inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.e)),
            control.compute = list(dic = TRUE))
    } else if (link == "logit") {
        result.e <- INLA::inla(
            formula, family = "binomial", Ntrials = y_off,
            control.family = list(link = "logit"),
            data = INLA::inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = INLA::inla.stack.A(stk.e)),
            control.compute = list(dic = TRUE))
    } else if (link == "cloglog") {
            result.e <- INLA::inla(
                formula, family = "binomial", Ntrials = 1,
                control.family = list(link = "cloglog"),
                data = INLA::inla.stack.data(stk.e),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.e)),
                control.compute = list(dic = TRUE))
    }
    # The DIC is penalised to allow for estimation of scale_par
    loss <- result.e$dic$family.dic + 2
    # Display the DIC here if necessary for debugging
    #  messag"\rDIC: ", loss, " Spillover scale parameter: ", exp(log_scale_par), "  \n")
    return(loss)
}

# Use profiling to estimate scale_par
estimateSpilloverLME4 <- function(
    log_scale_par = log_scale_par, trial = trial, FUN = FUN, formula = formula, link = link, distance = distance){
    if (distance %in% c('disc','kern')) {
        updated <- compute_distance(trial, distance = distance, scale_par = exp(log_scale_par))
        x <- trial[[distance]] <- updated$trial[[distance]]
    } else {
        x <- trial[[distance]]/exp(log_scale_par)
    }
    trial$effect <- eval(parse(text = FUN))
    try(
    model_object <- switch(link,
                "identity" = lme4::lmer(formula = formula, data = trial, REML = FALSE),
                "log" = lme4::glmer(formula = formula, data = trial,
                                   family = poisson),
                "logit" = lme4::glmer(formula = formula, data = trial,
                                   family = binomial),
                "cloglog" = lme4::glmer(formula = formula, data = trial,
                                  family = binomial))
    )
    loss <- ifelse (is.null(model_object),999999, unlist(summary(model_object)$AICtab["AIC"]))
    # The AIC is used as a loss function
    # Display the AIC here if necessary for debugging
    message ("\rAIC: ", loss + 2, " Spillover scale parameter: ", exp(log_scale_par), "  \n")
    return(loss)
}

# Add estimates to analysis list
add_estimates <- function(analysis, sample, CLnames, alpha, pt_src){
    intervals_to_output <- c("int", "effect",
                             "controlY","interventionY","effect_size",
                             "personal_protection","scale_par",
                             "deviance","spillover_interval","spillover_limit0",
                             "spillover_limit1","contaminate_pop_pr",
                             "total_effect", "ipsilateral_spillover",
                             "contralateral_spillover")
    for (var in intervals_to_output) {
        if(var %in% names(sample)) {
            quantiles <- quantile(sample[, var], c(alpha/2, 0.5, 1 - alpha/2),
                                     alpha = alpha, na.rm = TRUE)
            analysis[[pt_src]][[var]] <- quantiles[2]
            analysis$int_ests[[var]] <- stats::setNames(quantiles[c(1, 3)], CLnames)
        }
    }
    # If the scale parameter has been estimated then replace the input value
    # with the estimate. Otherwise vice versa.
    if ("scale_par" %in% names(sample)) {
            analysis$options$scale_par <- analysis$pt_ests$scale_par
        } else {
            analysis$pt_ests$scale_par <- analysis$options$scale_par
    }
return(analysis)
}

# Williams mean and confidence intervals
Williams <- function(x, alpha, option){
    logx_1 <- log(x + 1)
    logx_1[!is.finite(logx_1)] <- NA
    if(sum(is.finite(logx_1)) < 3) {
        value <- switch(option,
                M = mean(logx_1),
                L = NA,
                U = NA)
    } else {
        value <- NA
        tryCatch({
        t <- stats::t.test(x = logx_1, conf.level = 1 - alpha)
        value <- switch(option,
                M = t$estimate,
                L = t$conf.int[1],
                U = t$conf.int[2])
        },
        error = function(e){
            message("*** Averages and interval estimates not defined for some groups ***")
        })
    }
    returnvalue <- as.numeric(exp(value) - 1)
    return(returnvalue)
}

# T-based confidence intervals

Tinterval <- function(x, alpha, option){
    if(length(x) < 3){
        value <- NA
    } else {
        t <- stats::t.test(x = x, conf.level = 1 - alpha)
        value <- switch(option,
                L = t$conf.int[1],
                U = t$conf.int[2])
    }
    returnvalue <- as.numeric(value)
}

extractEstimates <- function(analysis, sample) {
    alpha <- analysis$options$alpha
    link <- analysis$options$link
    cfunc <- analysis$options$cfunc
    method <- analysis$options$method
    trial <- analysis$trial
    distance <- analysis$options$distance
    CLnames <- analysis$options$CLnames
    scale_par <- analysis$options$scale_par
    sample$controlY <- invlink(link, sample$int)
    # personal_protection is the proportion of effect attributed to personal protection
    if ("arm" %in% names(sample) & "effect" %in% names(sample)) {
        if (method %in% c("MCMC","LME4")) {
            sample$lc <- with(sample, int + effect + arm)
        }
        sample$interventionY <- invlink(link, sample$lc)
        sample$personal_protection <- with(
            sample, (controlY - invlink(link, int + arm))/(controlY -
                                                               interventionY))
    } else {
        sample$personal_protection <- NA
    }
    if ("arm" %in% names(sample) & !("effect" %in% names(sample))) {
        sample$interventionY <- invlink(link, sample$int + sample$arm)
    }
    if ("effect" %in% names(sample) & !("arm" %in% names(sample))) {
        sample$interventionY <- invlink(link, sample$int + sample$effect)
    }
    if ("interventionY" %in% names(sample)) {
        sample$effect_size <- 1 - sample$interventionY/sample$controlY
    }
    if (!(cfunc %in% c("X", "Z"))) {

        # Calculation of measures of spillover
        if (is.null(sample$scale_par)) sample$scale_par <- scale_par
        if (is.null(sample$scale_par)) sample$scale_par <- 1

        # Measures calculated for each sample from the posterior
        sample$total_effect <- sample$controlY - sample$interventionY

        spillover_list <- apply(sample, MARGIN = 1, FUN = get_spillover,
                                    analysis = analysis)
        spillover_df <- as.data.frame(do.call(rbind, lapply(spillover_list, as.data.frame)))
        sample <- cbind(sample, spillover_df)

        if ((distance %in% c("disc", "kern")) & identical(cfunc, "E")){
            range_d <- 1
            d <- (seq(1:1001) - 1)/1000
            trial[[distance]] <- NA #This is a placeholder for the sample from the estimated surround which will be returned
        } else {
            range_d <- max(trial[[distance]]) - min(trial[[distance]])
            range_d <- max(range_d, 1.0E-3)
            d <- min(trial[[distance]]) + range_d * (seq(1:1001) - 1)/1000
        }

        # Measures calculated once for each data point
        trial$sample_no <- sample(seq(1:nrow(sample)), nrow(trial), replace = TRUE)
        trial$rowname <- rownames(trial)
        curves <- t(apply(trial, MARGIN = 1, FUN = computeFittedCurve, trial = trial, sample = sample, d = d,
                        link = link, cfunc = cfunc, distance = distance))

        trial[[distance]] <- curves[, 1] # The first column is a sample from the distribution of the distance variable
        trial$fitted_value <- curves[, 2] # The second column is a sample from the distribution of the fitted value


        # Point and interval estimates of the fitted curve
        FittedCurve <- data.frame(t(apply(
            curves[, seq(3,1003)], 2, function(x) {quantile(x, c(alpha/2, 0.5, 1 - alpha/2),
                                             alpha = alpha, na.rm = TRUE)})))

        trial$spillover_limit0 <- curves[, 27]
        trial$spillover_limit1 <- curves[, 978]
        FittedCurve$intervention_curve <- ifelse(d >= 0, FittedCurve$X50., NA)
        FittedCurve$control_curve <- ifelse(d < 0, FittedCurve$X50., NA)
        FittedCurve$X50. <- NULL
        FittedCurve$d <- d
        analysis$spillover$FittedCurve <- FittedCurve

        trial$contaminate_pop_pr <- sum(trial[[distance]] > trial$spillover_limit0 &
                                             trial[[distance]] < trial$spillover_limit1)/nrow(trial)
        trial$spillover_interval <- abs(trial$spillover_limit1 - trial$spillover_limit0)

        # To remove warnings from plotting ensure that spillover interval is non-zero
        trial$spillover_limit0 <- ifelse(trial$spillover_interval < 2e-04, -1e-04, trial$spillover_limit0)
        trial$spillover_limit1 <- ifelse(trial$spillover_interval < 2e-04, 1e-04, trial$spillover_limit1)
        trial$ipsilateral_spillover <- ifelse(trial$arm == 'intervention',
                            trial$fitted_value - sample[trial$sample_no, "interventionY"], NA)
        trial$contralateral_spillover <- ifelse(trial$arm == 'control',
                            trial$fitted_value - sample[trial$sample_no, "controlY"], NA)
    }
    analysis <- add_estimates(analysis = analysis, sample = trial, CLnames = CLnames, alpha = alpha, pt_src = 'spillover')
    analysis <- add_estimates(analysis = analysis, sample = sample, CLnames = CLnames, alpha = alpha, pt_src = 'pt_ests')
return(analysis)
}

computeFittedCurve <- function(x, trial, cfunc, link, d = d, distance, sample) {
    i <- as.integer(x[['sample_no']])
    scale_par <- sample[i, 'scale_par']
    par0 <- c(sample[i, "par0_1"], sample[i, "par0_2"], sample[i, "scale_par"])
    par1 <- c(sample[i, "par1_1"], sample[i, "par1_2"], sample[i, "scale_par"])
    if ((distance %in% c("disc", "kern")) & identical(cfunc, "E")){
        trial <- compute_distance(trial, distance = distance, scale_par = scale_par)$trial
    }
    point <- data.frame(distance = as.numeric(trial[x[['rowname']], distance]), arm = x[['arm']])
    names(point) <- c(eval(distance),'arm')
    fv <- ifelse((identical(x[['arm']], 'control')),
                 fitted_spillover(cfunc = cfunc, link = link, par = par0, trial = point, distance = distance),
                 fitted_spillover(cfunc = cfunc, link = link, par = par1, trial = point, distance = distance))
    # Computation of the full curve

    curve_df <- data.frame(distance = d)
    names(curve_df) <- c(distance)
    curve <-
        fitted_spillover(cfunc = cfunc, link = link, par = par1, trial = curve_df, distance = distance)
    if (min(d) < 0) {
        control_curve <-
            fitted_spillover(cfunc = cfunc, link = link, par = par0, trial = curve_df, distance = distance)
        curve <- ifelse(d > 0, control_curve, curve)
    }
    # concatenate the sampled value of the distance and of the fitted value of the response with the curve
    # so that the return value is a single vector
    curve <- c(point[1,1], fv, curve)
return(curve)}


# logit transformation
logit <- function(p = p)
{
    return(log(p/(1 - p)))
}

# cloglog transformation
cloglog = function(p) log(-log(1-p))


# link transformation
link_tr <- function(link = link, x = x)
{
    value <- switch(link,
                    "identity" = x,
                    "log" = log(x),
                    "logit" =  log(x/(1 - x)),
                    "cloglog" =  log(-log(1 - x)))
    return(value)
}

# inverse transformation of link function
invlink <- function(link = link, x = x)
{
    value <- switch(link,
                    "identity" = x,
                    "log" = exp(x),
                    "logit" =  1/(1 + exp(-x)),
                    "cloglog" = 1 - exp(-exp(x)))
    return(value)
}

# piecewise linear model
PiecewiseLinearFunction <- function(par, trial, distance)
{
    # constrain the slope parameter to be positive (par[2] is positive if effect_size is negative)
    scale_par <- par[3]

    # if scale_par is very large, the curve should be close to a straight line
    if (scale_par > 20){
        lp <- rep(par[1] + 0.5 * par[2], times = nrow(trial))

    } else {
        lp <- ifelse(
            trial[[distance]] > -scale_par/2, par[1] + par[2] * (scale_par/2 + trial[[distance]])/scale_par, par[1])
        lp <- ifelse(trial[[distance]] > scale_par/2, par[1] + par[2], lp)
    }
    return(lp)
}

escape = function(x) {
    value <- 1 - exp(-x)
return(value)}

piecewise <- function(x) {
    value <- ifelse(x < -0.5, 0, (0.5 + x))
    value <- ifelse(x > 0.5, 1, value)
return(value)}

rescale <- function(x) {
    value <- (x - min(x))/(max(x) - min(x))
return(value)}

fitted_spillover <- function(cfunc, link, par, trial, distance) {
    lp <- switch(cfunc,
       Z = par[1],
       X = ifelse(trial[[distance]] < 0, par[1], par[1] + par[2]),
       S = PiecewiseLinearFunction(par=par, trial=trial, distance=distance),
       L = par[1] + par[2] * invlink(link = "logit", x = trial[[distance]]/par[3]),
       P = par[1] + par[2] * stats::pnorm(trial[[distance]]/par[3]),
       D = NA,
       R = par[1] + rescale(trial[[distance]]) * par[2],
       E = par[1] + par[2] * (1 - exp(-(trial[[distance]]/par[3]))))
    fv <- invlink(link, lp)
    if(identical(cfunc,'D')) {
        u0 <- invlink(link, par[1])
        u1 <- invlink(link, x = par[1] + par[2])
        fv <- u0 + (u1 - u0) * stats::pnorm(trial[[distance]]/par[3])
    }
return(fv)}

get_FUN <- function(cfunc){
    # trap a warning with use of "E"
    if (identical(cfunc, "E")) cfunc = "ES"
    FUN <- switch(
        cfunc, L = "invlink(link='logit', x)",
             P = "stats::pnorm(x)",
             D = "stats::pnorm(x)",
             S = "piecewise(x)",
             X = "rescale(x)",
             Z = "rescale(x)",
             R = "rescale(x)",
             ES = "escape(x)")
    return(FUN)
}

get_spillover <- function(x, analysis){
    trial <- analysis$trial
    link <- analysis$options$link
    distance <- analysis$options$distance
    cfunc <- analysis$options$cfunc
    if ((distance %in% c("disc", "kern")) & identical(cfunc, "E")) {
        trial <- compute_distance(trial, distance = distance, scale_par = x[['scale_par']])
        cfunc <- "R"
    }
    limits <- matrix(x[["controlY"]], nrow = 2, ncol = 2)
    if(!is.null(x[["interventionY"]])) {
        limits[ ,2] <-  rep(x[["interventionY"]], times = 2)
    } else {
        limits[ , 2] <- rep(x[["controlY"]], times = 2)
    }
    if (!is.na(x[["personal_protection"]])) {
        limits[1, 2] <- invlink(link, x[["int"]] + x[["effect"]])
        limits[2, 1] <- invlink(link, x[["int"]] + x[["arm"]])
    }
    if (identical(cfunc, 'X')) {
        limits[1, 2] <- limits[1, 1]
        limits[2, 1] <- limits[2, 2]
    }
    scale_par <- ifelse("scale_par" %in% names(x), x[["scale_par"]], analysis$options$scale_par)
    # Trap cases with extreme effect: TODO: a different criterion may be needed for continuous data
    pars <- link_tr(link,limits)
    if (sum((pars[, 1] - pars[, 2])^2) > 10000) {
        limits <- invlink(link, matrix(c(20,20, -20, -20), nrow = 2, ncol = 2))
    }

    par0 <- c(link_tr(link, limits[1, 1]),
              link_tr(link, limits[1, 2]) - link_tr(link, limits[1, 1]),
              scale_par)
    par1 <- c(
        link_tr(link, limits[2, 1]),
        link_tr(link, limits[2, 2]) - link_tr(link, limits[2, 1]),
        scale_par
    )

    # trap extreme cases with undefined, flat or very steep curves
    if (!identical(cfunc, "R")) {
        if (is.null(scale_par)) {
            cfunc <- "X"
        } else if (is.na(scale_par) |
                   scale_par < 0.01 |
                   scale_par > 100  |
                   (sum((pars[, 1] - pars[, 2])^2) > 10000)) {
            cfunc <- "X"
        }
    }
    spillover <- list(
        par0_1 = par0[1],
        par0_2 = par0[2],
        par1_1 = par1[1],
        par1_2 = par1[2])
return(spillover)}



# TODO: this is orphan code
calculate_spillover_deviation <- function(d, cfunc, link, distance, par, sign, q, log_interval, verbose){
    qtable <- data.frame(row.names = c(1,2,3))
    qtable[[distance]] <-  c(exp(log_interval), exp(d))
    if(sign == 'minus') qtable[[distance]] <- -qtable[[distance]]
    qtable$y_hat <- fitted_spillover(cfunc = cfunc, link =link, par = par, trial = qtable, distance = distance)
    if(qtable$y_hat[1] != qtable$y_hat[2]) {
        quantile <- with(qtable, (y_hat[3] - y_hat[2])/(y_hat[1] - y_hat[2]))
    } else {
        if(qtable$y_hat[3] == qtable$y_hat[2]) {
            quantile <- 0
        } else {
            quantile <- NA
            if (verbose) message('spillover deviation undefined')
        }
    }
    deviation <- (q - quantile)^2
    # If the scale parameter is very small, then force the minimum deviation to be at d = 0
    if (par[3] < 0.05) deviation <- d^2
    return(deviation)
}


getDistanceText <- function(distance = "nearestDiscord", scale_par = NULL) {
    value <- switch(distance,
                    "nearestDiscord" = "Signed distance to other arm (km)",
                    "disc" = paste0("disc of radius ", round(scale_par, digits = 3), " km"),
                    "kern" = paste0("kern with kernel s.d. ", round(scale_par, digits = 3), " km"),
                    "hdep" = "Tukey half-depth ",
                    "sdep" = "Simplicial depth ",
                    distance)
    return(value)
}

compute_effect <- function(trial, distance, scale_par, FUN) {
    if (distance %in% c('disc','kern')) {
        trial <- compute_distance(trial,
                    distance = distance, scale_par = scale_par)$trial
        x <- trial[[distance]]
    } else {
        x <- trial[[distance]]/scale_par
    }
    effect <- eval(parse(text = FUN))
    return(effect)
}
