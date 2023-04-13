#' Analysis of cluster randomized trial with contamination
#'
#' \code{CRTanalysis} carries out a statistical analysis of a cluster randomized trial (CRT).
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}) and outcome data (see details).
#' @param method statistical method with options:
#'  \tabular{ll}{
#' \code{"EMP"} \tab simple averages of the data   \cr
#' \code{"T"}   \tab comparison of cluster means by t-test \cr
#' \code{"GEE"} \tab Generalised Estimating Equations \cr
#' \code{"LME4"} \tab Generalized Linear Mixed-Effects Models \cr
#' \code{"INLA"}\tab Integrated Nested Laplace Approximation (INLA) \cr
#' \code{"MCMC"}\tab Markov chain Monte Carlo using \code{"JAGS"} \cr
#' }
#' @param cfunc transformation defining the contamination function. \cr
#' options are:
#' \tabular{llll}{
#' \code{"Z"} \tab\tab arm effects not considered\tab reference model\cr
#' \code{"X"} \tab\tab contamination not modelled\tab the only valid value of \code{cfunc} for methods \code{"EMP"}, \code{"T"} and \code{"GEE"}\cr
#' \code{"L"} \tab\tab inverse logistic (sigmoid)\tab the default for \code{"INLA"} and \code{"MCMC"} methods\cr
#' \code{"P"} \tab\tab inverse probit (error function)\tab available with \code{"INLA"} and \code{"MCMC"} methods\cr
#' \code{"S"} \tab\tab piecewise linear\tab only available with the \code{"MCMC"} method\cr
#' }
#' @param link link function. options are:
#'  \tabular{ll}{
#' \code{"logit"}\tab (the default). \code{numerator} has a binomial distribution with denominator \code{denominator}.\cr
#' \code{"log"}  \tab \code{numerator} is Poisson distributed with an offset of log(\code{denominator}).
#' With the \code{"INLA"} and \code{"MCMC"} methods 'iid' random effects are used to model extra-Poisson variation.\cr
#' \code{"identity"}\tab The outcome is \code{numerator/denominator} with a normally distributed error function.\cr
#' }
#' @param numerator string: name of numerator variable for outcome
#' @param denominator string: name of denominator variable for outcome data (if present)
#' @param excludeBuffer logical: indicator of whether any buffer zone (records with \code{buffer=TRUE}) should be excluded from analysis
#' @param alpha numeric: confidence level for confidence intervals and credible intervals
#' @param baselineOnly logical: indicator of whether required analysis is of effect size or of baseline only
#' @param baselineNumerator string: name of numerator variable for baseline data (if present)
#' @param baselineDenominator string: name of denominator variable for baseline data (if present)
#' @param personalProtection logical: indicator of whether the model includes local effects with no contamination
#' @param clusterEffects logical: indicator of whether the model includes cluster random effects
#' @param spatialEffects logical: indicator of whether the model includes spatial random effects
#' @param resamples integer: number of samples for sample-based intervals
#' @param requireMesh logical: indicator of whether spatial predictions are required
#' @param inla_mesh string: name of pre-existing INLA input object created by \code{new_mesh()}
#' @return list of class \code{CRTanalysis} containing the following results of the analysis:
#' \itemize{
#' \item \code{description} : description of the dataset
#' \item \code{method} : statistical method
#' \item \code{pt_ests} : point estimates
#' \item \code{int_ests} : interval estimates
#' \item \code{model_object} : object returned by the fitting routine
#' \item \code{contamination} : function values and statistics describing the estimated contamination
#' }
#' @details \code{CRTanalysis} is a wrapper for the statistical analysis packages:
#' [geepack](https://www.jstatsoft.org/article/view/v015i02),
#' [INLA](https://www.r-inla.org/),
#' [jagsUI](https://cran.r-project.org/web/packages/jagsUI/index.html),
#' and the [t.test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test)
#' function of package \code{stats}.
#' The wrapper does not provide an interface to the full functionality of these packages.
#' It is specific for typical analyses of cluster randomized trials with geographical clustering.\cr \cr
#'
#' The key results of the analyses can be extracted using a \code{summary()} of the output list.
#' The \code{model_object} in the output list is the usual output from the statistical analysis routine,
#' and can be also be inspected with \code{summary()}, or analysed using \code{stats::fitted()}
#' for purposes of evaluation of model fit etc..
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @export
#' @examples
#' {example <- readdata('exampleCRT.txt')
#' # Analysis of test dataset by t-test
#' exampleT <- CRTanalysis(example, method = "T")
#' summary(exampleT)
#' # Standard GEE analysis of test dataset ignoring contamination
#' exampleGEE <- CRTanalysis(example, method = "GEE")
#' summary(exampleGEE)
#' # LME4 analysis with error function contamination function
#' exampleLME4 <- CRTanalysis(example, method = "LME4", cfunc = "P")
#' summary(exampleLME4)
#' }
CRTanalysis <- function(
    trial, method = "GEE", cfunc = "L", link = "logit", numerator = "num",
    denominator = "denom", excludeBuffer = FALSE, alpha = 0.05,
    baselineOnly = FALSE, baselineNumerator = "base_num", baselineDenominator = "base_denom",
    personalProtection = FALSE, clusterEffects = TRUE, spatialEffects = FALSE,
    resamples = 10000, requireMesh = FALSE, inla_mesh = NULL)
    {
    # Test of validity of inputs
    if (!method %in% c("EMP", "T", "MCMC", "GEE", "INLA", "LME4"))
        {
        stop("*** Invalid value for statistical method ***")
        return(NULL)
    }

    if (!cfunc %in% c("S", "L", "P", "X", "Z"))
        {
        stop("*** Invalid contamination function ***")
        return(NULL)
    }

    if (!link %in% c("logit", "log", "identity"))
    {
        stop("*** Invalid link ***")
        return(NULL)
    }

    # MAIN FUNCTION CODE STARTS HERE

    cluster <- NULL

    if (identical(class(trial),"CRTsp")) trial <- trial$trial

    if ("buffer" %in% colnames(trial) & excludeBuffer)
        {
        trial <- trial[!trial$buffer, ]
    }

    # trial needs to be ordered for some analyses
    if(!is.null(trial$cluster)) trial <- trial[order(trial$cluster), ]

    # Some statistical methods do not allow for contamination
    if (method %in% c("EMP", "T", "GEE")) cfunc <- "X"

    # Some statistical methods only run if there are cluster effects
    if (method %in% c("LME4", "MCMC")) clusterEffects <- TRUE

    if (baselineOnly){
        # Baseline analyses are available only for GEE and INLA
        if (method %in% c("EMP", "T", "GEE", "MCMC", "LME4"))
            {
            method <- "GEE"
            cat("Analysis of baseline only, using GEE\n")
        } else if (method == "INLA")
        {
            cat("Analysis of baseline only, using INLA\n")
        }
        # cfunc='Z' is used to remove the estimation of effect size from the model
        cfunc <- "Z"
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
        trial$y0 <- trial[[denominator]] - trial[[numerator]]
        trial$y_off <- trial[[denominator]]

        # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant
        # coordinate
        if (is.null(trial$nearestDiscord)) trial$nearestDiscord <- get_nearestDiscord(trial)
    }

    # create model formula for display even though this is not used for MCMC models
    fterms <- switch(cfunc,
        Z = NULL,
        X = "arm",
        L = "pvar",
        P = "pvar",
        S = "pvar"
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
        if (identical(link,"log")) fterms <- c(fterms, "f(id, model = \'iid\')")
    }
    ftext <- paste(fterms, collapse = " + ")

    # create names for confidence limits for use throughout
    CLnames <- c(
        paste0(alpha/0.02, "%"),
        paste0(100 - alpha/0.02, "%")
    )

    # store options here- noting that the model formula depends on allowable values of other options
    options <- list(method = method, link = link, cfunc = cfunc,
                    alpha = alpha, baselineOnly = baselineOnly,
                    fterms = fterms, ftext = ftext,
                    CLnames = CLnames, clusterEffects = clusterEffects)
    # create scaffolds for lists
    pt_ests <- list(contamination_par = NA, personal_protection = NA, contamination_interval = NA)
    int_ests <- list(controlY = NA, interventionY = NA, effect_size = NA)
    model_object <- list()
    description <- get_description(trial=trial, link=link, baselineOnly)
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
           "MCMC" = MCMCanalysis(analysis)
    )
    if (!baselineOnly){
        fittedCurve <- get_curve(x = analysis$pt_ests, analysis = analysis)
        contamination <- get_contaminationStats(fittedCurve=fittedCurve, trial=analysis$trial)
        analysis <- tidyContamination(contamination, analysis, fittedCurve)
    }
    class(analysis) <- "CRTanalysis"
    return(analysis)
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
    # GEE analysis of cluster effects
    fterms <- c(switch(link,
                 "identity" = "y1/y_off ~ 1",
                 "log" = "y1 ~ 1 + offset(log(y_off))",
                 "logit" = "cbind(y1,y0) ~ 1"),
                fterms)
    formula <- stats::as.formula(paste(fterms, collapse = "+"))
    if (link == "log") {
        model_object <- geepack::geeglm(
            formula = formula, id = cluster, data = trial, family = poisson(link = "log"),
            corstr = "exchangeable", scale.fix = FALSE)
    } else if (link == "logit") {
        model_object <- geepack::geeglm(
            formula = formula, id = cluster, corstr = "exchangeable",
            data = trial, family = binomial(link = "logit"))
    } else if (link == "identity") {
        model_object <- geepack::geeglm(
            formula = formula, id = cluster, corstr = "exchangeable",
            data = trial, family = gaussian)
    }

    summary.fit <- summary(model_object)

    z <- -qnorm(alpha/2)  #standard deviation score for calculating confidence intervals
    lp_yC <- summary.fit$coefficients[1, 1]
    se_lp_yC <- summary.fit$coefficients[1, 2]

    clusterSize <- nrow(trial)/nlevels(as.factor(trial$cluster))


    # remove the temporary objects from the dataframe
    model_object$data$y1 <- model_object$data$y0 <- model_object$data$y_off <- NULL
    pt_ests$controlY <- invlink(link, lp_yC)
    int_ests$controlY <- namedCL(
        invlink(link, c(lp_yC - z * se_lp_yC, lp_yC + z * se_lp_yC)),
        alpha = alpha
    )

    # Intracluster correlation
    pt_ests$ICC <- noLabels(summary.fit$corr[1])  #with corstr = 'exchangeable', alpha is the ICC
    se_ICC <- noLabels(summary.fit$corr[2])
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
            model_object$geese$vbeta[1, 1] +
            model_object$geese$vbeta[2, 2] +
            2 * model_object$geese$vbeta[1,2]
        )

        int_ests$interventionY <- namedCL(
            invlink(link, c(lp_yI - z * se_lp_yI, lp_yI + z * se_lp_yI)),
            alpha = alpha)

        int_ests$effect_size <- estimateCLeffect_size(
            q50 = summary.fit$coefficients[, 1], Sigma = model_object$geese$vbeta,
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

INLAanalysis <- function(analysis, requireMesh = requireMesh, inla_mesh = inla_mesh){
    trial <- analysis$trial
    cfunc <- analysis$options$cfunc
    link <- analysis$options$link
    alpha <- analysis$options$alpha
    FUN <- get_FUN(cfunc, variant = 0)
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
    # If spatial predictions are not required a minimal mesh is generated
    pixel <- 0.5
    if (!requireMesh) pixel <- 2
    if (is.null(inla_mesh)) {
            inla_mesh <- new_mesh(
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
    log_beta <- NA
    FUN <- get_FUN(cfunc=cfunc, variant = 0)
    if (cfunc %in% c("L", "P", "S")) {
        if (identical(Sys.getenv("TESTTHAT"), "true")) {
            log_beta <- 2.0
        } else {
#               cat("Estimating scale parameter for contamination range", "\n")
            log_beta <- stats::optimize(
                f = estimateContaminationINLA, interval = c(-10, 10),
                trial = trial, FUN = FUN, inla_mesh = inla_mesh, formula = formula,
                tol = 0.1, link = link)$minimum
        }
        beta <- exp(log_beta)
        x <- trial$nearestDiscord * beta
        trial$pvar <- eval(parse(text = FUN))

        effectse$df$pvar <- trial$pvar

        x <- inla_mesh$prediction$nearestDiscord * beta
        inla_mesh$prediction$pvar <- eval(parse(text = FUN))
        effectsp$df$pvar <- inla_mesh$prediction$pvar

        # set up linear contrasts (not required for cfunc='X' or 'Z')
        if (grepl("pvar", ftext, fixed = TRUE))
            {
            lc <- INLA::inla.make.lincomb(int = 1, pvar = 1)
            if (grepl("arm", ftext, fixed = TRUE))
              {
              lc <- INLA::inla.make.lincomb(int = 1, pvar = 1, arm = 1)
            }
        } else if (grepl("arm", ftext, fixed = TRUE))
            {
            lc <- INLA::inla.make.lincomb(int = 1, arm = 1)
        }
    }
    # stack for estimation stk.e
    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla_mesh$A),
        effects = effectse
    )

    # stack for prediction stk.p
    stk.p <- INLA::inla.stack(
        tag = "pred", data = list(y1 = NA, y_off = NA),
        A = list(1, inla_mesh$Ap),
        effects = effectsp
    )

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
    }

    analysis$pt_ests$contamination_par <- exp(log_beta)
    # The DIC is penalised if a contamination parameter was estimated
    analysis$pt_ests$DIC <- ifelse(cfunc %in% c("L", "P", "S"),
                   model_object$dic$dic + 2, model_object$dic$dic)
    # Augment the inla results list with application specific quantities
    index <- INLA::inla.stack.index(stack = stk.full, tag = "pred")$data
    inla_mesh$prediction$prediction <-
            invlink(link, model_object$summary.linear.predictor[index, "0.5quant"])
    # Compute sample-based confidence limits for intervened outcome and effect_size
    # intervention effects are estimated
    q50 <- cov <- list()
    if (grepl("pvar", ftext, fixed = TRUE) |
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
return(analysis)}


group_data <- function(analysis){
    # define the limits of the curve both for control and intervention arms
        trial <- analysis$trial
        link <- analysis$options$link
        alpha <- analysis$options$alpha
        y_off <- y1 <- average <- upper <- lower <- NULL
        cats <- nearestDiscord <- NULL

        # categorisation of trial data for plotting
        range_d <- max(trial$nearestDiscord) - min(trial$nearestDiscord)
        trial$cats <- cut(
            trial$nearestDiscord, breaks =
                c(-Inf, min(trial$nearestDiscord) + seq(1:9) * range_d/10, Inf),labels = FALSE)

        # calculate for log link
        data <- data.frame(
            trial %>%
                group_by(cats) %>%
                dplyr::summarize(
                    positives = sum(y1),
                    total = sum(y_off),
                    d = median(nearestDiscord),
                    average = Williams(x=y1/y_off, alpha=alpha, option = 'M'),
                    lower = Williams(x=y1/y_off, alpha=alpha, option = 'L'),
                    upper = Williams(x=y1/y_off, alpha=alpha, option = 'U')))
        if (link == 'logit') {
            # overwrite with proportions and binomial confidence intervals by category
            data$average <- data$positives/data$total
            data$upper <- with(data, average -
                                   qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
            data$lower <- with(data, average +
                                   qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
        }
        if (link == 'identity') {
            # overall means and t-based confidence intervals by category
            data <- trial %>%
                group_by(cats) %>%
                dplyr::summarize(
                    positives = sum(y1),
                    total = sum(y_off),
                    d = median(nearestDiscord),
                    average = mean(x=y1/y_off),
                    lower = Tinterval(y1/y_off, alpha = alpha, option = 'L'),
                    upper = Tinterval(y1/y_off, alpha = alpha, option = 'U')
                )
        }
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

# logit transformation
logit <- function(p = p)
    {
    return(log(p/(1 - p)))
}

# link transformation
link_tr <- function(link = link, x = x)
{
    value <- switch(link,
                    "identity" = x,
                    "log" = log(x),
                    "logit" =  log(x/(1 - x)))
    return(value)
}

# inverse transformation of link function
invlink <- function(link = link, x = x)
    {
    value <- switch(link,
        "identity" = x,
        "log" = exp(x),
        "logit" =  1/(1 + exp(-x)))
    return(value)
}

# Minimal data description and crude effect_size estimate
get_description <- function(trial, link, baselineOnly)
    {
    if(baselineOnly){
        sum.numerators = sum(trial$y1)
        sum.denominators = sum(trial$y_off)
        description <- list(
            sum.numerators = sum.numerators,
            sum.denominators = sum.denominators,
            controlY = sum.numerators/sum.denominators,
            interventionY = NULL,
            effect_size = NULL,
            nclusters = max(as.numeric(as.character(trial$cluster))),
            locations = nrow(trial)
        )
    } else {
        sum.numerators <- tapply(trial$y1, trial$arm, FUN = sum)
        sum.denominators <- tapply(trial$y_off, trial$arm, FUN = sum)
        ratio <- sum.numerators/sum.denominators
        effect_size <- switch(link,
               "identity" = ratio[2] - ratio[1],
               "log" = 1 - ratio[2]/ratio[1],
               "logit" =  1 - ratio[2]/ratio[1])
        description <- list(
            sum.numerators = sum.numerators,
            sum.denominators = sum.denominators,
            controlY = ratio[1],
            interventionY = ratio[2],
            effect_size = effect_size,
            nclusters = max(as.numeric(as.character(trial$cluster))),
            locations = nrow(trial)
        )
    }
    return(description)
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

# functions for INLA analysis

#' \code{new_mesh} create objects required for INLA analysis of an object of class \code{"CRTsp"}.
#' @param trial an object of class \code{"CRTsp"} or a data frame containing locations in (x,y) coordinates, cluster
#'   assignments (factor \code{cluster}), and arm assignments (factor \code{arm}) and outcome.
#' @param offset (see \code{inla.mesh.2d} documentation)
#' @param max.edge  (see \code{inla.mesh.2d} documentation)
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
#' @details \code{new_mesh} carries out the computationally intensive steps required for setting-up an
#' INLA analysis of an object of class \code{"CRTsp"}, creating the preface_validitytion mesh and the projection matrices.
#' The mesh can be reused for different models fitted to the same
#' geography. The computational resources required depend largely on the resolution of the prediction mesh.
#' The prediction mesh is thinned to include only pixels centred at a distance less than
#' \code{maskbuffer} from the nearest point.\cr
#' A warning may be generated unless the \code{Matrix} library is loaded.
#' @export
#' @examples
#' {
#' # low resolution mesh for test dataset
#' library(Matrix)
#' example <- readdata('exampleCRT.txt')
#' exampleMesh=new_mesh(example, pixel = 0.5)
#' }
#' \dontrun{
#' # 50m mesh for test dataset
#' library(Matrix)
#' example <- readdata('exampleCRT.txt')
#' exampleMesh=new_mesh(example, pixel = 0.05)
#' }
new_mesh <- function(trial = trial, offset = -0.1, max.edge = 0.25,
                     inla.alpha = 2, maskbuffer = 0.5, pixel = 0.5)
    {
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
        x = pred.coords[, 1], y = pred.coords[, 2], nearestNeighbour = nearestNeighbour)
    prediction$id <- trial$id[nearestNeighbour]
    if (!is.null(trial$arm)) prediction$arm <- trial$arm[nearestNeighbour]
    if (!is.null(trial$cluster)) prediction$cluster <- trial$cluster[nearestNeighbour]
    prediction <- with(
        prediction, prediction[order(y, x), ]
    )
    prediction$shortestDistance <- apply(distM, 1, min)
    rows <- seq(1:nrow(prediction))
    if (!is.null(trial$arm)) {
        prediction$nearestDiscord <- sapply(rows,
                                            FUN = calcNearestDiscord,
                                            trial = trial,
                                            prediction = prediction,
                                            distM = distM)
    }
    inla_mesh <- list(
        prediction = prediction, A = A, Ap = Ap, indexs = indexs, spde = spde, pixel = pixel
    )
    cat(
        "Mesh created of ", nrow(prediction), " pixels of size ", pixel," km \n"
    )

    return(inla_mesh)
}

# Calculate the distance to the nearest discordant location
calcNearestDiscord <- function(x, trial , prediction , distM)
{
    discords <- (trial$arm != prediction$arm[x])
    nearestDiscord <- min(distM[x, discords])
    nearestDiscord <- ifelse(prediction$arm[x] == "control", -nearestDiscord, nearestDiscord)
    return(nearestDiscord)
}

# Use profiling to estimate beta
estimateContaminationINLA <- function(
    log_beta = log_beta, trial = trial, FUN = FUN, inla_mesh = inla_mesh, formula = formula, link = link){
    y1 <- y0 <- y_off <- NULL
    x <- trial$nearestDiscord * exp(log_beta)
    trial$pvar <- eval(parse(text = FUN))

    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla_mesh$A),
        effects = list(
            data.frame(
                int = rep(1, nrow(trial)),
                arm = ifelse(trial$arm == "intervention", 1, 0),
                pvar = trial$pvar, id = trial$id, cluster = trial$cluster
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
            control.compute = list(dic = TRUE)
        )
    }
    # The DIC is penalised to allow for estimation of beta
    loss <- result.e$dic$family.dic
    # Display the DIC here if necessary for debugging
    # cat("\rDIC: ", loss, " Contamination scale parameter: ", exp(log_beta), "  \r")
    return(loss)
}

# Use profiling to estimate beta
estimateContaminationLME4 <- function(
    log_beta = log_beta, trial = trial, FUN = FUN, formula = formula, link = link){
    x <- trial$nearestDiscord * exp(log_beta)
    trial$pvar <- eval(parse(text = FUN))
    try(
    model_object <- switch(link,
                "identity" = lme4::lmer(formula = formula, data = trial, REML = FALSE),
                "log" = lme4::glmer(formula = formula, data = trial,
                                   family = poisson),
                "logit" = lme4::glmer(formula = formula, data = trial,
                                   family = binomial))
    )
    loss <- ifelse (is.null(model_object),999999, unlist(summary(model_object)$AICtab["AIC"]))
    # The AIC is used as a loss function
    # Display the AIC here if necessary for debugging
    # cat("\rAIC: ", loss, " Contamination scale parameter: ", exp(log_beta), "  \n")
    return(loss)
}


# Add estimates to analysis list
add_estimates <- function(analysis, bounds, CLnames){
    bounds <- data.frame(bounds)
    for (variable in c("controlY","interventionY","effect_size",
                      "personal_protection","contamination_par",
                      "deviance","contamination_interval","contamination_limit0",
                      "contamination_limit1","contaminate_pop_pr")) {
        if (variable %in% colnames(bounds)) {
            analysis$pt_ests[[variable]] <- bounds[2, variable]
            analysis$int_ests[[variable]] <- stats::setNames(
                bounds[c(1, 3), variable], CLnames)
        }
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


#' Summary of the results of a statistical analysis of a CRT
#'
#' \code{summary.CRTanalysis} generates a summary of a \code{CRTanalysis} including the main results
#' @param object an object of class \code{"CRTanalysis"}
#' @param ... other arguments used by summary
#' @method summary CRTanalysis
#' @export
#'
summary.CRTanalysis <- function(object, ...) {
    defaultdigits <- getOption("digits")
    options(digits = 3)
    cat("=====================CLUSTER RANDOMISED TRIAL ANALYSIS =================\n")
    cat(
        "Analysis method: ", object$options$method, "\nLink function: ", object$options$link, "\n"
    )
    if (!is.null(object$options$ftext))
        cat("Model formula: ", object$options$ftext, "\n")
    cat(switch(object$options$cfunc,
               Z = "No comparison of arms \n",
               X = "No modelling of contamination \n",
               S = "Piecewise linear function for contamination\n",
               P = "Error function model for contamination\n",
               L = "Sigmoid (logistic) function for contamination\n"))
    CLtext <- paste0(" (", 100 * (1 - object$options$alpha), "% CL: ")
    cat(
        "Estimates:      Control: ", object$pt_ests$controlY,
        CLtext, unlist(object$int_ests$controlY),
        ")\n"
    )
    if (!is.null(object$pt_ests$effect_size))
    {
        if (!is.null(object$pt_ests$interventionY))
        {
            cat(
                "           Intervention: ", object$pt_ests$interventionY,
                CLtext, unlist(object$int_ests$interventionY),
                ")\n"
            )
            effect.measure <- ifelse(object$options$link == 'identity', "Effect size: ","   Efficacy: ")
            cat("           ",
                effect.measure, object$pt_ests$effect_size, CLtext, unlist(object$int_ests$effect_size),
                ")\n"
            )
        }
        if (!is.na(object$pt_ests$personal_protection))
        {
            cat(
                "Personal protection %  : ",
                object$pt_ests$personal_protection*100,
                CLtext, unlist(object$int_ests$personal_protection*100),
                    ")\n"
            )
            if (object$pt_ests$personal_protection < 0 | object$pt_ests$personal_protection >  1){
                cat(
                    "** Warning: different signs for main effect and personal protection effect:
            face validity check fails **\n")
            }
        }
        if (object$options$cfunc %in% c("L", "P", "S")) {
            if (!is.null(object$pt_ests$contamination_interval)){
                cat(
                    "Contamination range(km): ", object$pt_ests$contamination_interval,
                    CLtext, unlist(object$int_ests$contamination_interval),
                    ")\n"
                )
            }
            if (!is.null(object$contamination$contaminate_pop_pr)){
                cat(
                    "% locations contaminated:",
                    object$contamination$contaminate_pop_pr*100,
                    CLtext, unlist(object$int_ests$contaminate_pop_pr)*100,
                    " %)\n")
            }
        }
    }
    if (!is.null(object$pt_ests$ICC))
    {
        cat(
            "Intracluster correlation (ICC) : ", object$pt_ests$ICC,
            CLtext, unlist(object$int_ests$ICC),")\n"
        )
    }
    options(digits = defaultdigits)
    # goodness of fit
    if (!is.null(object$pt_ests$deviance)) cat("deviance: ", object$pt_ests$deviance, "\n")
    if (!is.null(object$pt_ests$DIC)) cat("DIC     : ", object$pt_ests$DIC)
    if (!is.null(object$pt_ests$AIC)) cat("AIC     : ", object$pt_ests$AIC)
    if (object$options$method %in% c("LME4","INLA") &
            object$options$cfunc %in% c("L", "P", "S")) {
        cat(" including penalty for the contamination scale parameter\n")
    } else {
        cat(" \n")
    }
# TODO: add the degrees of freedom to the output
    if (!is.null(object$pt_ests$p.value)){
        cat("P-value (2-sided): ", object$pt_ests$p.value, "\n")
    }

}

EMPanalysis <- function(analysis){
    description <- analysis$description
    pt_ests <- list()
    pt_ests$controlY <- unname(description$controlY)
    pt_ests$interventionY <- unname(description$interventionY)
    pt_ests$effect_size <- unname(description$effect_size)
    pt_ests$contamination_interval <- NA
    pt_ests$personal_protection <- NA
    analysis$pt_ests <- pt_ests
    return(analysis)
}

Tanalysis <- function(analysis) {
    y1 <- arm <- cluster <- y_off <- NULL
    trial <- analysis$trial
    link <- analysis$options$link
    alpha <- analysis$options$alpha
    clusterSum <- data.frame(
        trial %>%
            group_by(cluster) %>%
            dplyr::summarize(
                y = sum(y1),
                total = sum(y_off),
                arm = arm[1]
            )
    )
    clusterSum$lp <- switch(link,
                            "identity" = clusterSum$y/clusterSum$total,
                            "log" = log(clusterSum$y/clusterSum$total),
                            "logit" = logit(clusterSum$y/clusterSum$total))
    formula <- stats::as.formula("lp ~ arm")

    # Trap any non-finite values

    clusterSum$lp[!is.finite(clusterSum$lp)] <- NA

    model_object <- stats::t.test(
        formula = formula, data = clusterSum, alternative = "two.sided",
        conf.level = 1 - alpha, var.equal = TRUE
    )
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
    if (link %in% c("logit","log")){
        analysis$pt_ests$effect_size <- 1 - analysis$pt_ests$interventionY/
                                            analysis$pt_ests$controlY
        analysis$int_ests$effect_size <- 1 - exp(-unlist(model_object$conf.int))
    }
    analysis$pt_ests$t.statistic <- model_object$statistic
    analysis$pt_ests$df <- unname(model_object$parameter)
    analysis$pt_ests$p.value <- model_object$p.value
    # tidy up and consolidate the list of analysis
    analysis$model_object <- model_object
    return(analysis)
}

extractEstimates <- function(analysis, sample) {
    alpha <- analysis$options$alpha
    link <- analysis$options$link
    method <- analysis$options$method
    CLnames <- analysis$options$CLnames
    sample$controlY <- invlink(link, sample$int)
    # personal_protection is the proportion of effect attributed to personal protection
    if ("arm" %in% names(sample) & "pvar" %in% names(sample)) {
        if (identical(method,"LME4")) sample$lc <- with(sample, int + pvar + arm)
        sample$interventionY <- invlink(link, sample$lc)
        sample$personal_protection <- with(
            sample, (controlY - invlink(link, int + arm))/(controlY -
                                                               interventionY))
    } else {
        sample$personal_protection <- NA
    }
    if ("arm" %in% names(sample) & !("pvar" %in% names(sample))) {
        sample$interventionY <- invlink(link, sample$int + sample$arm)
    }
    if ("pvar" %in% names(sample) & !("arms" %in% names(sample))) {
        sample$interventionY <- invlink(link, sample$int + sample$pvar)
    }
    if ("interventionY" %in% names(sample)) {
        sample$effect_size <- 1 - sample$interventionY/sample$controlY
    }
    if ("beta" %in% names(sample)) {
        sample$contamination_par <- sample$beta
        contamination_list <- apply(sample, MARGIN = 1, FUN = get_contamination, analysis = analysis)
        contamination_df <- as.data.frame(do.call(rbind, lapply(contamination_list, as.data.frame)))
        sample <- cbind(sample, contamination_df)
    }
    bounds <- (apply(
        sample, 2, function(x) {quantile(x, c(alpha/2, 0.5, 1 - alpha/2),
                                         alpha = alpha, na.rm = TRUE)}))
    analysis <- add_estimates(analysis = analysis, bounds = bounds, CLnames = CLnames)

return(analysis)
}


LME4analysis <- function(analysis, cfunc, trial, link, fterms){
    trial <- analysis$trial
    link <- analysis$options$link
    cfunc <- analysis$options$cfunc
    FUN <- get_FUN(cfunc, variant = 0)
    alpha <- analysis$options$alpha
    fterms <- analysis$options$fterms
    # TODO replace the use of ftext with fterms
    ftext <- analysis$options$ftext
    log_beta <- NA
    contrasts <- NULL
    fterms = switch(link,
                    "identity" = c("y1/y_off ~ 1", fterms),
                    "log" = c("y1 ~ 1", fterms, "offset(log(y_off))"),
                    "logit" = c("cbind(y1,y0) ~ 1", fterms))
    formula <- stats::as.formula(paste(fterms, collapse = "+"))
    if (cfunc %in% c("L", "P", "S")) {
        if (identical(Sys.getenv("TESTTHAT"), "true")) {
            log_beta <- 2.0
        } else {
            tryCatch({
            #cat("Estimating scale parameter for contamination range\n")
            log_beta <- stats::optimize(
                f = estimateContaminationLME4, interval = c(-5, 5), maximum = FALSE,
                tol = 0.1, trial = trial, FUN = FUN, formula = formula, link = link)$minimum
            },
            error = function(e){
               message("*** Contamination scale parameter cannot be estimated ***")
               log_beta <- 5
            })
        }
        beta <- exp(log_beta)
        x <- trial$nearestDiscord * beta
        trial$pvar <- eval(parse(text = FUN))
    }
    model_object <- switch(link,
           "identity" = lme4::lmer(formula = formula, data = trial, REML = FALSE),
           "log" = lme4::glmer(formula = formula, data = trial,
                               family = poisson),
           "logit" = lme4::glmer(formula = formula, data = trial,
                                 family = binomial))
    analysis$pt_ests$contamination_par <- exp(log_beta)
    analysis$pt_ests$deviance <- unname(summary(model_object)$AICtab["deviance"])
    analysis$pt_ests$AIC <- unname(summary(model_object)$AICtab["AIC"])
    analysis$pt_ests$df <- unname(summary(model_object)$AICtab["df.resid"])
    # if the contamination parameter has been estimated then penalise the AIC and
    # adjust the degrees of freedom in the output

    cov <- q50 <- NULL
    if (grepl("pvar", ftext, fixed = TRUE)) {
        analysis$pt_ests$AIC <- analysis$pt_ests$AIC + 2
        analysis$pt_ests$df <- analysis$pt_ests$df - 1
        coefficients <- summary(model_object)$coefficients
        WaldP <- ifelse(ncol(coefficients) == 4, coefficients['pvar',4], NA)
    }
    if (grepl("pvar", ftext, fixed = TRUE) |
        grepl("arm", ftext, fixed = TRUE)) {
        q50 <- summary(model_object)$coefficients[,1]
        names(q50)[grep("Int",names(q50))] <- "int"
        names(q50)[grep("arm",names(q50))] <- "arm"
        names(q50)[grep("pvar",names(q50))] <- "pvar"
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

# MCMC Methods
MCMCanalysis <- function(analysis){
    trial <- analysis$trial
    link <- analysis$options$link
    cfunc <- analysis$options$cfunc
    alpha <- analysis$options$alpha
    fterms <- analysis$options$fterms
    clusterEffects<- analysis$options$clusterEffects
    nchains <- 2
    max.iter <- 50000
    burnin <- 5000

    datajags <- list(N = nrow(trial))
    if (!identical(cfunc, "Z")) datajags$d <- trial$nearestDiscord
    if (identical(link, 'identity')) {
        datajags$y <- trial$y1/trial$y_off
    } else {
        datajags$y1 <- trial$y1
        datajags$y_off <- trial$y_off
    }
    if (clusterEffects) {
        datajags$cluster <- as.numeric(as.character(trial$cluster))
        datajags$ncluster <- max(as.numeric(as.character(trial$cluster)))
    }
    # construct the rjags code by concatenating strings

    text1 <- "model{\n
          for(i in 1:N){\n"

    text2 <- switch(cfunc, S = "pr[i] <- ifelse(d[i] < -beta/2,0,
                                        ifelse(d[i] > beta/2,1,
                                        (beta/2 + d[i])/beta))\n",
                    P = "pr[i] <- pnorm(d[i]*beta,0, 1) \n",
                    L = "pr[i] <- 1/(1 + exp(-beta*d[i])) \n",
                    X = "pr[i] <- ifelse(d[i] > 0, 1 ,0) \n",
                    Z = NULL)

    text3 <- switch(link,
                    "identity" = "y[i] ~ dnorm(lp[i],tau1) \n",
                    "log" =  "gamma1[i] ~ dnorm(0,tau1) \n
                                  Expect_y[i] <- exp(lp[i] + gamma1[i]) * y_off[i] \n
                                  y1[i] ~ dpois(Expect_y[i]) \n",
                    "logit" = "logitp[i] <- lp[i]  \n
                                   p[i] <- 1/(1 + exp(-logitp[i])) \n
                                   y1[i] ~ dbin(p[i],y_off[i]) \n"
    )

    # construct linear predictor
    text4 <- "lp[i] <- int + pvar * pr[i]"
    if (identical(cfunc,'Z')) text4 <- "lp[i] <- int"
    text5 <- ifelse(clusterEffects,
                    " + gamma[cluster[i]] \n
            }\n
            for(ic in 1:ncluster) {\n
                gamma[ic] ~ dnorm(0, tau)\n
            }\n
            tau <- 1/(sigma * sigma) \n
            sigma ~ dunif(0, 2) \n
            ", "}\n")
    text6 <-
        "log_beta ~ dnorm(0, 1E-1) \n
            beta <- exp(log_beta) \n
            int ~ dnorm(0, 1E-2) \n"
    text7 <- ifelse(identical(cfunc,'Z'), "pvar <- 0 \n", "pvar ~ dnorm(0, 1E-2) \n")
    text8 <- switch(link,
                    "identity" = "tau1 <- 1/(sigma1 * sigma1) \n
                                  sigma1 ~ dunif(0, 2) } \n",
                    "log" = "tau1 <- 1/(sigma1 * sigma1) \n
                             sigma1 ~ dunif(0, 2) } \n",
                    "logit" = "} \n"
    )

    MCMCmodel <- paste0(text1, text2, text3, text4, text5, text6, text7, text8)
    parameters.to.save <- switch(cfunc, S = c("int", "pvar", "beta"),
                    P = c("int", "pvar", "beta"),
                    L = c("int", "pvar", "beta"),
                    X = c("int", "pvar"),
                    Z = c("int"))

    model_object <- jagsUI::autojags(data = datajags, inits = NULL,
                       parameters.to.save = parameters.to.save,
                       model.file = textConnection(MCMCmodel), n.chains = nchains,
                       iter.increment = 1000, n.burnin = burnin, max.iter=max.iter)
    sample <- data.frame(rbind(model_object$samples[[1]],model_object$samples[[2]]))
    analysis$model_object <- model_object
    analysis <- extractEstimates(analysis = analysis, sample = sample)
    analysis$pt_ests$DIC <- model_object$DIC
    analysis$model_object$MCMCmodel <- MCMCmodel
return(analysis)
}


# Contributions to the linear predictor for different contamination functions

StraightLine <- function(par, trial)
{
    par[2] <- par[3] <- -9
    lp <- par[1]
    return(lp)
}

# step function for the case with no contamination

StepFunction <- function(par, trial)
{
    par[3] <- -9
    lp <- ifelse(trial$nearestDiscord < 0, par[1], par[1] + par[2])
    return(lp)
}



# piecewise linear model
PiecewiseLinearFunction <- function(par, trial)
{
    # constrain the slope parameter to be positive (par[2] is positive if effect_size is negative)
    beta <- par[3]

    # if beta is very large, the curve should be close to a straight line
    if (beta > 20){
        lp <- par[1] + 0.5 * par[2]

    } else {
        lp <- ifelse(
            trial$nearestDiscord > -beta/2, par[1] + par[2] * (beta/2 + trial$nearestDiscord)/beta, par[1])
        lp <- ifelse(trial$nearestDiscord > beta/2, par[1] + par[2], lp)
    }
    return(lp)
}

piecewise <- function(x) {
    value <- ifelse(x < -0.5, 0, (0.5 + x))
    value <- ifelse(x > 0.5, 1, value)
return(value)}



# sigmoid (logit) function
InverseLogisticFunction <- function(par, trial)
{
    lp <- par[1] + par[2] * invlink(link = "logit", x = par[3] * trial$nearestDiscord)
    return(lp)
}

# inverse probit function
InverseProbitFunction <- function(par, trial)
{
    lp <- par[1] + par[2] * stats::pnorm(par[3] * trial$nearestDiscord)
    return(lp)
}


get_FUN <- function(cfunc,variant){
    # TODO: remove the duplication and simplify here
    # Specify the function used for calculating the linear predictor
    if (variant == 1) {
        LPfunction <- c(
            "StraightLine", "StepFunction", "PiecewiseLinearFunction",
            "InverseLogisticFunction", "InverseProbitFunction")[which(cfunc == c("Z", "X", "S", "L", "P"))]
        FUN <- eval(parse(text = LPfunction))
    } else {
        # specify functional form of sigmoid in distance from boundary 'L' inverse logit; 'P' inverse probit; 'X'
        # or 'Z' do not model contamination
        FUN <- switch(
            cfunc, L = "invlink(link='logit', x)",
                 P = "stats::pnorm(x)",
                 S = "piecewise(x)",
                 X = NULL,
                 Z = NULL)
    }
    return(FUN)
}

get_contamination <- function(x, analysis){
        # define the limits of the curve both for control and intervention arms
    fittedCurve <- get_curve(x = x, analysis = analysis)
    contamination <- get_contaminationStats(fittedCurve=fittedCurve, trial=analysis$trial)
return(contamination)
}

get_curve <- function(x, analysis) {
    trial <- analysis$trial
    link <- analysis$options$link
    cfunc <- analysis$options$cfunc
    limits <- c(x[["controlY"]], x[["controlY"]])
    if(!is.null(x[["interventionY"]])) limits[2] <- x[["interventionY"]]
    if(is.na(limits[2])) limits[2] <- x[["controlY"]]
    personal_protection <- x[["personal_protection"]]
    contamination_par <- x[["contamination_par"]]

    # Trap cases with extreme effect: TODO: a different criterion may be needed for continuous data
    pars <- link_tr(link,limits)
    if (abs(pars[1] - pars[2]) > 10000) {
        limits <- invlink(link,c(20,-20))
    }

    limits0 <- limits1 <- limits
    range_d <- max(trial$nearestDiscord) -
        min(trial$nearestDiscord)
    d <- min(trial$nearestDiscord) + range_d * (seq(1:1001) - 1)/1000

    Cp <- 1
    if (is.na(personal_protection)) {
        Cp <- 0
    } else if (0 <= personal_protection & personal_protection <= 1) {
        Cp <- personal_protection
        limits0 <- c(limits[1], Cp * limits[2] + (1 - Cp) * limits[1])
        limits1 <- c(Cp * limits[1] + (1 - Cp) * limits[2], limits[2])
    }

    if (identical(limits[1], limits[2])) {
        curve <- rep(limits[1],1001)
    } else {
        par0 <- c(link_tr(link, limits0[1]),
                  link_tr(link, limits0[2]) - link_tr(link, limits0[1]),
                  contamination_par)
        par1 <- c(
            link_tr(link, limits1[1]),
            link_tr(link, limits1[2]) - link_tr(link, limits1[1]),
            contamination_par
        )
        # trap extreme cases with undefined, flat or very steep ones
        if (is.null(contamination_par)) {
            cfunc <- "X"
        } else if (is.na(contamination_par) |
            contamination_par < 0.01 |
            contamination_par > 100  |
            (abs(pars[1] - pars[2]) > 10000)) {
            cfunc <- "X"
        }
        FUN1 <- get_FUN(cfunc, variant = 1)

        curve <- ifelse(d < 0,
                        invlink(link, FUN1(trial = data.frame(nearestDiscord = d), par = par0)),
                        invlink(link, FUN1(trial = data.frame(nearestDiscord = d), par = par1)))
    }
    fittedCurve <- list(d = d, contaminationFunction = curve, limits0 = limits0, limits1 = limits1)
    return(fittedCurve)
}

# This is called once for each row in the sample data frame (for obtaining interval estimates)
get_contaminationStats <- function(fittedCurve, trial) {
    # Compute the contamination range
    # The absolute values of the limits are used so that a positive range is
    # obtained even with negative effect_size
    limits0 <- fittedCurve$limits0
    limits1 <- fittedCurve$limits1
    curve <- fittedCurve$contaminationFunction
    d <- fittedCurve$d
    thetaL <- thetaU <- NA
    if (is.na(limits0[1]) | is.na(curve[1000])) {
        rabbit <- 6
    }
    if (abs(limits0[1] - curve[1000]) >
        0.025 * abs(limits0[1] - limits0[2]))
    {
        thetaL <- d[min(
            which(
                abs(limits0[1] - curve) >
                    0.025 * abs(limits0[1] - limits0[2])
            )
        )]
    }
    if (abs(limits1[2] - curve[1000]) <
        0.025 * abs(limits1[1] - limits1[2]))
    {
        thetaU <- d[min(
            which(
                abs(limits1[2] - curve) <
                    0.025 * abs(limits1[1] - limits1[2])
            )
        )]
    }
    if (is.na(thetaU))
        thetaU <- max(trial$nearestDiscord)
    if (is.na(thetaL))
        thetaL <- min(trial$nearestDiscord)

    # contamination range
    contamination_limits <- c(thetaL, thetaU)
    if (thetaL > thetaU)
        contamination_limits <- c(thetaU, thetaL)

    contaminate_pop_pr <- sum(trial$nearestDiscord > contamination_limits[1] &
                                  trial$nearestDiscord < contamination_limits[2])/nrow(trial)
    contamination_interval <- thetaU - thetaL
    if (identical(thetaU, thetaL)) {
        contamination_interval <- 0
        # To remove warnings from plotting ensure that contamination interval is non-zero
        contamination_limits <- c(-1e-04, 1e-04)
    }

    contamination <- list(
        contamination_interval = contamination_interval,
        contamination_limit0 = contamination_limits[1],
        contamination_limit1 = contamination_limits[2],
        contaminate_pop_pr = contaminate_pop_pr)
return(contamination)}


tidyContamination <- function(contamination, analysis, fittedCurve){
    contamination$contamination_limits <-
        with(contamination, c(contamination_limit0,contamination_limit1))
    if (analysis$options$cfunc %in% c("Z","X")) {
        contamination$contamination_interval <- NULL
        contamination$contaminate_pop_pr <- NULL
        contamination$contamination_limits <- c(-1.0E-4,1.0E-4)
    } else {
        if (is.na(analysis$pt_ests$contamination_par))
            analysis$pt_ests$contamination_par <- contamination$contamination_par
        if (is.na(analysis$pt_ests$contamination_interval))
            analysis$pt_ests$contamination_interval <-
                contamination$contamination_interval
    }
    contamination$contamination_limit0 <- contamination$contamination_limit1 <- NULL
    contamination$FittedCurve <- data.frame(d = fittedCurve$d,
                        contaminationFunction = fittedCurve$contaminationFunction)
    contamination$data <- group_data(analysis)
    analysis$contamination <- contamination
    return(analysis)}
