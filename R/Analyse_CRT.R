#' Analysis of cluster randomized trial with contamination
#'
#' \\code{Analyse_CRT} returns outputs from a statistical analysis of a cluster randomized trial (CRT).
#' @param trial trial dataframe including locations, clusters, arms, and outcomes
#' @param method statistical method used to analyse trial.
#' options are:
#' 'EMP'  : empirical,
#' 'T'   : comparison of cluster means by t-test
#' 'ML'   : maximum likelihood,
#' 'GEE'  : generalised estimating equations
#' 'INLA' : INLA
#' @param cfunc transformation defining the contamination function
#' options are:
#' 'S': piecewise linear (slope),
#' 'L': inverse logistic (sigmoid),
#' 'P': inverse probit,
#' 'X': contamination not modelled'
#' @param link string: link function- options are 'logit' (the default), 'log', and 'identity'
#' @param numerator string: name of numerator variable for efficacy data (if present)
#' @param denominator string: name of denominator variable for efficacy data (if present)
#' @param excludeBuffer logical: indicator of whether any buffer zone (records with buffer=TRUE) should be excluded from analysis
#' @param alpha numeric: confidence level for confidence intervals and credible intervals
#' @param requireBootstrap logical: indicator of whether bootstrap confidence intervals are required
#' @param baselineOnly logical: indicator of whether required analysis is of efficacy or of baseline only
#' @param baselineNumerator string: name of numerator variable for baseline data (if present)
#' @param baselineDenominator string: name of denominator variable for baseline data (if present)
#' @param localisedEffects logical: indicator of whether the model includes local effects with no contamination
#' @param clusterEffects logical: indicator of whether the model includes cluster random effects
#' @param spatialEffects logical: indicator of whether the model includes spatial random effects
#' @param resamples integer: number of bootstrap samples
#' @param requireMesh logical: indicator of whether spatial predictions are required
#' @param inla.mesh name of pre-existing INLA input object created by createMesh()
#' @return list containing the following results of the analysis
#' \\itemize{
#' \\item \\code{description}: Description of the trial dataset
#' \\item \\code{method}: statistical method
#' \\item \\code{pt.ests}: point estimates
#' \\item \\code{int.ests}: interval estimates
#' \\item \\code{model.object}: object returned by the fitting function
#' \\item \\code{contamination}: function values and statistics describing the estimated contamination
#' }
#' @importFrom grDevices rainbow
#' @importFrom stats binomial dist kmeans median na.omit qlogis qnorm quantile rbinom rnorm runif simulate
#' @importFrom utils head read.csv
#' @export
#'
#' @examples
#' # Standard GEE analysis of test dataset ignoring contamination
#' exampleGEE=Analyse_CRT(trial=readdata('test_Simulate_CRT.csv'),method='GEE')

Analyse_CRT <- function(
    trial, method = "GEE", cfunc = "L", link = "logit", numerator = "num",
    denominator = "denom", excludeBuffer = FALSE, alpha = 0.05, requireBootstrap = FALSE,
    baselineOnly = FALSE, baselineNumerator = "base_num", baselineDenominator = "base_denom",
    localisedEffects = FALSE, clusterEffects = FALSE, spatialEffects = FALSE,
    resamples = 1000, requireMesh = TRUE, inla.mesh = NULL)
    {
    # Test of validity of inputs
    if (!method %in% c("EMP", "T", "ML", "GEE", "INLA"))
        {
        cat("Error: Invalid value for statistical method\n")
        return(NULL)
    }
    if (!cfunc %in% c("S", "L", "P", "X"))
        {
        cat("Error: Invalid contamination function\n")
        return(NULL)
    }

    # MAIN FUNCTION CODE STARTS HERE
    cat(
        "\n=====================    ANALYSIS OF CLUSTER RANDOMISED TRIAL    =================\n"
    )

    cluster <- NULL

    if ("buffer" %in% colnames(trial) &
        excludeBuffer)
        {
        trial <- trial[!trial$buffer, ]
    }

    # create names for confidence limits for use throughout
    CLnames <- c(
        paste0(alpha/0.02, "%"),
        paste0(100 - alpha/0.02, "%")
    )

    # trial needs to be ordered for some analyses
    trial <- trial[order(trial$cluster), ]

    if (method %in% c("EMP", "T", "GEE"))
        {
        cat("** Note: statistical method does not allow for contamination **\n")
        cfunc <- "X"
    }
    if (baselineOnly){
        if (method %in% c("EMP", "ML", "GEE"))
            {
            method <- "GEE"
            cat("Analysis of baseline only, using GEE\n")
        } else if (method == "INLA")
        {
            cat("Analysis of baseline only, using INLA\n")
        }
        # cfunc='Z' is used to remove the estimation of effect size from the model
        cfunc <- "Z"
        trial$y1 <- trial[[baselineNumerator]]
        trial$y0 <- trial[[baselineDenominator]] - trial[[baselineNumerator]]
        trial$y_off <- trial[[baselineDenominator]]

    } else {
        trial$y1 <- trial[[numerator]]
        trial$y0 <- trial[[denominator]] - trial[[numerator]]
        trial$y_off <- trial[[denominator]]

        # if nearestDiscord is not provided augment the trial data frame with distance to nearest discordant
        # coordinate (specifyBuffer assigns a buffer only if a buffer width is > 0 is input)
        if (is.null(trial$nearestDiscord)){
            trial <- Specify_CRTbuffer(trial = trial, bufferWidth = 0)
        }
    }
    model.object <- list()
    pt.ests <- list(contamination.par = NA, pr.contaminated = NA, contaminationRange = NA)
    int.ests <- list(controlY = NA, interventionY = NA, efficacy = NA)
    sd <- 0.5/(qnorm(1 - alpha) * sqrt(2))  #initial value used in bootstrap calculations
    description <- ifelse(baselineOnly, list(), get_description(trial))
    # Specify the function used for calculating the linear predictor
    LPfunction <- c(
        "CalculateNoEffect", "CalculateNoContaminationFunction", "CalculatePiecewiseLinearFunction",
        "CalculateLogisticFunction", "CalculateProbitFunction")[which(cfunc == c("Z", "X", "S", "L", "P"))]
    FUN2 <- FUN1 <- eval(parse(text = LPfunction))
    # a second 'link' variable is required for the negative binomial case

    if (method == "EMP"){
        description <- get_description(trial)
        fit <- list(
            controlY = unname(description$ratios[1]),
            interventionY = unname(description$ratios[2]),
            efficacy = unname(description$efficacy),
            contaminationRange = NA)
        pt.ests$controlY <- fit$controlY
        pt.ests$interventionY <- fit$interventionY
        pt.ests$efficacy <- fit$efficacy
        if (requireBootstrap){
            boot_emp <- boot::boot(
                data = trial, statistic = BootEmpiricalAnalysis, R = resamples,
                sim = "parametric", ran.gen = rgen_emp, mle = fit
            )
            pt.ests$bootstrapMean_efficacy <- mean(boot_emp$t)
            int.ests$efficacy <- namedCL(
                quantile(boot_emp$t, c(alpha/2, 1 - alpha/2)),
                alpha = alpha
            )
        }
    } else if (method == "T"){
        y1 <- arm <- NULL
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
        model.object <- stats::t.test(
            formula = formula, data = clusterSum, alternative = "two.sided",
            conf.level = 1 - alpha, var.equal = TRUE
        )
        pt.ests$p.value <- model.object$p.value
        analysisC <- stats::t.test(
            clusterSum$lp[clusterSum$arm == "control"], conf.level = 1 - alpha)
        pt.ests$controlY <- invlink(link, analysisC$estimate[1])
        int.ests$controlY <- invlink(link, analysisC$conf.int)
        analysisI <- stats::t.test(
            clusterSum$lp[clusterSum$arm == "intervention"], conf.level = 1 - alpha)
        pt.ests$interventionY <- invlink(link, analysisI$estimate[1])
        int.ests$interventionY <- invlink(link, analysisI$conf.int)

        # Covariance matrix (note that two arms are independent so the off-diagonal elements are zero)
        Sigma <- matrix(
            data = c(analysisC$stderr^2, 0, 0, analysisI$stderr^2),
            nrow = 2, ncol = 2)
        pt.ests$efficacy <- 1 - pt.ests$interventionY/pt.ests$controlY
        int.ests$efficacy <- estimateCLEfficacy(
            mu = c(analysisC$estimate, analysisI$estimate),
            Sigma = Sigma, alpha = alpha, resamples = resamples, method = method,
            link = link)
    } else if (method == "GEE") {
        # GEE analysis of cluster effects

        if (link == "log")
        {
            fterms <- ifelse(
                cfunc == "Z", "y1 ~ 1 + offset(log(y_off))", "y1 ~ arm + offset(log(y_off))")
            formula <- stats::as.formula(fterms)
            fit <- geepack::geeglm(
                formula = formula, id = cluster, data = trial, family = poisson(link = "log"),
                corstr = "exchangeable", scale.fix = FALSE)
        } else if (link == "logit") {
            fterms <- ifelse(cfunc == "Z", "cbind(y1,y0) ~ 1", "cbind(y1,y0) ~ arm")
            formula <- stats::as.formula(fterms)
            fit <- geepack::geeglm(
                formula = formula, id = cluster, corstr = "exchangeable",
                data = trial, family = binomial(link = "logit"))
        } else if (link == "identity") {
            fterms <- ifelse(cfunc == "Z", "y1/y_off ~ 1", "y1/y_off ~ arm")
            formula <- stats::as.formula(fterms)
            fit <- geepack::geeglm(
                formula = formula, id = cluster, corstr = "exchangeable",
                data = trial, family = gaussian
            )
        }
        summary_fit <- summary(fit)

        z <- -qnorm(alpha/2)  #standard deviation score for calculating confidence intervals
        lp_yC <- summary_fit$coefficients[1, 1]
        se_lp_yC <- summary_fit$coefficients[1, 2]

        clusterSize <- nrow(trial)/nlevels(as.factor(trial$cluster))


        # remove the temporary objects from the dataframe
        fit$model.object$data$y1 <- fit$model.object$data$y0 <- fit$model.object$data$y_off <- NULL
        pt.ests$controlY <- invlink(link, lp_yC)
        int.ests$controlY <- namedCL(
            invlink(link, c(lp_yC - z * se_lp_yC, lp_yC + z * se_lp_yC)),
            alpha = alpha
        )

        pt.ests$model.object <- fit

        # Intracluster correlation
        pt.ests$ICC <- noLabels(summary_fit$corr[1])  #with corstr = 'exchangeable', alpha is the ICC
        se_ICC <- noLabels(summary_fit$corr[2])
        int.ests$ICC <- namedCL(
            noLabels(c(pt.ests$ICC - z * se_ICC, pt.ests$ICC + z * se_ICC)),
            alpha = alpha
        )
        pt.ests$DesignEffect <- 1 + (clusterSize - 1) * pt.ests$ICC  #Design Effect
        int.ests$DesignEffect <- 1 + (clusterSize - 1) * int.ests$ICC

        # Estimation of efficacy does not apply if analysis is of baseline only (cfunc='Z')
        pt.ests$efficacy <- 0
        if (cfunc == "X")
        {
            lp_yI <- summary_fit$coefficients[1, 1] + summary_fit$coefficients[2,
                1]
            se_lp_yI <- sqrt(
                fit$geese$vbeta[1, 1] + fit$geese$vbeta[2, 2] + 2 * fit$geese$vbeta[1,
                  2]
            )

            int.ests$interventionY <- namedCL(
                invlink(link, c(lp_yI - z * se_lp_yI, lp_yI + z * se_lp_yI)),
                alpha = alpha
            )
            int.ests$efficacy <- estimateCLEfficacy(
                mu = summary_fit$coefficients[, 1], Sigma = fit$geese$vbeta,
                alpha = alpha, resamples = resamples, method = method,
                link = link
            )

            pt.ests$interventionY <- invlink(link, lp_yI)
            pt.ests$efficacy <- (1 - invlink(link, lp_yI)/invlink(link, lp_yC))
        }

    } else if (method == "ML")
    {
        # ML Methods with contamination functions and logistic link

        par <- SingleTrialAnalysis(trial = trial, FUN2 = FUN2)
        fit <- FittingResults(trial, par = par, FUN1 = FUN1)
        pt.ests$controlY <- fit$controlY
        pt.ests$interventionY <- fit$interventionY
        pt.ests$efficacy <- fit$efficacy
        pt.ests$contamination.par <- fit$contamination.par
        pt.ests <- pt.ests[names(pt.ests) != "model.object"]
        pt.ests$pr.contaminated <- 1  #None of the ML models include local effects
        results <- list(
            description = description, method = method, pt.ests = pt.ests,
            int.ests = int.ests, model.object = model.object
        )
        if (requireBootstrap)
        {
            mle <- list(par = par, FUN1 = FUN1, link = "logit")
            boot_estimates <- data.frame(V1 = c(), V2 = c(), V3 = c(), V4 = c())
            # resampling can crash because of resamples containing data from only one arm to prevent this
            # crashing the whole program, bootstrapping is done in batches of 5 resamples using 'try' to avoid
            # crashing out
            resamples1 <- 5
            tries <- 0
            while (nrow(boot_estimates) <
                resamples)
                {
                # sample a random value each time through the loop so the seed is change
                boot_output <- try(
                  expr = boot::boot(
                    data = trial, statistic = SingleTrialAnalysis, R = resamples1,
                    sim = "parametric", ran.gen = rgen, mle = mle, FUN2 = FUN2
                ),
                  silent = TRUE
              )
                if (!substr(boot_output[1], 1, 5) ==
                  "Error")
                  {
                  new_estimates <- as.data.frame(
                    t(
                      matrix(
                        data = unlist(
                          apply(
                            as.data.frame(boot_output$t),
                            1, FittingResults, trial = trial, FUN1 = FUN1
                        )
                      ),
                        ncol = resamples1, nrow = length(fit)
                    )
                  )
                )

                  boot_estimates <- rbind(boot_estimates, new_estimates)
                }
                tries <- tries + 5
                cat(
                  "\r", nrow(boot_estimates),
                  " bootstrap samples analysed, out of", tries, " tries    \r"
              )
            }
            colnames(boot_estimates) <- names(fit)
            bounds <- (apply(
                boot_estimates, 2, function(x)
                  {
                  quantile(
                    x, c(alpha/2, 0.5, 1 - alpha/2),
                    alpha = alpha
                )
                }
            ))
            results <- add_estimates(results = results, bounds = bounds, CLnames = CLnames)
        }


        ############### INLA Methods #################
    } else if (method == "INLA")
    {
        trial <- dplyr::mutate(trial, id =  dplyr::row_number())
        if (is.null(inla.mesh))
            {
            inla.mesh <- createMesh(
                trial = trial, offset = -0.1, max.edge = 0.25, inla.alpha = 2,
                maskbuffer = 0.5, ncells = 50
            )
        }

        y_off <- NULL
        # specify functional form of sigmoid in distance from boundary 'L' inverse logit; 'P' inverse probit; 'X'
        # or 'Z' do not model contamination
        FUN <- switch(
            cfunc, L = "invlink(link='logit', x)", P = "stats::pnorm(x)", X = NULL, Z = NULL)

        # create model formula
        fterms <- switch(link,
            "identity" = "y1/y_off ~ 0",
            "log" = "y1 ~ 0",
            "logit" = "y1 ~ 0"
        )

        fterms <- c(fterms, switch(
            cfunc, Z = "b0",
            X = "b0 + b1",
            L = "b0 + pvar",
            P = "b0 + pvar"
        ))

        if (localisedEffects & cfunc != 'X')
            fterms <- c(fterms, "b1")
        if (clusterEffects)
            fterms <- c(fterms, "f(cluster, model = \"iid\")")
        if (spatialEffects)
            fterms <- c(fterms, "f(s, model = spde)")
        if (link == 'log')
            fterms <- c(fterms, "f(id, model = \"iid\")")

        # console display of the formula
        formula_as_text <- paste(fterms, collapse = " + ")
        cat(formula_as_text, "\n")
        formula <- stats::as.formula(formula_as_text)

        spde <- inla.mesh$spde

        effectse <- list(
            df = data.frame(
                b0 = rep(1, nrow(trial)),
                b1 = ifelse(trial$arm == "intervention", 1, 0),
                id = trial$id,
                cluster = trial$cluster
            ),
            s = inla.mesh$indexs
        )
        effectsp <- list(
            df = data.frame(
                b0 = rep(1, nrow(inla.mesh$prediction)),
                b1 = ifelse(inla.mesh$prediction$arm == "intervention", 1, 0),
                id = inla.mesh$prediction$id,
                cluster = inla.mesh$prediction$cluster
            ),
            s = inla.mesh$indexs
        )
        lc <- NULL
        beta2 <- NA
        if (cfunc %in% c("L", "P"))
            {
            cat("Estimating scale parameter for contamination range", "\n")
            beta2 <- stats::optimize(
                f = estimateContamination, interval = c(-10, 10),
                trial = trial, FUN = FUN, inla.mesh = inla.mesh, formula = formula,
                tol = 0.1, link = link)$minimum
            x <- trial$nearestDiscord * exp(beta2)
            trial$pvar <- eval(parse(text = FUN))

            effectse$df$pvar <- trial$pvar
            x <- inla.mesh$prediction$nearestDiscord * exp(beta2)
            inla.mesh$prediction$pvar <- ifelse(
                cfunc == "X", rep(NA, nrow(inla.mesh$prediction)),
                eval(parse(text = FUN))
            )
            effectsp$df$pvar <- inla.mesh$prediction$pvar
            # set up linear contrasts (not required for cfunc='X' or 'Z')
            if (grepl("pvar", formula_as_text, fixed = TRUE))
                {
                lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1)
                if (grepl("b1", formula_as_text, fixed = TRUE))
                  {
                  lc <- INLA::inla.make.lincomb(b0 = 1, pvar = 1, b1 = 1)
                }
            } else if (grepl("b1", formula_as_text, fixed = TRUE))
                {
                lc <- INLA::inla.make.lincomb(b0 = 1, b1 = 1)
            }
        }
        # stack for estimation stk.e
        stk.e <- INLA::inla.stack(
            tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
            A = list(1, A = inla.mesh$A),
            effects = effectse
        )

        # stack for prediction stk.p
        stk.p <- INLA::inla.stack(
            tag = "pred", data = list(y1 = NA, y_off = NA),
            A = list(1, inla.mesh$Ap),
            effects = effectsp
        )

        # stk.full comprises both stk.e and stk.p
        stk.full <- INLA::inla.stack(stk.e, stk.p)
        cat("INLA analysis                                                 \n")

        if (link == "identity") {
            inla.result <- INLA::inla(
                formula, family = "gaussian", lincomb = lc,
                control.family = list(link = "identity"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        } else if (link == "log") {
            inla.result <- INLA::inla(
                formula, family = "poisson", lincomb = lc,
                control.family = list(link = "log"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        } else if (link == "logit") {
            inla.result <- INLA::inla(
                formula, family = "binomial", Ntrials = y_off, lincomb = lc,
                control.family = list(link = "logit"),
                data = INLA::inla.stack.data(stk.full),
                control.fixed = list(correlation.matrix = TRUE),
                control.predictor = list(compute = TRUE, link = 1,
                                         A = INLA::inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE))
        }

        # Augment the inla results list with application specific quantities
        index <- INLA::inla.stack.index(stack = stk.full, tag = "pred")$data
        inla.mesh$prediction$prediction <- invlink(link, inla.result$summary.linear.predictor[index, "0.5quant"])

        # Compute sample-based confidence limits for intervened outcome and efficacy if intervention effects are
        # estimated
        if (grepl("pvar", formula_as_text, fixed = TRUE) |
            grepl("b1", formula_as_text, fixed = TRUE))
                {
            # Specify the means of the variables
            mu <- inla.result$summary.lincomb.derived$mean
            names(mu) <- rownames(inla.result$summary.lincomb.derived)
            # Specify the covariance matrix of the variables
            cov <- inla.result$misc$lincomb.derived.covariance.matrix
            sample <- as.data.frame(MASS::mvrnorm(n = 10000, mu = mu, Sigma = cov))
            sample$controlY <- invlink(link, sample$b0)
            # pr.contaminated is the proportion of effect subject to contamination
            if ("b1" %in% names(mu) &
                "pvar" %in% names(mu))
                  {
                sample$interventionY <- invlink(link, sample$lc)
                sample$pr.contaminated <- with(
                  sample, 1 - (controlY - invlink(link, b0 + b1))/(controlY -
                    interventionY)
              )
            } else if ("b1" %in% names(mu))
                {
                sample$interventionY <- invlink(link, sample$b0 + sample$b1)
                sample$pr.contaminated <- 0
            } else if ("pvar" %in% names(mu))
                {
                sample$interventionY <- invlink(link, sample$b0 + sample$pvar)
                sample$pr.contaminated <- 1
            }
            sample$efficacy <- 1 - sample$interventionY/sample$controlY
            bounds <- (apply(
                sample, 2, function(x)
                  {
                  quantile(
                    x, c(alpha/2, 0.5, 1 - alpha/2),
                    alpha = alpha
                )
                }
            ))
        } else
        {
            controlY <- unlist(
                invlink(
                  link, inla.result$summary.fixed["b0", c("0.025quant", "0.5quant", "0.975quant")]
              )
            )
            bounds <- data.frame(
                controlY = controlY, interventionY = controlY, efficacy = rep(0, 3),
                pr.contaminated = rep(0, 3)
            )
        }
        results <- list(
            model.object = inla.result, inla.mesh = inla.mesh, pt.ests = list(),
            int.ests = list(), method = method
        )

        results <- add_estimates(results = results, bounds = bounds, CLnames = CLnames)
        results$int.ests$pr.contaminated <- stats::setNames(
            bounds[c(1, 3),
                "pr.contaminated"], CLnames
        )
        results$pt.ests$pr.contaminated <- bounds[2, "pr.contaminated"]

        results$passed.face.validity.check <- TRUE
        if (results$pt.ests$pr.contaminated < 0 | results$pt.ests$pr.contaminated >
            1)
            {
            cat(
                "** Warning: different signs for main effect and contamination: face validity check fails **\n"
            )
            results$passed.face.validity.check <- FALSE
            results$pt.ests$pr.contaminated <- NA
            results$int.ests$pr.contaminated <- c(NA, NA)
        }
        results$pt.ests$contamination.par <- beta2
        results$description <- description
    }
    if (method %in% c("EMP", "ML", "T", "GEE"))
        {
        # tidy up and consolidate the list of results
        model.object <- pt.ests$model.object
        pt.ests <- pt.ests[names(pt.ests) !=
            "GAsolution3"]
        pt.ests <- pt.ests[names(pt.ests) !=
            "model.object"]
        results <- list(
            description = description, method = method, pt.ests = pt.ests,
            int.ests = int.ests, model.object = model.object
        )
    }
    if (cfunc != "Z")
    {
        results$contamination <- getContaminationCurve(trial = trial, pt.ests = results$pt.ests,
                                                       FUN1 = FUN1, link = link, alpha = alpha)
        results$pt.ests$contaminationRange <- results$contamination$contaminationRange
        results$contamination$contaminationRange <- NULL
    } else
    {
        results$pt.ests$contaminationRange <- NA
    }
    ## Output to screen

    cat(
        "Analysis model: ", method, "Link function: ", link, "Contamination option: ",
        cfunc, "\n"
    )
    CLtext <- paste0(" (", 100 * (1 - alpha), "% CL: ")
    cat(
        "Estimate-      Control: ", results$pt.ests$controlY,
        CLtext, unlist(results$int.ests$controlY),
        ")\n"
    )
    if (results$pt.ests$efficacy != 0)
    {
        if (!is.null(results$pt.ests$interventionY))
            {
            cat(
                "          Intervention: ", results$pt.ests$interventionY,
                CLtext, unlist(results$int.ests$interventionY),
                ")\n"
            )
            cat(
                "Efficacy: ", results$pt.ests$efficacy, CLtext, unlist(results$int.ests$efficacy),
                ")\n"
            )
        }
        if (!is.na(results$pt.ests$pr.contaminated))
            {
            cat(
                "Proportion of effect subject to contamination: ", results$pt.ests$pr.contaminated,
                CLtext, unlist(results$int.ests$pr.contaminated),
                ")\n"
            )
        }
    }
    if (!is.null(pt.ests$ICC))
        {
        cat(
            "Estimated intracluster correlation (ICC): ", results$pt.ests$ICC,
            CLtext, unlist(results$int.ests$ICC),
            ")\n"
        )
    }
    if (!is.na(results$pt.ests$contaminationRange))
        {
        cat(
            "Contamination Range: ", results$pt.ests$contaminationRange,
            "\n"
        )
    }
    if (!is.null(results$model.object$dic$dic) &
        cfunc %in% c("L", "P"))
            {
        # The contamination parameter is not estimated by INLA but should be considered in the DIC
        results$model.object$dic$dic <- results$model.object$dic$dic + 2
        cat(
            "DIC: ", results$model.object$dic$dic, " including penalty for the contamination parameter\n"
        )
    } else if (!is.null(results$model.object$dic$dic))
        {
        cat("DIC: ", results$model.object$dic$dic, "\n")
    }
    return(results)
}

# Add estimates to results list
add_estimates <- function(results, bounds, CLnames)
    {
    results$pt.ests$controlY <- bounds[2, "controlY"]
    results$pt.ests$interventionY <- bounds[2, "interventionY"]
    results$pt.ests$efficacy <- bounds[2, "efficacy"]

    # Extract interval estimates
    results$int.ests$controlY <- stats::setNames(
        bounds[c(1, 3),
            "controlY"], CLnames
    )
    results$int.ests$interventionY <- stats::setNames(
        bounds[c(1, 3),
            "interventionY"], CLnames
    )
    results$int.ests$efficacy <- stats::setNames(
        bounds[c(1, 3),
            "efficacy"], CLnames
    )
    return(results)
}

getContaminationCurve <- function(trial, pt.ests, FUN1, link, alpha)
    {

    y_off <- y1 <- average <- upper <- lower <- cats <- nearestDiscord <- NULL

    range_d <- max(trial$nearestDiscord) -
        min(trial$nearestDiscord)
    d <- min(trial$nearestDiscord) +
        range_d * (seq(1:1001) -
            1)/1000

    # define the limits of the curve both for control and intervention arms
    limits <- c(pt.ests$controlY, pt.ests$interventionY)
    limits0 <- limits1 <- limits
    Cp <- 1
    if (is.na(pt.ests$pr.contaminated))
        {
        Cp <- 0
    } else if (0 <= pt.ests$pr.contaminated & pt.ests$pr.contaminated <=
        1)
        {
        Cp <- pt.ests$pr.contaminated
        limits0 <- c(limits[1], Cp * limits[2] + (1 - Cp) * limits[1])
        limits1 <- c(Cp * limits[1] + (1 - Cp) * limits[2], limits[2])
    }

    par0 <- c(link_tr(link, limits0[1]),
              link_tr(link, limits0[2]) - link_tr(link, limits0[1]),
                pt.ests$contamination.par
    )
    par1 <- c(
        link_tr(link, limits1[1]),
        link_tr(link, limits1[2]) - link_tr(link, limits1[1]),
        pt.ests$contamination.par
    )
    curve <- ifelse(d < 0,
               invlink(link, FUN1(trial = data.frame(nearestDiscord = d), par = par0)),
               invlink(link, FUN1(trial = data.frame(nearestDiscord = d), par = par1)))

    # estimate contamination range The absolute values of the limits are used so that a positive range is
    # obtained even with negative efficacy
    thetaL <- thetaU <- NA
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
    contaminatedInterval <- c(thetaL, thetaU)
    if (thetaL > thetaU)
        contaminatedInterval <- c(thetaU, thetaL)
    contaminationRange <- thetaU - thetaL
    if (Cp == 0)
        contaminationRange <- NA
    # To remove warnings from plotting ensure that contamination interval is non-zero
    if (is.na(contaminationRange) ||
        contaminationRange == 0)
        {
        contaminatedInterval <- c(-1e-04, 1e-04)
    }
    # categorisation of trial data for plotting
    trial$cats <- cut(
        trial$nearestDiscord, breaks = c(
            -Inf, min(d) +
                seq(1:9) *
                  range_d/10, Inf
        ),
        labels = FALSE
    )

    data <- data.frame(
        trial %>%
            group_by(cats) %>%
            dplyr::summarize(
                positives = sum(y1),
                total = sum(y_off),
                d = median(nearestDiscord),
                average = Williams(x=y1/y_off, alpha=alpha, option = 'M'),
                lower = Williams(x=y1/y_off, alpha=alpha, option = 'L'),
                upper = Williams(x=y1/y_off, alpha=alpha, option = 'U')
            )
    )
    if (link == 'logit') {
        # proportions and binomial confidence intervals by category
        data$average <- data$positives/data$total
        data$upper <- with(data, average -
                               qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
        data$lower <- with(data, average +
                               qnorm(alpha/2) * (sqrt(average * (1 - average)/total)))
    }

    returnList <- list(
        FittedCurve = data.frame(d = d, contaminationFunction = curve),
        contaminationRange = contaminationRange, contaminatedInterval = contaminatedInterval,
        data = data
    )
    return(returnList)
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

# Minimal data description and crude efficacy estimate
get_description <- function(trial)
    {
    positives <- tapply(trial$y1, trial$arm, FUN = sum)
    totals <- tapply(trial$y_off, trial$arm, FUN = sum)
    ratios <- positives/totals
    efficacy <- 1 - ratios[2]/ratios[1]
    description <- list(
        positives = positives, totals = totals, ratios = ratios, efficacy = efficacy,
        method = "DE"
    )
    return(description)
}

# Log Likelihood to be maximized
LogLikelihood <- function(par, FUN = FUN, trial)
    {
    logitexpectP <- FUN(par, trial)
    transf <- 1/(1 + exp(-logitexpectP))  #inverse logit transformation

    # for binomial
    LogLikelihood <- sum(
        trial$y1 * log(transf) +
            trial$y0 * log(1 - transf)
    )

    return(LogLikelihood)
}


# transform the parameters into interpretable functions
FittingResults <- function(trial, FUN1, par)
    {
    controlY <- invlink(link = "logit", x = par[1])
    interventionY <- invlink(link = "logit", x = (par[2] + par[1]))
    efficacy <- (controlY - interventionY)/controlY
    pt.ests <- list(
        controlY = controlY, interventionY = interventionY, efficacy = efficacy,
        contamination.par = par[3]
    )
    return(pt.ests)
}

###########################################################################
# Contributions to the linear predictor for different contamination functions

CalculateNoEffect <- function(par, trial)
    {
    lp <- par[1]
    return(lp)
}

# step function for the case with no contamination

CalculateNoContaminationFunction <- function(par, trial)
    {
    lp <- ifelse(trial$nearestDiscord < 0, par[1], par[1] + par[2])
    return(lp)
}


# piecewise linear model (on the logit scale)
CalculatePiecewiseLinearFunction <- function(par, trial)
    {
    # constrain the slope parameter to be positive (par[2] is positive if efficacy is negative)
    theta <- exp(par[3])
    lp <- ifelse(
        trial$nearestDiscord > -theta, par[1] + par[2] * (theta + trial$nearestDiscord)/(2 *
            theta), par[1]
    )
    lp <- ifelse(trial$nearestDiscord > theta, par[1] + par[2], lp)
    return(lp)
}


# sigmoid (logit) function (on the logit scale)
CalculateLogisticFunction <- function(par, trial)
    {
    theta <- exp(par[3])
    lp <- par[1] + par[2] * invlink(link = "logit", x = theta * trial$nearestDiscord)
    return(lp)
}

# inverse probit function (on the logit scale)
CalculateProbitFunction <- function(par, trial)
    {
    theta <- exp(par[3])
    lp <- par[1] + par[2] * stats::pnorm(theta * trial$nearestDiscord)
    return(lp)
}


# Functions for T and GEE analysis

noLabels <- function(x)
    {
    xclean <- as.matrix(x)
    dimnames(xclean) <- NULL
    xclean <- as.vector(xclean)
    return(xclean)
}

estimateCLEfficacy <- function(mu, Sigma, alpha, resamples, method, link)
    {

    # Use resampling approach to avoid need for Taylor approximation use at least 10000 samples (this is very
    # cheap)
    resamples1 <- max(resamples, 10000, na.rm = TRUE)
    samples <- MASS::mvrnorm(n = resamples1, mu = mu, Sigma = Sigma)

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

# functions for analysis of Maximum Likelihood models

rgen <- function(data, mle, link)
    {
    par <- mle$par
    FUN1 <- mle$FUN1
    # simulate data for numerator y1
    modelp <- FUN1(par = par, trial = data)
    transf <- invlink(link = mle$link, x = modelp)
    data$y1 <- rbinom(
        length(transf),
        data$y_off, transf
    )  #simulate from binomial distribution
    return(data)
}

SingleTrialAnalysis <- function(trial, FUN2 = FUN2)
    {

    GA <- GA::ga(
        "real-valued", fitness = LogLikelihood, FUN = FUN2, trial = trial,
        lower = c(-10, -10, -10),
        upper = c(10, 10, 10),
        maxiter = 500, run = 50, optim = TRUE, monitor = FALSE
    )
    result <- GA@solution

    return(result)
}



############################################################################## functions for Empirical Analysis

rgen_emp <- function(data, mle)
    {

    description <- psych::describeBy(data$y1/data$y_off, group = data$arm)
    pChat <- description$control$mean
    pIhat <- description$intervention$mean

    # simulate data for numerator num
    modelp <- ifelse(
        as.numeric(data$arm) -
            1 > 0, pIhat, pChat
    )
    data$y1 <- rbinom(
        length(modelp),
        data$y_off, modelp
    )  #simulate from Binomial distribution

    return(data)
}

# standard non-model based analyses
BootEmpiricalAnalysis <- function(resampledData)
    {
    description <- psych::describeBy(resampledData$y1/resampledData$y_off, group = resampledData$arm)
    # reports summary statistic by a grouping variable
    pChat <- description$control$mean
    pIhat <- description$intervention$mean
    Eshat <- 1 - pIhat/pChat

    return(Eshat)
}

########## FUNCTIONS FOR SPATIAL PARTIAL DIFFERENTIAL EQUATION MODEL IN INLA

#' \\code{createMesh} Create prediction mesh and other inputs required for INLA analyis of a CRT.
#' @param trial trial dataframe including locations, clusters, arms, and binary outcomes
#' @param offset (see inla.mesh.2d documentation)
#' @param max.edge (see inla.mesh.2d documentation)
#' @param inla.alpha parameter related to the smoothness
#' @param maskbuffer (see inla.mesh.2d documentation)
#' @param ncells resolution of mesh in terms of maximum of linear dimension in pixels
#' @return list containing the mesh
#' \\itemize{
#' \\item \\code{prediction}: Data.table containing the prediction points and covariate values
#' \\item \\code{A}: projection matrix from the observations to the mesh nodes.
#' \\item \\code{Ap}: projection matrix from the prediction points to the mesh nodes.
#' \\item \\code{indexs}:  index set for the SPDE model
#' \\item \\code{spde}: SPDE model
#' }
#' @export
#'
#' @examples
#' # low resolution mesh for test dataset
#' exampleMesh=createMesh(trial = readdata('test_Simulate_CRT.csv'), ncells = 7)
createMesh <- function(
    trial = trial, offset = -0.1, max.edge = 0.25, inla.alpha = 2, maskbuffer = 0.5,
    ncells = 50)
    {
    cat(
        "Creating mesh for INLA analysis: resolution parameter= ", ncells,
        "\n"
    )
    # create an id variable if this does not exist
    if(is.null(trial$id)) trial <- dplyr::mutate(trial, id =  dplyr::row_number())

    # create buffer around area of points
    trial.coords <- base::matrix(
        c(trial$x, trial$y),
        ncol = 2
    )
    sptrial <- sp::SpatialPoints(trial.coords)
    buf1 <- rgeos::gBuffer(sptrial, width = maskbuffer, byid = TRUE)
    buffer <- rgeos::gUnaryUnion(buf1)

    # estimation mesh construction

    mesh <- INLA::inla.mesh.2d(
        boundary = buffer, offset = offset, cutoff = 0.05, max.edge = max.edge
    )

    # set up SPDE (Stochastic Partial Differential Equation) model
    spde <- INLA::inla.spde2.matern(mesh = mesh, alpha = inla.alpha, constr = TRUE)
    indexs <- INLA::inla.spde.make.index("s", spde$n.spde)
    A <- INLA::inla.spde.make.A(mesh = mesh, loc = trial.coords)

    # 8.3.6 Prediction data from https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html
    bb <- sp::bbox(buffer)
    x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = ncells)
    y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = ncells)
    pred.coords <- as.matrix(expand.grid(x, y))
    buf.coords <- buffer@polygons[[1]]@Polygons[[1]]@coords
    ind <- sp::point.in.polygon(
        pred.coords[, 1], pred.coords[, 2], buf.coords[, 1], buf.coords[,
            2]
    )
    # prediction locations
    pred.coords <- pred.coords[which(ind == 1),
        ]

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
    distM <- matrix(
        distP, nrow = nrow(pred.coords),
        ncol = nrow(trial),
        byrow = TRUE
    )
    nearestNeighbour <- apply(distM, 1, function(x) return(array(which.min(x))))
    armP <- trial$arm[nearestNeighbour]
    clusterP <- trial$cluster[nearestNeighbour]
    idP <- trial$id[nearestNeighbour]
    prediction <- data.frame(
        x = pred.coords[, 1], y = pred.coords[, 2], nearestNeighbour = nearestNeighbour,
        arm = armP, id = idP, cluster = clusterP
    )
    prediction <- with(
        prediction, prediction[order(y, x),
            ]
    )
    prediction$shortestDistance <- apply(distM, 1, min)
    rows <- seq(1:nrow(prediction))
    prediction$nearestDiscord <- sapply(rows,
                                        FUN = calcNearestDiscord,
                                        trial = trial,
                                        prediction = prediction,
                                        distM = distM)
    inla.mesh <- list(
        prediction = prediction, A = A, Ap = Ap, indexs = indexs, spde = spde
    )
    return(inla.mesh)
}

# Calculate the distance to the nearest discordant location
calcNearestDiscord <- function(x, trial , prediction , distM)
{
    discords <- (trial$arm != prediction$arm[x])
    nearestDiscord <- min(distM[x, discords])
    nearestDiscord <- ifelse(prediction$arm[x] == "control", -nearestDiscord, nearestDiscord)
    return(nearestDiscord)
}

# Use profiling to estimate beta2
estimateContamination <- function(
    beta2 = beta2, trial = trial, FUN = FUN, inla.mesh = inla.mesh, formula = formula, link = link){
    y1 <- y0 <- y_off <- NULL
    x <- trial$nearestDiscord * exp(beta2)
    trial$pvar <- eval(parse(text = FUN))

    stk.e <- INLA::inla.stack(
        tag = "est", data = list(y1 = trial$y1, y_off = trial$y_off),
        A = list(1, A = inla.mesh$A),
        effects = list(
            data.frame(
                b0 = rep(1, nrow(trial)),
                b1 = ifelse(trial$arm == "intervention", 1, 0),
                pvar = trial$pvar, id = trial$id, cluster = trial$cluster
            ),
            s = inla.mesh$indexs
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
    # The DIC is penalised to allow for estimation of beta2
    loss <- result.e$dic$family.dic + 2
    cat("\rDIC: ", loss, " Contamination parameter: ", beta2, "  \r")
    return(loss)
}

# Williams mean and confidence intervals
Williams <- function(x, alpha, option){
    logx_1 <- log(x + 1)
    if(length(logx_1) < 3) {
        value <- switch(option,
                M = mean(logx_1),
                L = NA,
                U = NA)
    } else {
        t <- stats::t.test(x = logx_1, conf.level = 1 - alpha)
        value <- switch(option,
                M = t$estimate,
                L = t$conf.int[1],
                U = t$conf.int[2])
    }
    returnvalue <- as.numeric(exp(value) - 1)
    return(returnvalue)
}

