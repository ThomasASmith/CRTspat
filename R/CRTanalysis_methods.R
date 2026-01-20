#' Summary of the results of a statistical analysis of a CRT
#'
#' \code{summary.CRTanalysis} generates a summary of a \code{CRTanalysis} including the main results
#' @param object an object of class \code{"CRTanalysis"}
#' @param ... other arguments used by summary
#' @method summary CRTanalysis
#' @export
#' @return No return value, writes text to the console.
#' @examples
#' {example <- readdata('exampleCRT.txt')
#' exampleT <- CRTanalysis(example, method = "T")
#' summary(exampleT)
#' }
summary.CRTanalysis <- function(object, ...) {
  defaultdigits <- getOption("digits")
  on.exit(options(digits = defaultdigits))
  options(digits = 3)
  scale_par <- object$options$scale_par
  alpha <- object$options$alpha
  cat("\n=====================CLUSTER RANDOMISED TRIAL ANALYSIS =================\n")
  cat(
    "Analysis method: ", object$options$method, "\nLink function: ", object$options$link, "\n")
  if(!identical(object$options$distance_type,"No fixed effects of distance ")){
    cat(paste0("Measure of distance or surround: ", getDistanceText(distance = object$options$distance,
                                                                    scale_par = scale_par),"\n"))
    if (!is.null(object$options$linearity)){
      cat(object$options$linearity, "\n")
    }
  }
  if (!is.null(object$options$ftext))
    cat("Model formula: ", object$options$ftext, "\n")
  cat(switch(object$options$cfunc,
             Z = "No comparison of arms \n",
             X = "No modelling of spillover \n",
             S = "Piecewise linear function for spillover\n",
             P = "Error function model for spillover\n",
             L = "Sigmoid (logistic) function for spillover\n",
             R = "Rescaled linear function for spillover\n"))
  dummy <- printrow(var = 'controlY', pt_src = 'pt_ests', int_src = 'int_ests',
                    object = object, text = "Estimates:       Control: ", alpha = alpha)
  dummy <- printrow(var = 'interventionY', pt_src = 'pt_ests', int_src = 'int_ests',
                    object = object, text = "            Intervention: ", alpha = alpha)
  effect_distance <- ifelse(object$options$link == 'identity', "             Effect size: ",
                                                               "                Efficacy: ")
  dummy <- printrow(var = 'effect_size', pt_src = 'pt_ests', int_src = 'int_ests',
                    object = object, text = effect_distance, alpha = alpha)

  dummy <- printrow(var = 'personal_protection', pt_src = 'pt_ests', int_src = 'int_ests',
                    object = object, text = "Personal protection     : ", alpha = alpha)

  if (!is.null(object$pt_ests$personal_protection)){
    if (!is.na(object$pt_ests$personal_protection)){
      if (object$pt_ests$personal_protection < 0 | object$pt_ests$personal_protection >  1){
          cat(
            "** Warning: different signs for main effect and personal protection effect:
                face validity check fails **\n")
      }
    }
  }
  dummy <- printrow(var = 'scale_par', pt_src = 'pt_ests', int_src = 'int_ests',
                    object = object, text = "         Scale parameter: ", alpha = alpha)
  if (!is.na(object$pt_ests$spillover_interval)){
    dummy <- printrow(var = 'spillover_interval', pt_src = 'spillover', int_src = 'int_ests',
                    object = object, text = "  Spillover interval(km): ", alpha = alpha)
    dummy <- printrow(var = 'midpoint', pt_src = 'spillover', int_src = 'int_ests',
                    object = object, text = "  Spillover midpoint(km): ", alpha = alpha)
    dummy <- printrow(var = 'contaminate_pop_pr', pt_src = 'spillover', int_src = 'int_ests',
                    object = object, text = "Pr in spillover interval: ", alpha = alpha)
    dummy <- printrow(var = 'total_effect', pt_src = 'spillover', int_src = 'int_ests',
                    object = object, text = "            Total effect: ", alpha = alpha)
    dummy <- printrow(var = 'ipsilateral_spillover', pt_src = 'spillover', int_src = 'int_ests',
                    object = object, text = "   Ipsilateral Spillover: ", alpha = alpha)
    dummy <- printrow(var = 'contralateral_spillover', pt_src = 'spillover', int_src = 'int_ests',
                    object = object, text = " Contralateral Spillover: ", alpha = alpha)
  }
  if (!is.null(object$description$cv_percent)) object$int_ests$cv_percent <-
                    c(object$description$cv_lower, object$description$cv_upper)
  dummy <- printrow(var = 'cv_percent', pt_src = 'description', int_src = 'int_ests',
                    object = object, text = "Coefficient of variation: ", alpha = alpha)
  dummy <- printrow(var = 'ICC', pt_src = 'pt_ests', int_src = 'int_ests',
                    object = object, text = "Intracluster correlation: ", alpha = alpha)
  options(digits = defaultdigits)
  # goodness of fit
  if (!is.null(object$pt_ests$deviance)) cat("deviance                : ", object$pt_ests$deviance, "\n")
  if (!is.null(object$pt_ests$DIC)) cat("DIC                     : ", object$pt_ests$DIC)
  if (!is.null(object$pt_ests$elpd)) cat("elpd                    : ", object$pt_ests$elpd)
  if (!is.null(object$pt_ests$AIC)) cat("AIC                     : ", object$pt_ests$AIC)
  if (object$options$penalty > 0) {
    cat(" including penalty for the spillover scale parameter\n")
  } else {
    cat(" \n")
  }
  # TODO: add the degrees of freedom to the output
  if (!is.null(object$pt_ests$p.value)){
    cat("P-value (2-sided): ", object$pt_ests$p.value, "\n")
  }
}

#' Extract model fitted values
#'
#' \code{fitted.CRTanalysis} method for extracting model fitted values
#' @param object CRTanalysis object
#' @param ... other arguments
#' @export
#' @return the fitted values returned by the statistical model run within the \code{CRTanalysis} function
#' @examples
#' {example <- readdata('exampleCRT.txt')
#' exampleGEE <- CRTanalysis(example, method = "GEE")
#' fitted_values <- fitted(exampleGEE)
#' }
fitted.CRTanalysis <- function(object, ...){
  value = fitted(object = object$model_object, ...)
  return(value)
}

#' Extract model coefficients
#'
#' \code{coef.CRTanalysis} method for extracting model fitted values
#' @param object CRTanalysis object
#' @param ... other arguments
#' @export
#' @return the model coefficients returned by the statistical model run within the \code{CRTanalysis} function
#' @examples
#' {example <- readdata('exampleCRT.txt')
#' exampleGEE <- CRTanalysis(example, method = "GEE")
#' coef(exampleGEE)
#' }
coef.CRTanalysis <- function(object, ...){
  value = coef(object = object$model_object, ...)
  return(value)
}

#' Extract model residuals
#'
#' \code{residuals.CRTanalysis} method for extracting model residuals
#' @param object CRTanalysis object
#' @param ... other arguments
#' @export
#' @return the residuals from the statistical model run within the \code{CRTanalysis} function
#' @examples
#' {example <- readdata('exampleCRT.txt')
#' exampleGEE <- CRTanalysis(example, method = "GEE")
#' residuals <- residuals(exampleGEE)
#' }
residuals.CRTanalysis <- function(object, ...){
  if ("gee" %in% class(object$model_object)) {
    value <- object$model_object$residuals
  } else {
    value <- residuals(object = object$model_object, ...)
  }
  return(value)
}

#' Model predictions
#'
#' \code{predict.CRTanalysis} method for extracting model predictions
#' @param object CRTanalysis object
#' @param ... other arguments
#' @export
#' @return the model predictions returned by the statistical model run within the \code{CRTanalysis} function
#' @examples
#' {example <- readdata('exampleCRT.txt')
#' exampleGEE <- CRTanalysis(example, method = "GEE")
#' predictions <- predict(exampleGEE)
#' }#'
predict.CRTanalysis <- function(object, ...){
  value = predict(object = object$model_object, ...)
  return(value)
}

printrow <- function(var, pt_src, int_src, object, text, alpha){
  CLtext <- paste0(" (", 100 * (1 - alpha), "% CL: ")
  pstring <- NULL
  if (!is.null(object[[pt_src]][[var]])){
    if (!is.na(object[[pt_src]][[var]])){
      pstring <- paste0(text, round(object[[pt_src]][[var]], digits = 3))
      if (!is.null(object[[int_src]][[var]])){
        if (!is.na(object[[int_src]][[var]][1])){
          pstring <- paste0(pstring, CLtext, round(unlist(object[[int_src]][[var]][1]), digits = 3)," ",
                            round(unlist(object[[int_src]][[var]][2]), digits = 3),")")
        }
      }
      cat(paste0(pstring, "\n"))
    }
  }
}


