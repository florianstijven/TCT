#' Analyze a list of trials
#'
#' The [analyze_trials()] function analyzes a list of data frames where each
#' data frame corresponds to a single clinical trial. The fitted-model object
#' for each trial is returned as a list of model-fit objects. That list can
#' further be used to extract the relevant estimates or p-values for each trial.
#'
#' @param list_data_trials A list of data frames
#' @param trial_function A function that takes the data frame as argument, fits
#'   a model, and returns the corresponding fitted-model object.
#' @param n_cores The number of cores used for parallel computing. If
#' `n_cores = 1`, no parallel computng is performed.
#' @param ... Additional parameters that are passed on to `trial_function()`.
#'
#' @return A list of fitted-model objects.
#' @export
#'
#' @examples
analyze_trials = function(list_data_trials,
                          trial_function,
                          n_cores = 1,
                          ...) {
  if (n_cores == 1) {
    fitted_models = lapply(
      X = list_data_trials,
      FUN = function(x, ...) {
        tryCatch(
          expr = trial_function(x, ...),
          error = function(e) {
            return(NA)
          }
        )
      },
      ...
    )
  }
  else {
    cl = parallel::makeCluster(n_cores)
    # parallel::clusterExport(cl = cl, varlist = env_vars)
    parallel::clusterEvalQ(cl = cl, library(tidyverse))
    parallel::clusterEvalQ(cl = cl, library(TCT))
    fitted_models = parallel::parLapply(
      cl = cl,
      X = list_data_trials,
      fun = function(x, ...) {
        gc()
        tryCatch(
          expr = trial_function(x, ...),
          error = function(e) {
            return(NA)
          }
        )
      },
      ...
    )
    parallel::stopCluster(cl)
  }
  return(fitted_models)
}


