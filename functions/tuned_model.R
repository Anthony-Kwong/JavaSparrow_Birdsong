#' tuned_model constructor
#' 
#' An S3 object to hold the results of model_tune
#'
#' @param model : The type of statistical algorithm used (e.g. logistic_regression)
#' @param best_params : The best set of tuning parameters based on cross validation.
#' @param tune_tibble : A tibble containing the cv accuracy for each set of tuning parameters used. 
#' @param comp_time : Total time taken to tune the model
#' @param fitted_model : A parsnip model. The model with the best set of tuning params, fitted on the
#' entire training dataset. 
#'
#' @return model_tune S3 object
#' @export
#'
#' @examples
tuned_model <- function(model, best_params, tune_tibble, comp_time, fitted_model){
  value <- list(model = model, 
                best_params = best_params,
                tune_tibble = tune_tibble, 
                comp_time = comp_time,
                fitted_model = fitted_model)
  attr(value, "class") <- "tuned_model"
  return (value)
}