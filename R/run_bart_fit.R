
#' Wrapper to run a bart model and return only posterior distributions
#'
#' @param y_train 
#' @param x_train 
#' @param x_test 
#' @param bart_params 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
run_bart_fit = function(y_train,
         x_train, 
         x_test, 
         bart_params,
         seed=NULL) {

  if (!is.null(seed)) set.seed(seed)
  
  fit = do.call(pbart2, c(list(x.train = x_train, 
                         y.train = y_train, 
                         x.test = x_test), bart_params))
  
  
  train_posterior = fit$yhat.train %>% 
    t() %>% 
    pnorm() %>%
    as_tibble() %>%
    set_names(sprintf("draw_%s", seq_along(.)))
  
  test_posterior = fit$yhat.test %>% 
    t() %>% 
    pnorm() %>%
    as_tibble() %>%
    set_names(sprintf("draw_%s", seq_along(.)))
  
  fit$yhat.train.mean = NULL
  fit$yhat.train = NULL
  fit$yhat.test.mean = NULL
  fit$yhat.test = NULL
  
  list(train_posterior = train_posterior,
       test_posterior = test_posterior,
       fit = fit)
}