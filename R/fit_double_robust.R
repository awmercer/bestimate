
fit_double_robust = function(y_train,
                             x_train,
                             propensity_train,
                             x_test, 
                             propensity_test,
                             bart_params,
                             seed=NULL
                             ) {
  run_bart_fit(y_train = y_train,
               x_train = bind_cols(x_train, propensity_train),
               x_test = bind_cols(x_test, propensity_test),
               bart_params = bart_params
                )
}