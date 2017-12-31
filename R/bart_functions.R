
#' Convenience wrapper for pbart from the BART package. 
#'
#' @param x.train Data frame containing the X variables for the training data
#' @param y.train Vector containing the Y variable for the training data
#' @param x.test  Xata frame containing the X variables for the test data
#' @param mc.cores Number of cores to use
#' @param ndpost Number of posterior samples to draw
#' @param keeptrainfits Keep samples for training data?
#' @param verbose Show pbart progress messages?
#' @param ... Additional BART parameters passed to \code{pbart}
#'
#' @return Returns a pbart object
#' @export
#'
pbart2 = function(x.train, 
                  y.train, 
                  x.test = NULL, 
                  mc.cores=1, 
                  ndpost = 1000, 
                  keeptrainfits = FALSE, 
                  verbose = FALSE,
                  ...) {
  
  if (keeptrainfits) {
    nkeeptrain = ndpost
  } else {
    nkeeptrain = 0
  }

  t = proc.time()
  x.train.mat = dummify(x.train) %>% as.matrix() %>% t()
  
  if (!is.null(x.test)) {
    x.test.mat = dummify(x.test) %>% as.matrix() %>%t()
  } else {
    x.test.mat = matrix(0.0,0,0)
  }
  
  if (!verbose) sink("/dev/null")  
  if (mc.cores > 1) {
    bm = BART::mc.pbart(x.train = x.train.mat, 
                        y.train = y.train, 
                        x.test = x.test.mat,
                        ndpost = ndpost,
                        transposed=TRUE, 
                        keeptrainfits = keeptrainfits, 
                        mc.cores=mc.cores,
                        ...)
  } else {
    bm = BART::pbart(x.train = x.train.mat, 
                     y.train = y.train, 
                     x.test = x.test.mat,
                     ndpost = ndpost,
                     transposed=TRUE, 
                     nkeeptrain=nkeeptrain, 
                     ...)
  }
  if (!verbose) sink()  
  dur = proc.time() - t
  bm$elapsed_time = dur[[3]]
  bm$x_var_names = names(x.train)
  return(bm)
}



### ----------------------------------------------------------




#' Get posterior samples from pbart object for newdata
#'
#' @param pbart_fit Object of class "\code{pbart}"
#' @param newdata Data frame containing X variables for use in generating
#' new predictions
#' @param mc.cores Number of cores to use 
#' @param return_posterior_mean If TRUE returns the posterior mean for each
#' observation in \code{newdata}. Otherwise returns the full set of posterior
#' draws.  Defaults to \code{FALSE}.
#'
#' @return Returns a tibble containing predicted values. Either a single
#' posterior mean or the full set of posterior samples.
#' @export
#'
pbart_posterior = function(pbart_fit, newdata, mc.cores=1, return_posterior_mean=FALSE) {
  # Convert data to dummified matrix format
  newdata.mat = select(newdata, one_of(pbart_fit$x_var_names)) %>%
    dummify() %>%
    as.matrix() 
  # Get predictions for newdata and convert from probits to probabilities
  sink("/dev/null")
  pos = predict(pbart_fit, newdata=newdata.mat, mc.cores=mc.cores) %>%
    t() # Transpose so that columns correspond to draws
  sink()
  
  pos = pnorm(pos)

  if (return_posterior_mean) {
    pos = rowMeans(pos)
  } else {
    pos = as_tibble(pos)
    names(pos) = sprintf("draw_%s", seq_along(pos))
  }
  pos
}



### ----------------------------------------------------------


#' Convert factors in data frame to binary variables
#'
#' @param df Data frame to convert
#' @param drop_unused Whether or not to drop levels that do not have any
#' corresponding observations.
#' @param sep Separator between variable and value names 
#'
#' @return Returns a tibble where the factors have been converted to binary
#' variables
#' @export
#'
dummify = function(df, drop_unused=TRUE, sep="_") {
  df = mutate_if(df, is.character, as.factor)
  if (drop_unused) {
    df = mutate_if(df, is.factor, forcats::fct_drop)
  }
  
  # Convert any factors with a single level to numeric with value of 1
  df = rename_if(df, 
                  ~is.factor(.) & length(levels(.))==1,
                  ~sprintf("%s_%s", ., map(df[.], levels))) %>%
    mutate_if(~is.factor(.) & length(levels(.))==1, ~1)

  df = rename_if(df, is.factor, function(n) sprintf("%s%s", n, sep))
  
  mm = model.matrix(~0+., data = df, 
                    contrasts.arg = df %>% 
                      select_if(is.factor) %>%
                      select_if(~length(levels(.))>2)%>% 
                      map(contrasts, contrasts=FALSE)
                    ) %>% 
    as.tibble()
  names(mm) = stringr::str_replace_all(names(mm), "[ -]+", "_")
  mm

}

