
pbart2 = function(x.train, y.train, mc.cores=1, nkeeptrain=0, ...) {
  
  keeptrainfits = FALSE
  if (nkeeptrain!=0) {
    keeptrainfits = TRUE
  }
  
  t = proc.time()
  x.train.mat = dummify(x.train)
    model.matrix(~0+., data=x.train) %>% as.matrix() %>% t()
  sink("/dev/null")  
  bm = BART::mc.pbart(x.train = x.train.mat, y.train = y.train, transposed=TRUE, keeptrainfits = keeptrainfits, mc.cores=mc.cores, ...)
  /sink()  
  #bm = BART::pbart(x.train = x.train.mat, y.train = y.train, transposed=TRUE, nkeeptrain=nkeeptrain, ...)

  dur = proc.time() - t
  bm$y.train = y.train
  bm$elapsed_time = dur[3]
  bm$x_var_names = names(x.train)
  return(bm)
}


# Extracts posterior draws form pbart model
# Function handles conversion to matrix consistent with
# pbart2 function, and converts predictions from 
# probit to probability
pbart_posterior = function(pbart_fit, newdata, mc.cores=1, type="response", return_posterior_mean=FALSE) {
  # Convert data to dummified matrix format
  newdata.mat = model.matrix(~0+., 
                             data=select(newdata, one_of(pbart_fit$x_var_names))) %>% 
    as.matrix() 
  # Get predictions for newdata and convert from probits to probabilities
  sink("/dev/null")
  pos = predict(pbart_fit, newdata=newdata.mat, mc.cores=mc.cores) %>%
    t() # Transpose so that columns correspond to draws
  sink()
  
  if (type=="response") {
    pos = pnorm(pos)
  }
  
  if (return_posterior_mean) {
    pos = rowMeans(pos)
  } else {
    pos = as_tibble(pos)
    names(pos) = sprintf("draw_%s", seq_along(pos))
  }
  pos
}

# return tibble frame with factors converted to dummies
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
  
  mm = model.matrix(~0+., data = df) %>% 
    as.tibble()
  names(mm) = stringr::str_replace_all(names(mm), "[ -]+", "_")
  mm
  
}