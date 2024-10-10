
cv_xgb <- function(df, y, weigths = NULL,  X = "X") {
  
  dtrain <- xgb.DMatrix(data = as.matrix(df), label = y, weight = weights)
  
  if (all(y %in% c(0, 1))) {
    
    params <- list(objective = "binary:logistic", eval_metric = "logloss")
  } else {
    
    params <- list(objective = "reg:squarederror", eval_metric = "rmse")
  }
  
  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    prediction = TRUE,
    verbose = FALSE
  )
  
  xgb <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = cv$best_iteration,
    verbose = FALSE
  )
  attr(xgb, "binary") <- binary
  
  xgb
}

pred_xgb <- function(xgb, df_test, intervention = NULL, X = "X") {
  
  if (!is.null(intervention)) {
    
    df_test[[X]] <- intervention
  }
  
  predict(xgb, as.matrix(df_test))
}

est_yx1wx0 <- function(sclass, n, seed, nested_mean = c("wregr", "refit")) {
  
  nested_mean <- match.arg(nested_mean, c("wregr", "refit"))
  set.seed(seed)
  
  # get the data and SFM
  c(data, SFM, gt) %<-% gen_from_scm(sclass, n)
  c(X, Z, W, Y) %<-% SFM
  
  # get the one-step debiased estimator
  
  # split into K folds
  K <- 10
  folds <- sample(x = rep(1:K, each = n / K))
  
  # for each fold, get out-of-fold fitted values
  px <- y0 <- y1 <- rep(0, n)
  y <- data[[Y]]
  x <- data[[X]]
  pseudo_out <- rep(NA, n)
  for (i in seq_len(K)) {
    
    tst <- folds == i
    dev <- folds %in% setdiff(seq_len(K), i)[1:6]
    val <- folds %in% setdiff(seq_len(K), i)[7:9]
    
    # develop models on dev
    mod_x_z <- cv_xgb(df[dev, Z], df[dev, X])
    mod_x_zw <- cv_xgb(df[dev, c(Z, W)], df[dev, X])
    mod_y_xzw <- cv_xgb(df[dev, c(X, Z, W)], df[dev, Y])
    
    # get the nested mean
    if (nested_mean == "wregr") {
      
      px_z_val <- pred_xgb(mod_x_z, df[val, Z])
      px_zw_val <- pred_xgb(mod_x_zw, df[val, c(Z, W)])
      
      weights <- ifelse(x == 1, px_z_val / (1 - px_z_val) * (1 - px_zw_val) / px_zw_val,
                        px_z_val / (1 - px_z_val))
      
      mod_nested <- cv_xgb(df[val, Z], df[val, Y], weights = weights)
      ey_nest_tst <- pred_xgb(mod_nested, df[val, Z])
    } else if (nested_mean == "refit") {
      
      y_tilde <- pred_xgb(mod_y_xzw, df[val, c(X, Z, W)], intervention = 1)
      mod_nested <- cv_xgb(df[val, c(X, Z)], y_tilde)
      ey_nest_tst <- pred_xgb(mod_nested, df[tst, c(X, Z)], intervention = 0)
    }
    
    # get the pseudo-outcome on the test fold
    ey_x1zw_tst <- pred_xgb(mod_y_xzw, df[tst, c(X, Z, W)], intervention = 1)
    px_zw_tst <- pred_xgb(mod_x_zw, df[tst, c(Z, W)])
    px_z_tst <- pred_xgb(mod_x_zw, df[tst, Z])

    pseudo_out[tst] <-
      x[tst] * (y[tst] - ey_x1zw_tst) * (1 - px_zw_tst) / px_zw_tst * 1 / (1 - px_z_tst) +
      (1-x)[tst] / (1 - px_z_tst) * (ey_x1zw_tst - ey_nest_tst) +
      ey_nest_tst
  }
  
  # compute the debiased estimator
  psi_osd <- mean(pseudo_out)
  lwr <- psi_osd - 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
  upr <- psi_osd + 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
  
  # get the ground truth
  gt <- compute_PO(sclass, c(0, 1), 0, 1, n = n)
  
  data.frame(#psi_pi = psi_pi, 
             psi_osd = psi_osd, 
             lwr = lwr, upr = upr, # confidence intervals
             gt = gt, seed = seed)
}