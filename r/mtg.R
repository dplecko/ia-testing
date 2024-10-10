
ia_vals_mtg <- function(data, X, Z, W, Y, x0 = 0, x1 = 1, 
                        pw_model = c("logistic", "xgb"), 
                        thr = 0.5, quant = FALSE, nboot = 1L, nboot2 = 100L, 
                        norm_method = "adapt") {
  
  if (nboot > 1L) {
    
    mtgs <- lapply(
      seq_len(nboot),
      function(i) {
        
        set.seed(i)
        boot_idx <- if (i == 1) TRUE else sample(nrow(data), replace = TRUE)
        boot_data <- data[boot_idx, ]
        ret <- ia_vals_mtg(boot_data, X, Z, W, Y, x0, x1, pw_model, thr, quant, 
                           nboot = 1L, norm_method = norm_method)
        ret$outer_boot <- i
        ret
      }
    )
    return(do.call(rbind, mtgs))
  }
  
  pw_model <- match.arg(pw_model, c("logistic", "xgb"))
  
  trim_to_01 <- function(x) pmax(0, pmin(x, 1))
  thresh <- function(x, threshold) as.integer(x >= threshold)
  wgh_sum <- function(x, wgh) sum(x * wgh) / sum(wgh)
  
  compute_po <- function(fx = 0, wx = 0, zx = 0, ingr, thresh = FALSE, 
                         norm_method = "adapt") {
    
    c(probx1, probx0, px, px_z, px_zw) %<-% ingr
    if (is.na(zx)) zx <- -1
    
    if (wx == 0) wgh <- (1 - px_zw) / (1 - px_z) else wgh <- px_zw / px_z
    if (zx == 0) {
      wgh <- wgh * (1 - px_zw) / (1 - px)
    } else if (zx == 1) wgh <- wgh * px_z / px
    
    if (fx == 0) po_samp <- probx0 else po_samp <- probx1
    if (norm_method == "adapt") norm_const <- sum(wgh) else 
      norm_const <- length(po_samp)
    if (thresh) po_samp <- as.integer(po_samp > thr)
    list(po_samp = po_samp, wgh = wgh, norm_const = norm_const)
  }
  
  eval_po <- function(po) sum(po$po_samp * po$wgh) / po$norm_const
  diff_po <- function(po1, po2) eval_po(po1) - eval_po(po2)
  samp_imp <- function(po1, po2) {
    
    if (po1$norm_const != po2$norm_const) return(NA)
    po1$po_samp - po2$po_samp
  }
  
  if (is.logical(data[[Y]])) data[[Y]] <- as.integer(data[[Y]])
  binary <- if (all(data[[Y]] %in% c(0, 1))) TRUE else FALSE
  
  if (binary) { # use the logistic loss for binary classification
    
    params <- list(objective = "binary:logistic", eval_metric = "logloss",
                   nthread = 1)
  } else { # use the squared loss for regression
    
    params <- list(objective = "reg:squarederror", eval_metric = "rmse",
                   nthread = 1)
  }

  # build a predictor (xgboost) for the outcome Y
  dtrain <- xgb.DMatrix(data = as.matrix(data[, c(X, Z, W)]), 
                        label = data[[Y]])

  # Cross-validation to find optimal nrounds for each fold
  cv_fold <- xgb.cv(
    params = params, data = dtrain, nfold = 5, nrounds = 100,
    early_stopping_rounds = 10, verbose = 0, prediction = TRUE 
  )
  prob <- cv_fold$pred
  optimal_nrounds <- cv_fold$best_iteration
  
  # K-fold cross-validation setup
  K <- 5
  folds <- sample(rep(1:K, length.out = nrow(data)))
  probx0 <- probx1 <- rep(NA, nrow(data))
  
  # Loop over K folds
  for (k in 1:K) {
    
    # Split into training and validation sets
    train_idx <- which(folds != k)
    valid_idx <- which(folds == k)
    
    dtrain <- xgb.DMatrix(data = as.matrix(data[train_idx, c(X, Z, W)]), 
                          label = data[train_idx, ][[Y]])
    
    # Train the model on the current fold
    mod <- xgb.train(params = params, data = dtrain, nrounds = optimal_nrounds)
    
    # Prepare counterfactual datasets for validation fold
    data_do_x0 <- data_do_x1 <- data[valid_idx, ]
    data_do_x0[[X]] <- 0
    data_do_x1[[X]] <- 1
    
    # Get counterfactual probabilities
    probx0[valid_idx] <- predict(mod, xgb.DMatrix(data = as.matrix(data_do_x0[, c(X, Z, W)])))
    probx1[valid_idx] <- predict(mod, xgb.DMatrix(data = as.matrix(data_do_x1[, c(X, Z, W)])))
  }
  
  if (binary) {
    
    probx0 <- trim_to_01(probx0)
    probx1 <- trim_to_01(probx1)
  }
  
  if (quant) {
    
    thr <- quantile(prob, probs = thr)
  }
  
  px <- mean(data[[X]])
  if (pw_model == "logistic") {
    
    # get the propensity scores - regress X on Z, Z+W
    mod_xz <- glm(paste(X, "~", paste(Z, collapse = "+")), family = "binomial",
                  data = data)
    px_z <- mod_xz$fitted.values
    mod_xzw <- glm(paste(X, "~", paste(c(Z, W), collapse = "+")), 
                   family = "binomial", data = data)
    px_zw <- mod_xzw$fitted.values
  } else if (pw_model == "xgb") {
    
    params <- list(objective = "binary:logistic", eval_metric = "logloss",
                   nthread = 1)
    
    #' * out-of-bag xgboost probabilities P(X | Z, W) *
    dtrain <- xgb.DMatrix(data = as.matrix(data[, c(Z, W)]), label = data[[X]])
    cv_zw <- xgb.cv(
      params = params, data = dtrain, nfold = 5, nrounds = 100,
      early_stopping_rounds = 10, verbose = 0, prediction = TRUE
    )
    px_zw <- cv_zw$pred
    
    #' * out-of-bag xgboost probabilities P(X | Z) *
    dtrain <- xgb.DMatrix(data = as.matrix(data[, c(Z)]), label = data[[X]])
    cv_z <- xgb.cv(
      params = params, data = dtrain, nfold = 5, nrounds = 100,
      early_stopping_rounds = 10, verbose = 0, prediction = TRUE
    )
    px_z <- cv_z$pred
  }
  
  
  res <- c()
  for (i in seq_len(nboot2)) { 
    
    b_idx <- if (i == 1) TRUE else sample(nrow(data), replace = TRUE)
    ingr <- list(probx1[b_idx], probx0[b_idx], px, px_z[b_idx], px_zw[b_idx])
    
    fx1wx0_x0 <- compute_po(fx = 1, wx = 0, zx = 0, ingr, FALSE, norm_method)
    fx0wx0_x0 <- compute_po(fx = 0, wx = 0, zx = 0, ingr, FALSE, norm_method)
    fx1wx1_x0 <- compute_po(fx = 1, wx = 1, zx = 0, ingr, FALSE, norm_method)
    fx0wx1_x0 <- compute_po(fx = 0, wx = 1, zx = 0, ingr, FALSE, norm_method)
    
    fx1wx0_x1 <- compute_po(fx = 1, wx = 0, zx = 1, ingr, FALSE, norm_method)
    fx0wx0_x1 <- compute_po(fx = 0, wx = 0, zx = 1, ingr, FALSE, norm_method)
    fx1wx1_x1 <- compute_po(fx = 1, wx = 1, zx = 1, ingr, FALSE, norm_method)
    fx0wx1_x1 <- compute_po(fx = 0, wx = 1, zx = 1, ingr, FALSE, norm_method)
    
    # TE-SE
    tese <- diff_po(fx1wx1_x0, fx0wx0_x0) - diff_po(fx1wx1_x1, fx0wx0_x1)
    
    # DE-IE
    deie_x0 <- diff_po(fx1wx0_x0, fx0wx0_x0) - diff_po(fx1wx1_x0, fx0wx1_x0)
    deie_x1 <- diff_po(fx1wx0_x1, fx0wx0_x1) - diff_po(fx1wx1_x1, fx0wx1_x1)
    
    # DE-SE
    dese <- diff_po(fx1wx0_x0, fx0wx0_x0) - diff_po(fx1wx0_x1, fx0wx0_x1)
    
    # IE-SE
    iese <- diff_po(fx1wx0_x1, fx1wx1_x1) - diff_po(fx1wx0_x0, fx1wx1_x0)
    
    # DE-IE-SE
    deiese <- deie_x0 - deie_x1
    
    res <- rbind(
      res,
      data.frame(value = tese, inner_boot = i, ia = "TE x SE"),
      data.frame(value = deie_x0, inner_boot = i, ia = "DE x IE"),
      data.frame(value = dese, inner_boot = i, ia = "DE x SE"),
      data.frame(value = iese, inner_boot = i, ia = "IE x SE"),
      data.frame(value = deiese, inner_boot = i, ia = "DE x IE x SE")
    )
  }

  res
}
