
cv_xgb <- function(df, y, weights = NULL,  X = "X") {
  
  dtrain <- xgb.DMatrix(data = as.matrix(df), label = y, weight = weights)
  
  binary <- all(y %in% c(0, 1))
  if (binary) {
    
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

ia_tests <- function(data, X, Z, W, Y, nested_mean = c("refit", "wregr"),
                     log_scale = FALSE) {
  
  n <- nrow(data)
  nested_mean <- match.arg(nested_mean, c("refit", "wregr"))
  
  # split into K folds
  K <- 10
  folds <- sample(x = rep(1:K, each = n / K))
  
  y <- data[[Y]]
  x <- data[[X]]
  pso <- list(list(list(list(), list()), list(list(), list())),
              list(list(list(), list()), list(list(), list())))
  
  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
    
    pso[[xz+1]][[xw+1]][[xy+1]] <- rep(NA, n)
  }
  
  # get the simplest pseudo-outcomes (no cross-fitting needed)
  for (xz in c(0, 1)) {
    
    pso[[xz+1]][[xz+1]][[xz+1]] <- (x == xz) / mean(x == xz) * y
  }
  
  for (i in seq_len(K)) {
    
    # split into dev, val, tst
    tst <- folds == i
    dev <- folds %in% setdiff(seq_len(K), i)[1:6]
    val <- folds %in% setdiff(seq_len(K), i)[7:9]
    
    # develop models on dev
    mod_x_z <- cv_xgb(data[dev, Z], data[dev, X])
    mod_x_zw <- cv_xgb(data[dev, c(Z, W)], data[dev, X])
    mod_y_xz <- cv_xgb(data[dev, c(X, Z)], data[dev, Y])
    mod_y_xzw <- cv_xgb(data[dev, c(X, Z, W)], data[dev, Y])
    
    # get the val set predictions
    px_zw_val <- pred_xgb(mod_x_zw, data[val, c(Z, W)])
    px_zw_val <- list(1 - px_zw_val, px_zw_val)
    px_z_val <- pred_xgb(mod_x_z, data[val, Z])
    px_z_val <- list(1 - px_z_val, px_z_val)
    y_xzw_val <- list(
      pred_xgb(mod_y_xzw, data[val, c(X, Z, W)], intervention = 0, X = X),
      pred_xgb(mod_y_xzw, data[val, c(X, Z, W)], intervention = 1, X = X)
    )
    
    # get the test set values
    px_zw_tst <- pred_xgb(mod_x_zw, data[tst, c(Z, W)])
    px_zw_tst <- list(1 - px_zw_tst, px_zw_tst)
    px_z_tst <- pred_xgb(mod_x_z, data[tst, Z])
    px_z_tst <- list(1 - px_z_tst, px_z_tst)
    y_xzw_tst <- list(
      pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)], intervention = 0, X = X),
      pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)], intervention = 1, X = X)
    )
    y_xz_tst <- list(
      pred_xgb(mod_y_xz, data[tst, c(X, Z)], intervention = 0, X = X),
      pred_xgb(mod_y_xz, data[tst, c(X, Z)], intervention = 1, X = X)
    )
    
    # get the total-effect pseudo-outcomes (ETT-like)
    for (xz in c(0, 1)) {
      
      xw <- xy <- 1 - xz
      pso[[xz+1]][[xw+1]][[xy+1]][tst] <- 
        ((x[tst] == xw) / mean(x == xz) * px_z_tst[[xz+1]] / px_z_tst[[xw+1]] * (y[tst] - y_xz_tst[[xw+1]]) + 
           (x[tst] == xz) / mean(x == xz) * y_xz_tst[[xw+1]])
    }
    
    # get the nested mean
    ey_nest_tst <- list(NULL, NULL)
    for (xw in c(0, 1)) {
      
      xy <- 1 - xw
      if (nested_mean == "wregr") {
        
        weights <- ifelse(
          x[val] == xy, 
          px_z_val[[xy+1]] / px_z_val[[xw+1]] * 
            px_zw_val[[xw+1]] / px_zw_val[[xy+1]],
          px_z_val[[xy+1]] / px_z_val[[xw+1]]
        )
        
        mod_nested <- cv_xgb(data[val, Z], data[val, Y], weights = weights)
        ey_nest_tst[[xy+1]] <- pred_xgb(mod_nested, data[val, Z])
      } else if (nested_mean == "refit") {
        
        y_tilde <- pred_xgb(mod_y_xzw, data[val, c(X, Z, W)], 
                            intervention = xy, X = X)
        mod_nested <- cv_xgb(data[val, c(X, Z)], y_tilde)
        ey_nest_tst[[xy+1]] <- pred_xgb(mod_nested, data[tst, c(X, Z)], 
                                        intervention = xw, X = X)
      }
    }
    
    for (xz in c(0, 1)) {
      
      for (xw in c(0, 1)) {
        
        xy <- 1 - xw
        pso[[xz+1]][[xw+1]][[xy+1]][tst] <- 
          # Term T1
          (x[tst] == xy) * (y[tst] - y_xzw_tst[[xy+1]]) * 
          px_zw_tst[[xw+1]] / px_zw_tst[[xy+1]] * 
          px_z_tst[[xz+1]] / px_z_tst[[xw+1]] * 1 / mean(x == xz) +
          # Term T2
          (x[tst] == xw) / mean(x == xz) * px_z_tst[[xz+1]] / px_z_tst[[xw+1]] *  
          (y_xzw_tst[[xy+1]] - ey_nest_tst[[xy+1]]) +
          # Term T3
          (x[tst] == xz) / mean(x == xz) * ey_nest_tst[[xy+1]] 
      }
    }
  }
  
  ias <- list(
    tese = list(
      sgn = c(1, -1, -1, 1),
      spc = list(
        c(0, 1, 1), c(0, 0, 0), c(1, 1, 1), c(1, 0, 0)
      ),
      ia = "TE x SE"
    ),
    deie0 = list(
      sgn = c(1, -1, -1, 1),
      spc = list(
        c(0, 0, 1), c(0, 0, 0), c(0, 1, 1), c(0, 1, 0)
      ),
      ia = "DE x IE"
    ),
    # deie1 = list(
    #   sgn = c(1, -1, -1, 1),
    #   spc = list(
    #     c(1, 0, 1), c(1, 0, 0), c(1, 1, 1), c(1, 1, 0)
    #   ),
    #   ia = "DE x IE"
    # ),
    dese = list(
      sgn = c(1, -1, -1, 1),
      spc = list(
        c(0, 0, 1), c(0, 0, 0), c(1, 0, 1), c(1, 0, 0)
      ),
      ia = "DE x SE"
    ),
    iese = list(
      sgn = c(1, -1, -1, 1),
      spc = list(
        c(0, 1, 0), c(0, 0, 0), c(1, 1, 0), c(1, 0, 0)
      ),
      ia = "IE x SE"
    ),
    deiese = list(
      sgn = c(c(1, -1, -1, 1), -c(1, -1, -1, 1)),
      spc = list(
        c(0, 0, 1), c(0, 0, 0), c(0, 1, 1), c(0, 1, 0),
        c(1, 0, 1), c(1, 0, 0), c(1, 1, 1), c(1, 1, 0)
      ),
      ia = "DE x IE x SE"
    ),
    xde = list(
      sgn = c(1, -1),
      spc = list(c(0, 0, 1), c(0, 0, 0)),
      ia = "x-DE"
    ),
    xie = list(
      sgn = c(1, -1),
      spc = list(c(0, 0, 1), c(0, 1, 1)),
      ia = "x-IE"
    ),
    xse = list(
      sgn = c(1, -1),
      spc = list(c(0, 1, 1), c(1, 1, 1)),
      ia = "x-SE"
    ),
    ett = list(
      sgn = c(1, -1),
      spc = list(c(0, 1, 1), c(0, 0, 0)),
      ia = "x-TE"
    )
  )
  
  res <- c()
  for (i in seq_along(ias)) {
    
    if (log_scale) {
      
      scale <- "log"
      psi_osd <- var_psi_osd <- 0
      for (j in seq_along(ias[[i]]$sgn)) {
        
        xz <- ias[[i]]$spc[[j]][1]
        xw <- ias[[i]]$spc[[j]][2]
        xy <- ias[[i]]$spc[[j]][3]
        
        # update the \hat\psi
        psi_i <- mean(pso[[xz+1]][[xw+1]][[xy+1]])
        # update log-scale \psi
        psi_osd <- psi_osd + ias[[i]]$sgn[j] * log(psi_i)
        # variance update (log-scale)
        var_psi_osd <- var_psi_osd + var(pso[[xz+1]][[xw+1]][[xy+1]]) / psi_i^2 
      }
      
      var_psi_osd <- var_psi_osd / length(pso[[xz+1]][[xw+1]][[xy+1]])
      dev <- sqrt(var_psi_osd)
      lwr <- psi_osd - 1.96 * dev
      upr <- psi_osd + 1.96 * dev
    } else {
      
      scale <- "difference"
      pseudo_out <- 0
      for (j in seq_along(ias[[i]]$sgn)) {
        
        xz <- ias[[i]]$spc[[j]][1]
        xw <- ias[[i]]$spc[[j]][2]
        xy <- ias[[i]]$spc[[j]][3]
        pseudo_out <- pseudo_out + ias[[i]]$sgn[j] * pso[[xz+1]][[xw+1]][[xy+1]]
      }
      psi_osd <- mean(pseudo_out)
      dev <- sqrt(var(pseudo_out) / length(pseudo_out))
      lwr <- psi_osd - 1.96 * dev 
      upr <- psi_osd + 1.96 * dev
    }

    res <- rbind(
      res,
      data.frame(psi_osd = psi_osd, lwr = lwr, upr = upr, dev = dev, 
                 ia = ias[[i]]$ia, scale = scale)
    )
  }
  
  # res <- c()
  # for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
  #   
  #   gt <- compute_PO(sclass, xz, xw, xy, n = n)
  #   psi_osd <- mean(pso[[xz+1]][[xw+1]][[xy+1]])
  #   lwr <- psi_osd - 1.96 * sqrt(var(pso[[xz+1]][[xw+1]][[xy+1]]) / n)
  #   upr <- psi_osd + 1.96 * sqrt(var(pso[[xz+1]][[xw+1]][[xy+1]]) / n)
  #   res <- rbind(
  #     res,
  #     data.frame(xz = xz, xw = xw, xy = xy, psi_osd = psi_osd, lwr = lwr, 
  #                upr = upr, gt = gt, seed = seed)
  #   )
  # }
  
  as.data.table(res)
}

# cv_xgb <- function(df, y, df_test, X = "X") {
#   
#   dtrain <- xgb.DMatrix(data = as.matrix(df), label = y)
#   
#   if (all(y %in% c(0, 1))) {
#     
#     params <- list(objective = "binary:logistic", eval_metric = "logloss")
#   } else {
#     
#     params <- list(objective = "reg:squarederror", eval_metric = "rmse")
#   }
#   
#   cv <- xgb.cv(
#     params = params,
#     data = dtrain,
#     nrounds = 1000,
#     nfold = 5,
#     early_stopping_rounds = 10,
#     prediction = TRUE,
#     verbose = FALSE
#   )
#   
#   if (!is.null(df_test)) {
#     
#     xgb <- xgb.train(
#       params = params,
#       data = dtrain,
#       nrounds = cv$best_iteration,
#       verbose = FALSE
#     )
#     
#     if (X %in% names(df)) {
#       
#       df_x0 <- df_x1 <- df_test
#       df_x0[[X]] <- 0
#       df_x1[[X]] <- 1
#       return(
#         list(
#           po = predict(xgb, as.matrix(df_test)),
#           po0 = predict(xgb, as.matrix(df_x0)),
#           po1 = predict(xgb, as.matrix(df_x1))
#         )
#       )
#     } else {
#       
#       return(list(po = predict(xgb, as.matrix(df_test))))
#     }
#   }
#   
#   return(cv$pred)
# }
# 
# est_yx1_x0 <- function(sclass, n, seed) {
#   
#   set.seed(seed)
#   
#   # get the data and SFM
#   c(data, SFM, gt) %<-% gen_from_scm(sclass, n)
#   c(X, Z, W, Y) %<-% SFM
#   
#   # get the \hat P estimator (plug-in)
#   pos <- cv_xgb(data[, c(X, Z)], data[[Y]], data[, c(X, Z)])
#   psi_pi <- (
#     mean(pos[["po1"]][data[[X]] == 0]) -  mean(pos[["po0"]][data[[X]] == 0])) - (
#       mean(pos[["po1"]][data[[X]] == 1]) - mean(pos[["po0"]][data[[X]] == 1])
#     )
#   
#   # get the one-step debiased estimator
#   
#   # split into K folds
#   K <- 5
#   folds <- sample(x = rep(1:K, each = n / K))
#   
#   # for each fold, get out-of-fold fitted values
#   px <- y0 <- y1 <- rep(0, n)
#   y <- data[[Y]]
#   x <- data[[X]]
#   for (i in seq_len(K)) {
#     
#     px[folds == i] <- cv_xgb(data[folds != i, Z], data[folds != i, X],
#                              data[folds == i, Z])[["po"]]
#     
#     y_pos <- cv_xgb(data[folds != i, c(X, Z)], data[folds != i, Y],
#                     data[folds == i, c(X, Z)])
#     y0[folds == i] <- y_pos[["po0"]]
#     y1[folds == i] <- y_pos[["po1"]]
#   }
#   
#   # compute the debiased estimator
#   pseudo_out <- ( 
#     # E[Y_{x_1} | x_0]
#     (x / (1 - mean(x)) * (1 - px) / px * (y - y1) + (1 - x) / (1 - mean(x)) * y1) -
#       (1 - x) / (1 - mean(x)) * y # E[Y_{x_0} | x_0]
#   ) - 
#     ( 
#       x / mean(x) * y # E[Y_{x_1} | x_1]
#       - # E[Y_{x_0} | x_1]
#         ((1 - x) / mean(x) * px / (1 - px) * (y - y0) + x / mean(x) * y0)
#     ) 
#   
#   psi_osd <- mean(pseudo_out)
#   lwr <- psi_osd - 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
#   upr <- psi_osd + 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
#   
#   # get the ground truth
#   gt <- (compute_PO(sclass, 0, 1, 1, n = n) -
#            compute_PO(sclass, 0, 0, 0, n = n)) - (
#              compute_PO(sclass, 1, 1, 1, n = n) -
#                compute_PO(sclass, 1, 0, 0, n = n)
#            )
#   
#   data.frame(psi_pi = psi_pi, psi_osd = psi_osd, 
#              lwr = lwr, upr = upr, # confidence intervals
#              gt = gt, seed = seed)
# }

# est_yx1wx0 <- function(sclass, n, seed, nested_mean = c("wregr", "refit")) {
#   
#   nested_mean <- match.arg(nested_mean, c("wregr", "refit"))
#   set.seed(seed)
#   
#   # get the data and SFM
#   c(data, SFM, gt) %<-% gen_from_scm(sclass, n)
#   c(X, Z, W, Y) %<-% SFM
#   
#   # get the one-step debiased estimator
#   
#   # split into K folds
#   K <- 10
#   folds <- sample(x = rep(1:K, each = n / K))
#   
#   # for each fold, get out-of-fold fitted values
#   px <- y0 <- y1 <- rep(0, n)
#   y <- data[[Y]]
#   x <- data[[X]]
#   pseudo_out <- rep(NA, n)
#   for (i in seq_len(K)) {
#     
#     tst <- folds == i
#     dev <- folds %in% setdiff(seq_len(K), i)[1:6]
#     val <- folds %in% setdiff(seq_len(K), i)[7:9]
#     
#     # develop models on dev
#     mod_x_z <- cv_xgb(data[dev, Z], data[dev, X])
#     mod_x_zw <- cv_xgb(data[dev, c(Z, W)], data[dev, X])
#     mod_y_xzw <- cv_xgb(data[dev, c(X, Z, W)], data[dev, Y])
#     
#     # get the nested mean
#     if (nested_mean == "wregr") {
#       
#       px_z_val <- pred_xgb(mod_x_z, data[val, Z])
#       px_zw_val <- pred_xgb(mod_x_zw, data[val, c(Z, W)])
#       
#       weights <- ifelse(x[val] == 1, px_z_val / (1 - px_z_val) * (1 - px_zw_val) / px_zw_val,
#                         px_z_val / (1 - px_z_val))
#       
#       mod_nested <- cv_xgb(data[val, Z], data[val, Y], weights = weights)
#       ey_nest_tst <- pred_xgb(mod_nested, data[val, Z])
#     } else if (nested_mean == "refit") {
#       
#       y_tilde <- pred_xgb(mod_y_xzw, data[val, c(X, Z, W)], intervention = 1)
#       mod_nested <- cv_xgb(data[val, c(X, Z)], y_tilde)
#       ey_nest_tst <- pred_xgb(mod_nested, data[tst, c(X, Z)], intervention = 0)
#     }
#     
#     # get the pseudo-outcome on the test fold
#     ey_x1zw_tst <- pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)], intervention = 1)
#     px_zw_tst <- pred_xgb(mod_x_zw, data[tst, c(Z, W)])
#     px_z_tst <- pred_xgb(mod_x_z, data[tst, Z])
# 
#     pseudo_out[tst] <-
#       x[tst] * (y[tst] - ey_x1zw_tst) * (1 - px_zw_tst) / px_zw_tst * 1 / (1 - px_z_tst) +
#       (1-x)[tst] / (1 - px_z_tst) * (ey_x1zw_tst - ey_nest_tst) +
#       ey_nest_tst
#   }
#   
#   # compute the debiased estimator
#   psi_osd <- mean(pseudo_out)
#   lwr <- psi_osd - 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
#   upr <- psi_osd + 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
#   
#   # get the ground truth
#   gt <- compute_PO(sclass, c(0, 1), 0, 1, n = n)
#   
#   data.frame(#psi_pi = psi_pi, 
#              psi_osd = psi_osd, 
#              lwr = lwr, upr = upr, # confidence intervals
#              gt = gt, seed = seed)
# }
