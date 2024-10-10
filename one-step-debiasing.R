
ricu:::init_proj()

cv_xgb <- function(df, y, df_test, X = "X") {
  
  dtrain <- xgb.DMatrix(data = as.matrix(df), label = y)
  
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
  
  if (!is.null(df_test)) {
    
    xgb <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = cv$best_iteration,
      verbose = FALSE
    )
    
    if (X %in% names(df)) {
      
      df_x0 <- df_x1 <- df_test
      df_x0[[X]] <- 0
      df_x1[[X]] <- 1
      return(
        list(
          po = predict(xgb, as.matrix(df_test)),
          po0 = predict(xgb, as.matrix(df_x0)),
          po1 = predict(xgb, as.matrix(df_x1))
        )
      )
    } else {
      
      return(list(po = predict(xgb, as.matrix(df_test))))
    }
  }
  
  return(cv$pred)
}

est_yx1_x0 <- function(sclass, n, seed) {
  
  set.seed(seed)
  
  # get the data and SFM
  c(data, SFM, gt) %<-% gen_from_scm(sclass, n)
  c(X, Z, W, Y) %<-% SFM
  
  # get the \hat P estimator (plug-in)
  pos <- cv_xgb(data[, c(X, Z)], data[[Y]], data[, c(X, Z)])
  psi_pi <- (
    mean(pos[["po1"]][data[[X]] == 0]) -  mean(pos[["po0"]][data[[X]] == 0])) - (
      mean(pos[["po1"]][data[[X]] == 1]) - mean(pos[["po0"]][data[[X]] == 1])
    )
  
  # get the one-step debiased estimator
  
  # split into K folds
  K <- 5
  folds <- sample(x = rep(1:K, each = n / K))
  
  # for each fold, get out-of-fold fitted values
  px <- y0 <- y1 <- rep(0, n)
  y <- data[[Y]]
  x <- data[[X]]
  for (i in seq_len(K)) {
    
    px[folds == i] <- cv_xgb(data[folds != i, Z], data[folds != i, X],
                             data[folds == i, Z])[["po"]]
    
    y_pos <- cv_xgb(data[folds != i, c(X, Z)], data[folds != i, Y],
                    data[folds == i, c(X, Z)])
    y0[folds == i] <- y_pos[["po0"]]
    y1[folds == i] <- y_pos[["po1"]]
  }
  
  # compute the debiased estimator
  pseudo_out <- ( 
    # E[Y_{x_1} | x_0]
    (x / (1 - mean(x)) * (1 - px) / px * (y - y1) + (1 - x) / (1 - mean(x)) * y1) -
      (1 - x) / (1 - mean(x)) * y # E[Y_{x_0} | x_0]
  ) - 
    ( 
      x / mean(x) * y # E[Y_{x_1} | x_1]
       - # E[Y_{x_0} | x_1]
        ((1 - x) / mean(x) * px / (1 - px) * (y - y0) + x / mean(x) * y0)
    ) 
  
  psi_osd <- mean(pseudo_out)
  lwr <- psi_osd - 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
  upr <- psi_osd + 1.96 * sqrt(var(pseudo_out) / length(pseudo_out))
  
  # get the ground truth
  gt <- (compute_PO(sclass, 0, 1, 1, n = n) -
    compute_PO(sclass, 0, 0, 0, n = n)) - (
      compute_PO(sclass, 1, 1, 1, n = n) -
        compute_PO(sclass, 1, 0, 0, n = n)
    )
  
  data.frame(psi_pi = psi_pi, psi_osd = psi_osd, 
             lwr = lwr, upr = upr, # confidence intervals
             gt = gt, seed = seed)
}

nrep <- 100
n <- 5000
sclass <- "B"
res <- c()
for (i in seq_len(nrep)) {
  
  res <- rbind(res, est_yx1wx0(sclass, n, i, nested_mean = "refit"))
  cat("\r", i)
}
res <- as.data.table(res)
truth <- (compute_PO(sclass, 0, 1, 1, n = 10^6) -
         compute_PO(sclass, 0, 0, 0, n = 10^6)) - (
           compute_PO(sclass, 1, 1, 1, n = 10^6) -
             compute_PO(sclass, 1, 0, 0, n = 10^6)
         )

res[, truth := truth]
cat("Coverage =", mean(res$lwr < truth & res$upr > truth), "\n")

# visualize the confidence intervals
ggplot(res, aes(x = psi_osd, y = seed)) +
  geom_point() + geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
  theme_bw() + geom_vline(xintercept = unique(res$truth), color = "red")

# visualize the spread
plt <- melt(res, id.vars = c("seed", "truth"))
ggplot(plt, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) + theme_bw() +
  geom_vline(xintercept = unique(plt$truth), linetype = "dashed", color = "red")
