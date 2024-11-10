
gen_from_scm <- function(sclass, n, X_to_W = "random", X_to_Y = "random") {
  
  # Function to generate or fix X based on the setting
  gen_X <- function(Z, option) {
    if (option == "random") {
      return(rbinom(n, 1, prob = plogis(Z %*% c(0.3, -0.2, 0.5) + 0.2 * (Z[,1]^2))))
    } else if (option == 0) {
      return(rep(0, n))
    } else if (option == 1) {
      return(rep(1, n))
    }
  }
  
  if (sclass == "A") {
    # SCM Class A:
    Z <- matrix(rnorm(n * 3), n, 3)
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    
    # W: 3-dimensional, linear interaction of Z and X_W
    W <- Z %*% matrix(c(0.4, 0.1, -0.3, 0.2, -0.1, 0.3, 0.3, -0.2, 0.1), 3, 3) + 
      X_W * matrix(c(0.5, 0.4, 0.3), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    
    # Y: influenced by W, Z, X_Y, includes interactions of X_Y-W and X_Y-Z
    Y <- W %*% c(0.5, 0.4, 0.3) + Z %*% c(0.2, 0.1, 0.4) + X_Y * 0.7 + 
      X_Y * W[,1] * 0.2 + rnorm(n)
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )
    
  } else if (sclass == "B") {
    # SCM Class B:
    Z <- cbind(rexp(n, rate = 1), rnorm(n, 5, 1), runif(n, -2, 2))
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else X_W <- gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else X_Y <- gen_X(Z, X_to_Y)
    
    # W: cubic interaction of Z and X_W
    W <- Z %*% matrix(c(0.3, -0.5, 0.2, -0.1, 0.3, 0.1, 0.2, 0.2, -0.3), 3, 3) + 
      X_W * matrix(c(0.2, 0.1, 0.4), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    
    # Y: influenced by W, Z, X_Y, no X-W interaction but X-Y interaction
    Y <- W %*% c(0.4, 0.3, 0.2) + Z %*% c(0.1, 0.3, 0.2) + X_Y * 0.4 + rnorm(n) + 
      (W %*% c(0.1, -0.3, -0.3)) * (Z %*% c(0.1, -0.2, 0.2))
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = TRUE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )
    
  } else if (sclass == "C") {
    # SCM Class C:
    Z <- matrix(rnorm(n * 3), n, 3)
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else X_W <- gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else X_Y <- gen_X(Z, X_to_Y)
    
    # W: quadratic relationship with Z and X_W
    W <- matrix(c((Z[,1]^2) * 0.3, Z[,2] * 0.5, X_W * 0.4), n, 3) + matrix(rnorm(n * 3), n, 3)
    
    # Y: influenced by W, Z, X_Y, no X-W or X-Z interaction
    Y <- W %*% c(0.3, 0.2, 0.1) + Z %*% c(0.2, 0.1, 0.3) + rnorm(n)
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = FALSE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )
  } else if (sclass == "D") {
    # SCM Class D:
    Z <- matrix(rnorm(n * 3), n, 3)
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else X_W <- gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else X_Y <- gen_X(Z, X_to_Y)
    
    # W: 3-dimensional, influenced by Z and X_W
    W <- Z %*% matrix(c(0.4, 0.1, -0.3, 0.2, -0.1, 0.3, 0.3, -0.2, 0.1), 3, 3) + 
      X_W * matrix(c(0.5, 0.4, 0.3), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    
    # Y: includes X_Y-Z interaction but no X-W interaction
    Y <- W %*% c(0.5, 0.4, 0.3) + Z %*% c(0.2, 0.1, 0.4) + 
      X_Y * Z[,1] * 0.3 + rnorm(n)
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )
  } else if (sclass == "E") {
    # SCM Class E:
    Z <- matrix(rnorm(n * 3), n, 3)
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else X_W <- gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else X_Y <- gen_X(Z, X_to_Y)
    
    # W: 3-dimensional, quadratic relationship with Z and X_W
    W <- matrix(c((Z[,1]^2) * 0.3, Z[,2] * 0.5, X_W * 0.4), n, 3) + 
      matrix(rnorm(n * 3), n, 3)
    
    # Y: includes interaction between X_Y, Z, and W
    Y <- W %*% c(0.4, 0.3, 0.2) + Z %*% c(0.2, 0.1, 0.3) + 
      X_Y * Z[,1] * W[,3] * 0.5 + Z[,2] * W[,3] * (-0.4) + 
      X_Y * Z[, 3] + rnorm(n)
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = TRUE),
      data.frame(measure = "DE x IE x SE", gt = TRUE)
    )
  } else if (sclass == "F") {
    
    Z <- matrix(rnorm(n * 3), n, 3)
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    
    # W: linear relationship with Z and X_W
    W <- Z %*% matrix(c(0.3, -0.2, 0.4, 0.5, -0.3, 0.1, 0.2, 0.1, -0.2), 3, 3) + 
      X_W * matrix(c(0.4, 0.3, 0.5), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    
    # Y: logistic outcome influenced by W, Z, X_Y
    logit_Y <- W %*% c(0.3, 0.4, 0.2) + Z %*% c(0.1, 0.3, 0.2) + X_Y * 0.6
    prob_y <- plogis(logit_Y)
    Y <- rbinom(n, 1, prob = prob_y)
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = TRUE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )
  } else if (sclass == "G") {
    
    Z <- matrix(rnorm(n * 3), n, 3)
    
    # Generate the real X
    X <- gen_X(Z, "random")
    
    # Generate X -> W
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    
    # Generate X -> Y
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    
    # W: quadratic relationship with Z and X_W
    W <- matrix(c((Z[,1]^2) * 0.2, Z[,2] * 0.4, X_W * 0.5), n, 3) + matrix(rnorm(n * 3), n, 3)
    
    # Y: binary outcome influenced by risk ratio scale
    risk_ratio <- exp(W %*% c(0.2, 0.3, 0.1) + Z %*% c(0.2, -0.1, 0.3) + X_Y * 0.5)
    prob_y <- risk_ratio / (1 + risk_ratio)
    Y <- rbinom(n, 1, prob = prob_y)
    
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = FALSE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = TRUE)
    )
  }
  
  # Create a data frame with the variables
  data <- data.frame(X = X, Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], 
                     W1 = W[,1], W2 = W[,2], W3 = W[,3], Y = Y)
  
  # Return data, mapping of variables, and ground truth
  ret <- list(
    data = data, 
    mapping = list(X = "X", Z = c("Z1", "Z2", "Z3"), W = c("W1", "W2", "W3"), 
                   Y = "Y"), 
    gt = gt
  )
  
  if (exists("prob_y")) ret$prob_y <- prob_y
  return(ret)
}

compute_PO <- function(sclass, xz, xw, xy, n = 10^5, log_risk = FALSE) {
  
  gen <- gen_from_scm(sclass, n, X_to_W = xw, X_to_Y = xy)
  dat <- gen$data
  if (is.element(sclass, c("F", "G"))) dat$Y <- gen$prob_y
  if (log_risk) dat$Y <- log(dat$Y)
  mean(dat[dat$X %in% xz,]$Y)
}

ia_gt <- function(sclass, n = 10^5, log_risk = FALSE) {
  
  scale <- if (log_risk) "log-risk" else "difference"
  
  measures <- list(
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
    )
  )
  
  res <- NULL
  for (i in seq_along(measures)) {
    
    psi_osd <- 0
    for (j in seq_along(measures[[i]]$sgn)) {
      
      xz <- measures[[i]]$spc[[j]][1]
      xw <- measures[[i]]$spc[[j]][2]
      xy <- measures[[i]]$spc[[j]][3]
      
      # update the \hat\psi
      psi_i <- compute_PO(sclass, xz, xw, xy, n = n, log_risk = log_risk)
      # update log-scale \psi
      psi_osd <- psi_osd + measures[[i]]$sgn[j] * psi_i
    }
    
    res <- rbind(
      res,
      data.frame(gt_value = psi_osd, measure = measures[[i]]$ia, scale = scale,
                 scm_class = sclass)
    )
  }
  
  res
}
