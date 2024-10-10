
ci_validity <- function(res, alpha, bootstrap = c("percent", "bc")) {
  
  cov_ind <- function(value, bootstrap) {
    
    if (bootstrap == "percent") {
      
      c(lwr, upr) %<-% quantile(value, probs = c(alpha, 1 - alpha)) 
    } else if (bootstrap == "bc") {
      
      val <- value[1]
      z0 <- qnorm(mean(value[-1] < val[1]))
      z_alpha <- qnorm(alpha)
      z_1malpha <- qnorm(1 - alpha)
      alpha_p <- pnorm(2 * z0 + z_alpha)
      inv_alpha_p <- pnorm(2 * z0 + z_1malpha)
      
      c(lwr, upr) %<-% quantile(value[-1], probs = c(alpha_p, inv_alpha_p)) 
    }
    
    lwr < 0 & upr > 0
  }
  
  # use only non-interactions
  ret <- res[gt == FALSE]
  
  # define the by_args
  by_args <- c("scm_class", "ia", "sample_size", "method", "rep")
  
  ret <- ret[, list(cov = cov_ind(value, bootstrap)), by = by_args]

  ret
}
