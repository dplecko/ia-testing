
ia_vals <- function(data, X, Z, W, Y, x0, x1, method = c("medDML", "MTG"),
                    model, tune_params = FALSE,
                    ret = c("p-values", "values"),
                    nboot1 = 5L, nboot2 = 100L, ...) {
  
  ret <- match.arg(ret, c("p-values", "values"))
  method <- match.arg(method, c("medDML", "MTG"))
  
  y <- as.numeric(data[[Y]]) - is.factor(data[[Y]]) # need to check (!)
  
  # handle empty Z/W sets
  Z <- if (length(Z) == 0 | identical(Z, "")) NULL else Z
  W <- if (length(W) == 0 | identical(W, "")) NULL else W
  
  if (method == "medDML") {
    
    # coerce X to a factor with levels x0, x1
    data[[X]] <- factor(data[[X]], levels = c(x0, x1))
    
    params <- res <- list()
    
    for (rep in seq_len(nboot1)) {
      
      res[[rep]] <- ia_vals_mdml(data, X, Z, W, Y, x0, x1, method, model = model, 
                                 rep, nboot = nboot2, tune_params = tune_params,
                                 params = params)
      
      params <- attr(res[[rep]], "params")
      if (rep == 1) {
        pw <- attr(res[[rep]], "pw")
      }
      
      res[[rep]]$outer_boot <- rep
    }
    
    res <- do.call(rbind, res)
  } else if (method == "MTG") {
    
    res <- ia_vals_mtg(data, X, Z, W, Y, pw_model = model, nboot = nboot1)
  }
  
  res <- as.data.table(res)
  if (ret == "values") return(res)
  res[, list(pval = 2 * min(mean(value > 0), mean(value < 0))), by = "ia"]
}

