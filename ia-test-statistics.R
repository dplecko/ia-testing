
ricu:::init_proj()
sclass <- "D"
sample_size <- 5000
ia_tests <- c()
for (i in 1:10) {
  
  set.seed(i)
  c(data, SFM, gt) %<-% gen_from_scm(sclass, n = sample_size)
  c(X, Z, W, Y) %<-% SFM
  iter <- ia_vals(data, X = X, Z = Z, W = W, Y = Y, x0 = 0, x1 = 1, 
                  method = method, model = model, ret = "values")
  iter[, seed := i]
  ia_tests <- rbind(ia_tests, iter)
  cat("\r", i)
}

ggplot(
  ia_tests[ia == "DE x IE"],
  aes(x = value, fill = factor(seed)
      )
) +
  geom_density(alpha = 0.4) +
  theme_bw()
