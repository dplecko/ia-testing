
ricu:::init_proj()

nrep <- 100
n <- 5000
sclass <- "B"
res <- c()
for (i in seq_len(nrep)) {
  
  res <- rbind(res, ia_tests(sclass, n, i, nested_mean = "refit"))
  cat("\r", i)
}
res <- as.data.table(res)

true_vals <- c()
for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) { 
  
  true_val <- compute_PO(sclass, xz, xw, xy, n = 10^6)
  true_vals <- rbind(
    true_vals,
    data.table(xz = xz, xw = xw, xy = xy, truth = true_val)
  )
}

res <- merge(res, true_vals, by = c("xz", "xw", "xy"))

res[, list(cov = mean(lwr < truth & upr > truth)), by = c("xz", "xw", "xy")]

# visualize the confidence intervals
ggplot(res, aes(x = psi_osd, y = seed)) +
  geom_point() + geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
  theme_bw() + geom_vline(aes(xintercept = truth), color = "red") +
  facet_wrap(~interaction(xz, xw, xy))

# visualize the spread
plt <- melt(res, id.vars = c("seed", "truth"))
ggplot(plt, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) + theme_bw() +
  geom_vline(xintercept = unique(plt$truth), linetype = "dashed", color = "red")
