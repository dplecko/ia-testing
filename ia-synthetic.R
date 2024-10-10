
# nohup Rscript ia-synthetic.R > parallel-synth.log 2>&1 &
ricu:::init_proj()
set.seed(2024)

scm_class <- LETTERS[1:5]
nrep <- 100
res <- c()
for (sclass in scm_class) {
  
  # for (sample_size in c(500, 750)) {
  for (sample_size in c(500, 750, 1500, 3000, 5000, 8000)) {
    
    cat("Class:", sclass, ", n = ", sample_size, "parallelized...\n")
    
    ia_chunk <- parallel::mclapply(
      # lapply(
      seq_len(nrep),
      function(rep) {
        
        c(data, SFM, gt) %<-% gen_from_scm(sclass, n = sample_size)
        c(X, Z, W, Y) %<-% SFM
        
        ret <- c()
        for (method in c("MTG")) {
          
          model <- if (method == "medDML") "ranger" else "logistic"
          ia_iter <- ia_vals(data, X = X, Z = Z, W = W, Y = Y, x0 = 0, x1 = 1, 
                             method = method, model = model,
                             ret = "values")
          ia_iter <- merge(ia_iter, gt, by = "ia") # merge-in ground truth
          ia_iter[, scm_class := sclass]
          ia_iter[, method := method]
          ia_iter[, sample_size := sample_size]
          ia_iter[, rep := rep]
          ret <- rbind(ret, ia_iter)
        }
        return(ret)
      }, mc.cores = 64 # parallel::detectCores() / 2
    )
    res <- rbind(res, do.call(rbind, ia_chunk))
  }
}

# analyze and plot
res <- as.data.table(res)
save(res, file = file.path("results", "ia-synthetic-stats-glm.RData"))
# load(file.path("results", "ia-synthetic-stats.RData"))

#' * bias correctio n / symmetry inspection * 
ggplot(
  res[gt == FALSE & method == "MTG" & inner_boot == 1 & outer_boot == 1][,
      mean(value > 0), by = c("scm_class", "method", "ia", "sample_size")],
  aes(x = V1)
) +
  geom_histogram() + theme_bw() +
  facet_grid(rows = vars(sample_size))

# 
# # zoom in on a specific instance
# sing_dat <- res[scm_class == "B" & ia == "DE x IE" &
#                   inner_boot == 1 & outer_boot == 1]
# ggplot(
#   sing_dat,
#   aes(x = value, fill = method)
# ) +
#   geom_density(alpha = 0.4) +
#   geom_vline(xintercept = 0, color = "red", linetype = "dashed") + theme_bw() +
#   facet_grid(rows = vars(sample_size)) +
#   theme(legend.position = "bottom")#+
#   #coord_cartesian(xlim = quantile(sing_dat$value, probs = c(0.025, 0.975)))
# 
# 
# res[scm_class == "D" & ia == "DE x IE" & method == "MTG" & gt == FALSE][,
#   list(lwr = quantile(value, probs = c(0.05)),
#        upr = quantile(value, probs = c(0.95))), 
#   by = c("sample_size", "rep")
# ][, mean(lwr < 0 & upr > 0), by = c("sample_size")]
# 

#' * confidence interval validity *
ret <- ci_validity(res[method == "MTG"], alpha = 0.05, bootstrap = "bc")

ggplot(
  ret[, list(cov = mean(cov)), by = c("scm_class", "ia", "sample_size", "method")],
  aes(x = cov)
) +
  geom_histogram() + theme_bw() +
  facet_grid(rows = vars(sample_size))


# res[, ia_ind := ifelse(gt, "IA", "No IA")]
# 
# # analysis I:
# ggplot(res[method == "MTG"], aes(x = pval, color = factor(scm_class), 
#                                     linetype = factor(method),
#                 linewidth = sqrt(sample_size), 
#                 group = interaction(scm_class, method, sample_size))) +
#   stat_ecdf() + theme_bw() +
#   facet_grid(rows = vars(ia_ind), cols = vars(ia)) +
#   geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed") +
#   scale_linetype_manual(values = c("solid", "dotted")) +
#   scale_linewidth_continuous(range = c(0.3, 3))
# 
# 
# # analysis II:
# alpha <- 0.05
# res[, h0_true := !gt]
# res[, reject := pval < alpha / 2]
# res[h0_true == TRUE, err := (reject == 1)]
# res[h0_true == FALSE, err := (reject == 0)]
# 
# # Type I error
# type_err <- res[, list(err = mean(err), h0_true = h0_true), 
#                 by = c("method", "scm_class", "sample_size", "ia")]
# type_err[, type := ifelse(h0_true == TRUE, "Type I", "Type II")]
# ggplot(type_err[method == "MTG"], 
#        aes(y = err, x = sample_size, color = scm_class, linetype = method)) +
#   geom_line() + theme_bw() +
#   facet_grid(rows = vars(type), cols = vars(ia)) +
#   xlab("Sample size") + ylab("Type I / Type II Error")
