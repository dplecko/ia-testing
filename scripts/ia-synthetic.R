
# nohup Rscript scripts/ia-synthetic.R > parallel-synth.log 2>&1 &
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
        
        ia_iter <- ia_tests(data, X, Z, W, Y, nested_mean = "refit")
        ia_iter <- merge(ia_iter, gt, by = "ia") # merge-in ground truth
        ia_iter[, scm_class := sclass]
        ia_iter[, method := "one-step"]
        ia_iter[, sample_size := sample_size]
        ia_iter[, rep := rep]
        return(ia_iter)
      }, mc.cores = parallel::detectCores() / 2
    )
    res <- rbind(res, do.call(rbind, ia_chunk))
  }
}

res <- as.data.table(res)
# save(res, file = file.path("results", "ia-synthetic-stats-one-step.RData"))
# load(file.path("results", "ia-synthetic-stats-one-step.RData"))

#' * symmetry *
ggplot(
  res[gt == FALSE][, mean(psi_osd > 0), by = c("scm_class", "ia", "sample_size")],
  aes(x = V1)
) +
  geom_histogram() + theme_bw() +
  facet_grid(rows = vars(sample_size))


#' * confidence interval validity *
res <- res[, cov := lwr < 0 & upr > 0]
ggplot(
  res[gt == FALSE, list(cov = mean(cov)), 
      by = c("scm_class", "ia", "sample_size", "method")],
  aes(x = cov)
) +
  geom_histogram() + theme_bw() +
  facet_grid(rows = vars(sample_size))


#' * Analysis I: distribution of p-values * 
res[, pval := 2 * pnorm(-abs(psi_osd / dev))]
res[, ia_ind := ifelse(gt, "Interaction", "No Interaction")]
res[, ia_ind := factor(ia_ind, levels = c("No Interaction", "Interaction"))]
ggplot(res, aes(x = pval, color = factor(scm_class),
                linewidth = sqrt(sample_size),
                group = interaction(scm_class, sample_size))) +
  stat_ecdf() + theme_bw() +
  facet_grid(rows = vars(ia_ind), cols = vars(ia)) +
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed") +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_linewidth_continuous(range = c(0.3, 3),
                             name = tex("\\sqrt{sample size}")) +
  scale_color_discrete(
    name = "SCM", labels = sapply(1:5, function(i) tex(paste0("$M_", i, "$")))
  ) +
  ylab("Empirical Cumulative Distribution Function") +
  xlab("p-value") +
  theme(
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("results/p-value-distr.png", width = 10, height = 5)

#' * Analysis II: type I and II error rates *
alpha <- 0.05
res[, h0_true := !gt]
res[, reject := pval < alpha / 2]
res[h0_true == TRUE, err := (reject == 1)]
res[h0_true == FALSE, err := (reject == 0)]

type_err <- res[, list(err = mean(err), h0_true = h0_true),
                by = c("scm_class", "sample_size", "ia")]
type_err[, type := ifelse(h0_true == TRUE, "Type II", "Type I")]
type_err[, type := factor(type, levels = c("Type II", "Type I"))]
ggplot(type_err,
       aes(y = err, x = sample_size, color = scm_class)) +
  geom_line() + theme_bw() +
  facet_grid(rows = vars(type), cols = vars(ia),
             scales = "free") +
  xlab("Sample size") + ylab("Testing Error") +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(data = subset(type_err, type == "Type II"), 
             aes(yintercept = 0.05), linetype = "dashed", color = "gray") +
  scale_color_discrete(
    name = "SCM", labels = sapply(1:5, function(i) tex(paste0("$M_", i, "$")))
  ) +
  theme(
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("results/test-errors.png", width = 10, height = 5)
