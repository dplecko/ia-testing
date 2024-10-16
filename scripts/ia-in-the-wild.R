
ricu:::init_proj()

ia_inspect <- function(dataset) {
  
  p_val <- function(m, s) return(2 * (1 - pnorm(abs(m / s))))
  
  ltx_row <- function(p, dataset) {
    
    ds_names <- list(
      compas = "COMPAS",
      census = "Census 2018",
      credit = "UCI Credit",
      mimic = "MIMIC-IV",
      heloc = "HELOC",
      adult = "UCI Adult",
      wine = "UCI Wine Quality",
      bank = "UCI Bank Marketing",
      obesity = "UCI Obesity",
      diabetes = "BRFSS Diabetes"
    )
    
    paste(c(ds_names[[dataset]], ifelse(p < 0.05, "\\textbullet", "")), 
          collapse = " & ")
  }
  
  dat <- load_data(dataset)
  dat$data <- as.data.frame(dat$data)
  dat$data <- dat$data[, unlist(dat$sfm[c("X", "Z", "W", "Y")])]
  dat <- preproc_data(dat$data, dat$sfm)
  data <- dat$data
  sfm <- dat$sfm
  
  # narrow down the quantities
  target_ia <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
  # check if the outcome is binary
  binary <- length(unique(data[[sfm$Y]])) == 2
  pvals <- rep(0, length(target_ia))
  
  for (log_scale in unique(c(FALSE, binary))) {
    
    iat <- ia_tests(data, sfm$X, sfm$Z, sfm$W, sfm$Y, log_scale = log_scale)
    iat[, pval := 2 * pnorm(-abs(psi_osd / dev))]
    
    # see if there is any evidence of not rejecting
    pvals <- pmax(pvals, iat[ia %in% target_ia]$pval)
  }
  
  cat(ltx_row(pvals, dataset), " \\\\ \\hline \n")
}

datasets <- c("compas", "census", "credit", "mimic", "heloc", "adult", "wine", 
              "bank", "obesity", "diabetes")

# run the interaction testing on the datasets
for (ds in datasets) ia_inspect(ds)
