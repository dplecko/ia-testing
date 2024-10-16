
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
  
  # cat("Dataset", dataset, "\n")
  # 
  # cat(
  #   vapply(data, class, character(1L)), "\n"
  # )
  # 
  # cat("Data has", ncol(data), "names; SFM has", 
  #     length(unlist(sfm[c("X", "Z", "W", "Y")])), "\n")
  # 
  # if (setequal(names(data), unlist(sfm[c("X", "Z", "W", "Y")]))) {
  #   
  #   cat("Names match\n")
  #   
  # } else cat("Names do not match\n")
  
  # run the testing at difference scale
  
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
  
  # iat <- faircause:::ia_testing(
  #   dat$data, X = dat$sfm$X, Z = dat$sfm$Z, W = dat$sfm$W, 
  #   Y = dat$sfm$Y, method = "medDML", model = "ranger",
  #   x0 = dat$sfm$x0, x1 = dat$sfm$x1, nboot1 = 1
  # )
  # 
  # iat <- as.data.table(iat)
  # iat <- iat[, list(mean_val = mean(value), sd_val = sd(value)), by = "ia"]
  # 
  cat(ltx_row(pvals, dataset), " \\\\ \\hline \n")
}

# dataset <- "census"
# dts <- load_data(dataset)
# data <- preproc_data(dts$data, dts$sfm)
# 
# datasets <- "diabetes"
datasets <- c("compas", "census", "credit", "mimic", "heloc", "adult", "wine", 
              "bank", "obesity", "diabetes")

# run the interaction testing on the datasets
for (ds in datasets) ia_inspect(ds)
