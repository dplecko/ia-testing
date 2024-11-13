
ricu:::init_proj()
set.seed(2024)

ia_inspect <- function(dataset) {
  
  if (length(dataset) > 1) {
    
    return(
      do.call(rbind, lapply(dataset, ia_inspect))
    )
  }
  
  p_val <- function(m, s) return(2 * (1 - pnorm(abs(m / s))))
  
  dat <- load_data(dataset)
  dat$data <- as.data.frame(dat$data)
  dat$data <- dat$data[, unlist(dat$sfm[c("X", "Z", "W", "Y")])]
  dat <- preproc_data(dat$data, dat$sfm$X, dat$sfm$Z, dat$sfm$W, dat$sfm$Y)
  data <- dat$data
  sfm <- dat$sfm
  
  # narrow down the quantities
  target_ia <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
  # check if the outcome is binary
  binary <- length(unique(data[[sfm$Y]])) == 2
  pvals <- rep(0, length(target_ia))
  
  itw <- NULL
  for (log_risk in unique(c(FALSE, binary))) {
    
    iat <- one_step_debias(data, sfm$X, sfm$Z, sfm$W, sfm$Y, log_risk = log_risk)
    iat <- as.data.table(iat)
    iat <- iat[measure %in% target_ia]
    iat[, pval := 2 * pnorm(-abs(value / sd))]  
    itw <- rbind(itw, iat)
  }
  
  itw[, dataset := dataset]
  return(itw)
}

ia_itw_latex <- function(itw) {
  
  ltx_row <- function(p, n, dataset) {
    
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
    
    paste(c(ds_names[[dataset]], paste0(n), ifelse(p < 0.05, "\\textbullet", "")), 
          collapse = " & ")
  }
  
  for (dts in unique(itw$dataset)) {

    pvals <- itw[dataset == dts][, max(pval), by = c("measure", "dataset")]$V1
    n_samp <- nrow(load_data(dts)$data)
    cat(ltx_row(pvals, n_samp, dts), " \\\\ \\hline \n")
  }
}

datasets <- c("compas", "census", "credit", "mimic", "heloc", "adult", "wine", 
              "bank", "obesity", "diabetes")

# run the interaction testing on the datasets
itw <- ia_inspect(datasets)

# get the latex table rows
ia_itw_latex(itw)
