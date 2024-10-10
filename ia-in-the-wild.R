
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
  iat <- faircause:::ia_testing(
    dat$data, X = dat$sfm$X, Z = dat$sfm$Z, W = dat$sfm$W, 
    Y = dat$sfm$Y, method = "medDML", model = "ranger",
    x0 = dat$sfm$x0, x1 = dat$sfm$x1, nboot1 = 1
  )
  
  iat <- as.data.table(iat)
  iat <- iat[, list(mean_val = mean(value), sd_val = sd(value)), by = "ia"]
  
  cat(ltx_row(p_val(iat$mean_val, iat$sd_val), dataset), " \\\\ \\hline \n")
}

datasets <- "diabetes" 
# c("compas", "census", "credit", "mimic", "heloc", "adult", "wine", "bank")

# run the interaction testing on the datasets
for (ds in datasets) ia_inspect(ds)
