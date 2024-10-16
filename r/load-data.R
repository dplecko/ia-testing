
load_data <- function(dataset) {
  
  switch(
    dataset,
    compas = load_compas(),
    census = load_census(),
    credit = load_credit(),
    mimic = load_mimic(),
    heloc = load_heloc(),
    adult = suppressMessages(load_adult()),
    wine = load_wine(),
    bank = load_bank(),
    obesity = load_obesity(),
    diabetes = load_diabetes()
  )
}

preproc_data <- function(data, SFM) {
  
  for (cmp in c("X", "Z", "W", "Y")) {
    
    
    
    for (var in SFM[[cmp]]) {
      
      tryCatch(class(data[[var]]), error = function(e) browser())
      
      if (is.numeric(data[[var]]) || is.integer(data[[var]])) next
      if (is.logical(data[[var]])) {
        
        data[[var]] <- as.integer(data[[var]])
        next
      } 
      
      if (is.character(data[[var]])) data[[var]] <- as.factor(data[[var]])
      
      # only factors at this point
      if (length(levels(data[[var]])) == 2) {
        
        data[[var]] <- as.integer(data[[var]] == levels(data[[var]])[1])
      } else {
        
        enc_mat <- model.matrix(~ ., data = data.frame(var = data[[var]]))[, -1]
        colnames(enc_mat) <- paste0(var, seq_len(ncol(enc_mat)))
        
        # remove the factor column
        SFM[[cmp]] <- setdiff(SFM[[cmp]], var)
        data[[var]] <- NULL
        
        # add the encoded data
        SFM[[cmp]] <- c(SFM[[cmp]], colnames(enc_mat))
        data <- cbind(data, enc_mat)
      }
    }
  }
  
  list(data = data, sfm = SFM)
}

load_compas <- function() {
  
  list(data = faircause::compas, sfm = faircause::SFM_proj("compas"))
}

load_census <- function() {
  
  list(data = faircause::gov_census[1:20000,], sfm = faircause::SFM_proj("census"))
}

load_credit <- function() {
  
  df <- read.csv("data/uci-credit.csv")
  
  names(df)[names(df) == "default.payment.next.month"] <- "DEFAULT"
  df$ID <- NULL
  
  sfm <- list(
    X = "SEX",
    Y = "DEFAULT",
    W = c("MARRIAGE", "EDUCATION", "LIMIT_BAL", "PAY_0", "PAY_2", "PAY_3", 
          "PAY_4", "PAY_5", "PAY_6", "BILL_AMT1", "BILL_AMT2", "BILL_AMT3", 
          "BILL_AMT4", "BILL_AMT5", "BILL_AMT6", "PAY_AMT1", "PAY_AMT2", 
          "PAY_AMT3", "PAY_AMT4", "PAY_AMT5", "PAY_AMT6"),
    Z = "AGE",
    x0 = 0, x1 = 1
  )
  
  df$SEX <- 2 - df$SEX
  
  list(data = df, sfm = sfm)
}

load_mimic <- function() {
  
  # get root
  root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
  
  # load the mortality data
  if (!file.exists(file.path(root, "data", "pci_or_death.RData"))) {
    patient_ids <- id_col(load_concepts("adm_episode", "miiv", 
                                        verbose = FALSE)[adm_episode == 1])
    dat <- load_concepts(c("acu_24", "adm_diag", "age", "sex", "charlson",
                           "lact_24", "pafi_24", "ast_24",
                           "race", "death", "los_icu"), "miiv", 
                         patient_ids = patient_ids,
                         verbose = FALSE)
    dat <- dat[race %in% c("Caucasian", "African American")]
    dat[, c(index_var(dat)) := NULL]
    dat[, pci := los_icu >= 10]
    dat[, pci_or_death := death | pci]
    dat[, c("pci", "los_icu", "death") := NULL]
    
    # remove factors
    dat[, race := as.integer(race == "Caucasian")]
    dat[, sex := (race == "Male")]
    diag_mat <- model.matrix(~ . - 1, dat[, "adm_diag"])[, -1]
    dat[, adm_diag := NULL]
    dat <- cbind(dat, diag_mat)
    
    imp_lst <- list(
      age = 65,
      acu_24 = 0,
      charlson = 0,
      lact_24 = 1,
      ast_24 = 20,
      pafi_24 = 500,
      pci_or_death = 0
    )
    
    for (i in seq_len(ncol(dat))) {
      
      var <- names(dat)[i]
      if (any(is.na(dat[[var]])) & !is.null(imp_lst[[var]]))
        dat[is.na(get(var)), c(var) := imp_lst[[var]]]
    }
    
    save(dat, diag_mat, file = file.path(root, "data", "pci_or_death.RData"))
  } else {
    
    load(file.path(root, "data", "pci_or_death.RData"))
  }
  
  dat <- dat[, c(id_vars(dat)) := NULL]
  
  sfm <- list(
    X = "race", Z = c("age", "sex"),
    W = c("acu_24", colnames(diag_mat), "charlson", "lact_24", "pafi_24", 
          "ast_24"),
    Y = "pci_or_death", x0 = 0, x1 = 1
  )
  
  list(data = dat, sfm = sfm)
}

load_heloc <- function() {
  
  root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
  df <- read.csv(file.path(root, "data", "heloc_dataset_v1.csv"))
  df$PercentTradesNeverDelq <- as.integer(df$PercentTradesNeverDelq > 90)
  
  sfm <- list(
    X = "PercentTradesNeverDelq",
    Z = c(
      "MSinceOldestTradeOpen", "AverageMInFile", # tracking period
      "NumTotalTrades", "NumTradesOpeninLast12M", "PercentInstallTrades", 
      "NumSatisfactoryTrades", "NumTrades60Ever2DerogPubRec", 
      "NumTrades90Ever2DerogPubRec", "MSinceMostRecentDelq", 
      "MaxDelq2PublicRecLast12M", "MaxDelqEver" # trading history
    ),
    W = c(
      "NetFractionRevolvingBurden", "NetFractionInstallBurden", # credit utilization
      "NumRevolvingTradesWBalance", "NumInstallTradesWBalance", 
      "NumBank2NatlTradesWHighUtilization", "PercentTradesWBalance",
      "MSinceMostRecentInqexcl7days", "NumInqLast6M", "NumInqLast6Mexcl7days" # inquiries
    ),
    Y = "ExternalRiskEstimate", x0 = 0, x1 = 1
  )
  
  list(data = df, sfm = sfm)
}

load_adult <- function() {
  
  root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
  
  adult1 <- readr::read_csv(file.path(root, "data", "adult_data.csv"),
                            col_names = FALSE)
  
  adult2 <- readr::read_csv(file.path(root, "data", "adult_test.csv"),
                            col_names = FALSE)
  
  adult <- rbind(adult1, adult2)
  adult <- data.frame(adult)
  
  colnames(adult) <- c('age', 'workclass', 'fnlwgt', 'educatoin',
                       'educatoin_num', 'marital_status', 'occupation',
                       'relationship', 'race', 'sex', 'capital_gain',
                       'capital_loss', 'hours_per_week', 'native_country',
                       'income')
  adult[, c("educatoin", "relationship", "fnlwgt",
            "capital_gain", "capital_loss")] <- NULL
  factor.columns <- c("workclass", "marital_status", "occupation", "race",
                      "sex", "native_country", "income")
  
  for(i in factor.columns) {
    
    adult[, i] <- as.factor(adult[, i])
    
  }
  
  combine_levels <- function(factor_var, levs, newLabel) {
    levels(factor_var)[levels(factor_var) %in% levs] <- newLabel
    factor_var
  }
  
  ### recode the factors appropriately
  
  # variable workclass
  adult$workclass <- combine_levels(adult$workclass,
                             levs = c("?", "Never-worked", "Without-pay"),
                             newLabel = "Other/Unknown")
  
  adult$workclass <- combine_levels(adult$workclass,
                             levs = c("Self-emp-inc", "Self-emp-not-inc"),
                             newLabel = "Self-Employed")
  
  adult$workclass <- combine_levels(adult$workclass,
                             levs = c("Federal-gov", "Local-gov",
                                      "State-gov"),
                             newLabel = "Government")
  
  # variable marital_status
  adult$marital_status <- combine_levels(adult$marital_status,
                                  levs = c("Married-AF-spouse",
                                           "Married-civ-spouse",
                                           "Married-spouse-absent"),
                                  newLabel = "Married")
  
  adult$marital_status <- combine_levels(adult$marital_status,
                                  levs = c("Divorced", "Never-married",
                                           "Separated", "Widowed"),
                                  newLabel = "Not-Married")
  
  # variable native_country
  non.US <- setdiff(levels(adult$native_country), "United-States")
  adult$native_country <- combine_levels(adult$native_country,
                                  levs = non.US,
                                  newLabel = "Not-United-States")
  
  # variable race
  adult$race <- combine_levels(adult$race,
                        levs = c("Amer-Indian-Eskimo",
                                 "Asian-Pac-Islander", "Other"),
                        newLabel = "Other")
  
  # variable income
  levels(adult$income) <- c("<=50K", "<=50K", ">50K", ">50K")
  adult$income <- as.integer(adult$income) - 1L
  adult$educatoin_num <- as.integer(adult$educatoin_num > 12)
  
  sfm <- list(
    X = "educatoin_num",
    W = c("marital_status", "workclass", "occupation",
          "hours_per_week"),
    Z = c("sex", "age", "race", "native_country"),
    Y = "income", x0 = 0, x1 = 1
  )
  
  list(data = adult, sfm = sfm)
}

load_wine <- function() {
  
  red <- read.csv("data/wine_quality/winequality-red.csv", sep = ";")
  wht <- read.csv("data/wine_quality/winequality-white.csv", sep = ";")
  red$color <- "red"
  wht$color <- "white"
  df <- rbind(red, wht)
  df$alcohol <- as.integer(df$alcohol > 10)
  
  sfm <- list(
    X = c("alcohol"),
    Z = c(
      "fixed.acidity", "volatile.acidity", "citric.acid", "chlorides", "pH",
      "color"
    ),
    W = c(
      "residual.sugar", "free.sulfur.dioxide", "total.sulfur.dioxide", 
      "sulphates",
      "density"
    ),
    Y = "quality", x0 = 0, x1 = 1
  )
  
  list(data = df, sfm = sfm)
}

load_bank <- function() {
  
  df <- read.csv("data/bank-marketing.csv", sep = ";")
  df$y <- as.integer(df$y == "yes")
  
  sfm <- list(
    X = "housing",
    Z = c(
      "age", "job", "marital", "education"
    ),
    W = c(
      "default", "balance", "loan",
      "contact", "day", "month",
      "campaign", "pdays", "previous", "poutcome"
    ),
    Y = "y", x0 = "yes", x1 = "no" 
  )
  
  list(data = df, sfm = sfm)
}

load_obesity <- function() {
  
  df <- read.csv("data/uci-obesity.csv")
  df$BMI <- df$Weight / df$Height^2
  df$Obese <- as.integer(df$BMI > 30)
  
  sfm <- list(
    X = "family_history_with_overweight",
    Z = c("Gender", "Age"),
    W = c(
      "FAVC", "FCVC", "NCP", "CAEC", "SMOKE", 
      "CH2O", "SCC", "FAF", "TUE", "CALC", "MTRANS"
    ),
    Y = "BMI", x0 = "no", x1 = "yes"
  )
  
  list(data = df, sfm = sfm)
}

load_diabetes <- function() {
  
  df <- read.csv("data/brfss-diabetes.csv")
  set.seed(2024)
  df <- df[sample(nrow(df), 20000),]
  df$Diabetes_012 <- as.integer(df$Diabetes_012 > 0)
  
  sfm <- list(
    X = c("PhysActivity"),
    Z = c("Sex", "Age", "Education", "Income", "Smoker", "Fruits", "Veggies", 
          "HvyAlcoholConsump"),
    W = c(
      "BMI", "HighBP", "HighChol", "CholCheck", "Stroke", "HeartDiseaseorAttack", 
      "DiffWalk",
      "AnyHealthcare", "NoDocbcCost",
      "GenHlth", "MentHlth", "PhysHlth"
    ),
    Y = "Diabetes_012", x0 = 0, x1 = 1
  )
  
  list(data = df, sfm = sfm)
}
