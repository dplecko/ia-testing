
Sys.setenv(DATATABLE_NUM_THREADS=4)
Sys.setenv(OMP_NUM_THREADS=4)
Sys.setenv(RANGER_XGB_CORES=1)
library(data.table)
setDTthreads(1)
library(zeallot)
library(xgboost)
library(ggplot2)
tex <- latex2exp::TeX
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))
library(ricu)
