
tex <- function(x) {
  
  if (length(x) == 1) return(latex2exp::TeX(x))
  
  sapply(x, function(xi) latex2exp::TeX(xi))
}
