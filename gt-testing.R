
test_gt <- function(sclass) {
  
  res <- NULL
  for (i in 1:100) {
    
    # order Z, W, Y
    po_000 <- compute_PO(sclass, 0, 0, 0)
    po_001 <- compute_PO(sclass, 0, 0, 1)
    po_010 <- compute_PO(sclass, 0, 1, 0)
    po_011 <- compute_PO(sclass, 0, 1, 1)
    po_100 <- compute_PO(sclass, 1, 0, 0)
    po_101 <- compute_PO(sclass, 1, 0, 1)
    po_110 <- compute_PO(sclass, 1, 1, 0)
    po_111 <- compute_PO(sclass, 1, 1, 1)
    
    tese <- (po_011 - po_000) - (po_111 - po_100)
    
    deie0 <- (po_001 - po_000) - (po_011 - po_010)
    deie1 <- (po_101 - po_100) - (po_111 - po_110)
    
    dese <- (po_001 - po_000) - (po_101 - po_100)
    
    iese <- (po_010 - po_000) - (po_110 - po_100)
    
    deiese <- deie0 - deie1
    
    res <- rbind(
      res,
      data.frame(value = tese, inner_boot = i, ia = "TE x SE"),
      data.frame(value = deie0, inner_boot = i, ia = "DE x IE"),
      data.frame(value = dese, inner_boot = i, ia = "DE x SE"),
      data.frame(value = iese, inner_boot = i, ia = "IE x SE"),
      data.frame(value = deiese, inner_boot = i, ia = "DE x IE x SE")
    )
    cat("\r", i)
  }
  cat("\n")
  gt <- as.data.table(gen_from_scm(sclass, 10)$gt)
  res <- as.data.table(res)
  gt <- merge(gt, res[, list(test_stat = mean(value > 0)), by = "ia"], by = "ia")
  for (i in seq_len(nrow(gt))) {
    
    if (gt[i]$gt) {
      
      if (gt[i]$test_stat > 0.4 & gt[i]$test_stat < 0.6) {
        
        cat("SCM class", sclass, "interaction", gt[i]$ia, "may actually be FALSE.\n")
      }
    } else if (!gt[i]$gt) {
      
      if (gt[i]$test_stat < 0.4 || gt[i]$test_stat > 0.6) {
        
        cat("SCM class", sclass, "interaction", gt[i]$ia, "may actually be TRUE.\n")
      }
    }
  }
  
  return(NULL)
}

for (scl in LETTERS[1:5]) test_gt(scl)
