
#Likelihood ratio test

lrtest <- function(model1, model2, TestNu = TRUE, Siglevel = 0.05){

  L1 <- model1$loglikvalue
  L2 <- model2$loglikvalue
  dof <- abs(length(model1$coefficients)-length(model2$coefficients))

  if (L1<L2) {
    L1new <- L1
    L2new <- L2} else {
      L1new <- L2
      L2new <- L1
    }
  LR <- (-2)*(L1new-L2new)

  Nucutoff <- qchibarsq(q = 1-Siglevel, df = 1, mix = 0.5)
  regparcutoff <- qchisq(p = 1-Siglevel, df = dof, ncp = 0, lower.tail = TRUE, log.p = FALSE)

  if(TestNu){
    temp <- ifelse(LR>Nucutoff, 1, 0)
  }else{
    temp <- ifelse(LR>regparcutoff, 1, 0)
  }

  pvalue <- ifelse(TestNu, pchibarsq(p = LR, df = 1, mix = 0.5, lower.tail = FALSE, log.p = FALSE),
                   pchisq(q = LR, df = dof, lower.tail = FALSE, log.p = FALSE))

  output <- list("statistic"= LR, "df" = dof, "pvalue" = pvalue, "frailty" = TestNu)
  class(output) <- "lrtest"

  return(output)
}


print.lrtest <- function(x, digits = max(3, getOption("digits") - 3), ...){
  if (x$frailty){
    cat("Likelihood Ratio test for frailty variance nu\n")
  }else{
    cat("Likelihood Ratio test for regression parameters\n")
  }
  cat("Test statistic:\n")
  print(x$statistic, digits = digits)

  cat("Degree of freedom:\n")
  print(x$df)

  cat("p-value:\n")
  print(x$pvalue, digits = digits)
  class(x) <- "lrtest"
  invisible(x)
}
