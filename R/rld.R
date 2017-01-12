
rld <- function(formula, data, na.action, initial = NULL, lower = NULL, upper = NULL, frailty = TRUE){

  #call parameters
  sdata <- data
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  argnames <- c("formula", "data", "na.action")
  indx <- match(argnames, names(mf), nomatch = 0)
  if (indx[1] == 0) stop("A formula argument is required")
  varnames <- names(sdata)
  indxdose <- match("dose", varnames, nomatch = 0)
  if(indxdose==0) stop("Variable 'dose' is required in the dataset")
  temp <- mf[c(1, indx)]
  temp[[1]] <- as.name('model.frame')
  m <- eval(temp, parent.frame())
  Terms <- attr(m, 'terms')
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  X <- model.matrix(Terms, m)


  #data structure
  reptime <- as.numeric(Y[,1])
  repdelta <- as.numeric(Y[,2])
  repid <- rle(sdata$id)
  time <- repid$lengths
  k <- 1
  delta <- c()
  repeat{
    len <- reptime[time[k]]
    delta[k] <- repdelta[len]
    reptime <- reptime[-c(1:len)]
    repdelta <- repdelta[-c(1:len)]
    k <- k+1
    if (length(reptime)==0){break}
  }

  if(is.null(lower)){
    if(frailty){
      lower <- c(rep(-Inf, ncol(X)), 0.01)
    }else{
      lower <- rep(-Inf, ncol(X))
    }
  }

  if(is.null(upper)){
    if(frailty){
      upper <- rep(Inf, ncol(X)+1)
    }else{
      upper <- rep(Inf, ncol(X))
    }
  }

  if(is.null(initial)){
    if(frailty){
      initial <- rep(0.1, ncol(X)+1)
    }else{
      initial <- rep(0.1, ncol(X))
    }
  }

  result <- try(rld.fit(X = X, C = time, delta = delta, initial = initial, lower = lower,
                        upper = upper, frailty = frailty))

  VEexpr <- formula[-2]
  coefficients <- result$coefficients
  hessian <- result$hessian
  loglikvalue <- result$LikFunValue
  output <- list(coefficients, hessian, call, X, frailty, VEexpr, loglikvalue, data)
  names(output) <- c("coefficients", "hessian", "call", "X", "frailty", "VEexpr", "loglikvalue", "augdata")
  class(output) <- "rld"

  return(output)
}


print.rld <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  estimate <- x$coefficients
  if (x$frailty) {
    names(estimate) <- c(colnames(x$X), "\u03BD")
  } else {
    names(estimate) <- colnames(x$X)
  }
  print(estimate, digits = digits)
  class(x)<-"rld"
  invisible(x)
}


summary.rld <- function(object,...){

  coef <- object$coefficients
  se <- sqrt(diag(ginv(-object$hessian)))
  tval <- coef/se

  if (object$frailty){
    fitnonu <- rld(formula = object$call$formula, data = object$augdata,
                   frailty = FALSE)
    testnonu <- lrtest(model1 = object, model2 = fitnonu, TestNu = TRUE,
                       Siglevel = 0.05)
    LRTstat <- testnonu$statistic
    LRTpvalue <- testnonu$pvalue

    TAB1 <- cbind(coef[-length(coef)], se[-length(coef)], tval[-length(coef)],
                  2*pnorm(-abs(tval))[-length(coef)])
    colnames(TAB1) <- c("Estimate", "Std.Error", "t value", "P(>|t|)")
    rownames(TAB1) <- colnames(object$X)[-length(coef)]

    TAB2 <- cbind(coef[length(coef)], se[length(coef)], LRTstat, LRTpvalue)
    colnames(TAB2) <- c("Estimate", "Std.Error", "chisq", "p-value")
    rownames(TAB2) <- c("\u03BD")

    res <- list("call" = object$call, "coefficients1" =  round(TAB1, digits = 4),
                "coefficients2" =  round(TAB2, digits = 4), "frailty" = object$frailty)
  }else {
    TAB1 <- cbind(coef, se, tval, 2*pnorm(-abs(tval)))
    colnames(TAB1) <- c("Estimate", "Std.Error", "t value", "P(>|t|)")
    rownames(TAB1) <- colnames(object$X)

    res <- list("call" = object$call, "coefficients1" = round(TAB1, digits = 4), "frailty" = object$frailty)
  }
  class(res) <- "summary.rld"
  return(res)
}


print.summary.rld <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Coefficients:\n")
  if(x$frailty){
    printCoefmat(x$coefficients1)
    cat("-----------------------------\n")
    printCoefmat(x$coefficients2)
  }else{
    printCoefmat(x$coefficients1)
  }
    cat("-----------------------------\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  class(x)<-"summary.rld"
  invisible(x)
}
