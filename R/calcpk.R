
calcpk <- function(object, predlevel, CIlevel = 0.95){

  if(class(object)!="rld") stop("'object' class should be 'rld' ")

  param <- object$coefficients
  formulaRHS <- object$VEexpr
  data <- object$augdata
  frailty <- object$frailty
  Hessian <- object$hessian
  ndlevel <- max(object$augdata$dose)

  VarTheta <- ginv(-Hessian)
  critval <- abs(qnorm((1-CIlevel)/2))

  varnames <- names(data)
  indxtime <- match("time", varnames, nomatch = 0)
  indxdelta <- match("delta", varnames, nomatch = 0)
  indxid <- match("id", varnames, nomatch = 0)
  sdata <- data[,-c(indxid, indxtime, indxdelta)]

  doselevels <- c(1:ndlevel)

  indxdose <- match("dose", varnames, nomatch = 0)
  fixnames <- varnames[-c(indxid, indxtime, indxdelta, indxdose)]

  namespoi <- match(fixnames, names(predlevel), nomatch = 0)
  if(any(namespoi==0)) stop("Variable names in the vector should match with those in the data frame")

  predlevel <- predlevel[namespoi]

  if (length(predlevel)>1){
    specgrp <- cbind(doselevels, t(replicate(ndlevel, predlevel)))
  } else {
    specgrp <- cbind(doselevels, replicate(ndlevel, predlevel))
  }
  colnames(specgrp) <- names(sdata)

  tempdata <- rbind(specgrp, sdata)
  tempdesignmat <- model.matrix(formulaRHS, data = tempdata)
  designmat <- tempdesignmat[1:ndlevel,]

  if(frailty){
    regpar <- param[-length(param)]
    Nu <- param[length(param)]

    #vaccine efficacy estimation per challenge

    temprob <- 1-(Nu*exp(designmat%*%regpar)+1)^(-1/Nu)
    tempspecgrp <- temprob[1:ndlevel]
    names(tempspecgrp) <- paste("Vaccine_pk", 1:ndlevel)

    #variance-covariance matrix of pk

    D1 <- exp(designmat%*%regpar)
    Deparcoeff1 <- -(Nu*D1+1)^(-1/Nu-1)*D1
    DeNu1 <- ((1/Nu^2)*log(Nu*D1+1)-(1/Nu)*D1/(Nu*D1+1))*(Nu*D1+1)^(-1/Nu)
    Deparmat1 <- cbind(c(Deparcoeff1)*designmat, DeNu1)

    Varp1 <- c()
    Varlogp1 <- c()
    lwr1 <- c()
    upr1 <- c()

    for (k in 1:ndlevel){
      Dep1k <- -Deparmat1[k,]
      Varp1[k] <- t(Dep1k)%*%VarTheta%*%Dep1k
      Varlogp1[k] <- (tempspecgrp[k]*(1-tempspecgrp[k]))^(-2)*Varp1[k]
      lwr1[k] <- log(tempspecgrp[k]/(1-tempspecgrp[k]))-critval*sqrt(Varlogp1[k])
      upr1[k] <- log(tempspecgrp[k]/(1-tempspecgrp[k]))+critval*sqrt(Varlogp1[k])
    }
  }else{
    #vaccine efficacy estimation per challenge

    tempprob <- 1-exp(-exp(designmat%*%param))
    tempspecgrp <- tempprob[1:ndlevel]
    names(tempspecgrp) <- paste("Vaccine_pk", 1:ndlevel)

    #variance-covariance matrix of pk

    D1 <- exp(designmat%*%param)
    Deparcoeff1 <- -exp(-D1)*D1
    Deparmat1 <- c(Deparcoeff1)*designmat

    Varp1 <- c()
    Varlogp1 <- c()
    lwr1 <- c()
    upr1 <- c()

    for (k in 1:ndlevel){
      Dep1k <- -Deparmat1[k,]
      Varp1[k] <- t(Dep1k)%*%VarTheta%*%Dep1k
      Varlogp1[k] <- (tempspecgrp[k]*(1-tempspecgrp[k]))^(-2)*Varp1[k]
      lwr1[k] <- log(tempspecgrp[k]/(1-tempspecgrp[k]))-critval*sqrt(Varlogp1[k])
      upr1[k] <- log(tempspecgrp[k]/(1-tempspecgrp[k]))+critval*sqrt(Varlogp1[k])
    }
  }
  output <- list("pk" = tempspecgrp, "pk_se" = sqrt(Varp1), "lwr" = exp(lwr1)/(1+exp(lwr1)),
                 "upr" = exp(upr1)/(1+exp(upr1)))

  class(output) <- "calcpk"
  return(output)
}


print.calcpk <- function(x, digits = max(3, getOption("digits") - 3),...){

  cat("Per-challenge probability of infection:\n")
  print(c(x$pk), digits = digits)

  cat("Standard error estimates:\n")
  print(c(x$pk_se), digits = digits)

  cat("Lower bound of confidence interval for log odds pk:\n")
  print(c(x$lwr), digits = digits)

  cat("Upper bound of confidence interval for log odds pk:\n")
  print(c(x$upr), digits = digits)

  class(x) <- "calcpk"
  invisible(x)
}


summary.calcpk <- function(object,...){

  pk <- object$pk
  Vaccine_se <- object$pk_se
  Vaccine_lwr <- object$lwr
  Vaccine_upr <- object$upr
  ndlevel <- length(pk)
  TAB <- cbind(pk, Vaccine_se, Vaccine_lwr, Vaccine_upr)
  colnames(TAB) <- c("Estimate", "Std.Error", " lwr(log odds)", " upr(log odds)")
  rownames(TAB) <- paste("dose", 1:ndlevel, sep = "")

  res <- list("pkest" = round(TAB, digits = 4))
  class(res) <- "summary.calcpk"
  return(res)
}


print.summary.calcpk <- function(x,...){

  cat("Per-challenge probability of infection estimates (pk):\n")
  printCoefmat(x$pkest)
  class(x) <- "summary.calcpk"
  invisible(x)
}
