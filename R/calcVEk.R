
calcVEk <- function(object, newdata, CIlevel = 0.95){

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

  if(is.list(newdata)==FALSE) stop("Argument 'newdata' should be a list")

  indxcontrast <- match("contrgroup", names(newdata), nomatch = 0)
  if(indxcontrast==0) stop("Name 'contrgroup' is required in the newdata list")

  indxref <- match("refgroup", names(newdata), nomatch = 0)
  if(indxref==0) stop("Name 'refgroup' is required in the newdata list")

  contrlevel <- newdata$contrgroup
  reflevel <- newdata$refgroup

  indxdose <- match("dose", varnames, nomatch = 0)
  fixnames <- varnames[-c(indxid, indxtime, indxdelta, indxdose)]

  namespoicontr <- match(fixnames, names(contrlevel), nomatch = 0)
  if(any(namespoicontr==0)) stop("Variable names in the newdata list should match with those in the data frame")

  namespoiref <- match(fixnames, names(reflevel), nomatch = 0)
  if(any(namespoiref==0)) stop("Variable names in the newdata list should match with those in the data frame")

  contrlevel <- contrlevel[namespoicontr]
  reflevel <- reflevel[namespoiref]

  if (length(contrlevel)>1){
    contrastgrp <- cbind(doselevels, t(replicate(ndlevel, contrlevel)))
  } else {
    contrastgrp <- cbind(doselevels, replicate(ndlevel, contrlevel))
  }
  colnames(contrastgrp) <- names(sdata)

  if (length(reflevel)>1){
    refgrp <- cbind(doselevels, t(replicate(ndlevel, reflevel)))
  } else {
    refgrp <- cbind(doselevels, replicate(ndlevel, reflevel))
  }

  colnames(refgrp) <- names(sdata)
  tempdata <- rbind(contrastgrp, refgrp, sdata)
  tempdesignmat <- model.matrix(formulaRHS, data = tempdata)
  designmat <- tempdesignmat[1:(2*ndlevel),]

  if (frailty){
    regpar <- param[-length(param)]
    Nu <- param[length(param)]

    #vaccine efficacy estimation per challenge

    tempratio <- 1-(Nu*exp(designmat%*%regpar)+1)^(-1/Nu)
    tempcontr <- tempratio[1:ndlevel]
    tempref <- tempratio[(ndlevel+1):(2*ndlevel)]
    VE <- 1-tempcontr/tempref
    names(VE) <- paste("dose", 1:ndlevel)

    #variance-covariance matrix of VEk

    contrdesign <- designmat[1:ndlevel,]
    refdesign <- designmat[(ndlevel+1):(2*ndlevel),]

    D1 <- exp(contrdesign%*%regpar)
    D0 <- exp(refdesign%*%regpar)

    Deparcoeff1 <- -(Nu*D1+1)^(-1/Nu-1)*D1
    Deparcoeff0 <- -(Nu*D0+1)^(-1/Nu-1)*D0

    DeNu1 <- ((1/Nu^2)*log(Nu*D1+1)-(1/Nu)*D1/(Nu*D1+1))*(Nu*D1+1)^(-1/Nu)
    DeNu0 <- ((1/Nu^2)*log(Nu*D0+1)-(1/Nu)*D0/(Nu*D0+1))*(Nu*D0+1)^(-1/Nu)

    Deparmat1 <- cbind(c(Deparcoeff1)*contrdesign, DeNu1)
    Deparmat0 <- cbind(c(Deparcoeff0)*refdesign, DeNu0)

    VarVE <- c()
    VarlogVE <- c()
    lwr <- c()
    upr <- c()

    for (k in 1:ndlevel){
      tempDe <- -Deparmat1[k,]*(1-(Nu*D0[k]+1)^(-1/Nu))+(1-(Nu*D1[k]+1)^(-1/Nu))*Deparmat0[k,]
      DeVEk <- -tempDe/(1-(Nu*D0[k]+1)^(-1/Nu))^2

      VarVE[k] <- t(DeVEk)%*%VarTheta%*%DeVEk
      VarlogVE[k] <- VarVE[k]/(1-VE[k])^2
      #lwr[k] <- VE[k]-critval*sqrt(VarVE[k])
      #upr[k] <- VE[k]+critval*sqrt(VarVE[k])
      lwr[k] <- log(1-VE[k])-critval*sqrt(VarlogVE[k])
      upr[k] <- log(1-VE[k])+critval*sqrt(VarlogVE[k])
    }
  }else{
    tempratio <- 1-exp(-exp(designmat%*%param))
    tempcontr <- tempratio[1:ndlevel]
    tempref <- tempratio[(ndlevel+1):(2*ndlevel)]
    VE <- 1-tempcontr/tempref
    names(VE) <- paste("dose", 1:ndlevel)

    #variance-covariance matrix of VEk

    contrdesign <- designmat[1:ndlevel,]
    refdesign <- designmat[(ndlevel+1):(2*ndlevel),]

    D1 <- exp(contrdesign%*%param)
    D0 <- exp(refdesign%*%param)

    Deparcoeff1 <- -exp(-D1)*D1
    Deparcoeff0 <- -exp(-D0)*D0

    Deparmat1 <- c(Deparcoeff1)*contrdesign
    Deparmat0 <- c(Deparcoeff0)*refdesign

    VarVE <- c()
    VarlogVE <- c()
    lwr <- c()
    upr <- c()

    for (k in 1:ndlevel){
      tempDe <- -Deparmat1[k,]*(1-exp(-D0[k]))+(1-exp(-D1[k]))*Deparmat0[k,]
      DeVEk <- -tempDe/(1-exp(-D0[k]))^2

      VarVE[k] <- t(DeVEk)%*%VarTheta%*%DeVEk
      VarlogVE[k] <- VarVE[k]/(1-VE[k])^2
      #lwr[k] <- VE[k]-critval*sqrt(VarVE[k])
      #upr[k] <- VE[k]+critval*sqrt(VarVE[k])
      lwr[k] <- log(1-VE[k])-critval*sqrt(VarlogVE[k])
      upr[k] <- log(1-VE[k])+critval*sqrt(VarlogVE[k])
    }
  }

  output <- list("VE" = VE, "se" = sqrt(VarVE), "lwr" = 1-exp(upr), "upr" = 1-exp(lwr), "ndlevel" = ndlevel)
  class(output) <- "calcVEk"

  return(output)
}


print.calcVEk <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("Perchallenge VE estimates:\n")
  print(x$VE, digits = digits)

  cat("Standard error estimates:\n")
  print(x$se, digits = digits)

  cat("Lower bound of confidence interval for log(1-VE):\n")
  print(x$lwr, digits = digits)

  cat("Upper bound of confidence interval for log(1-VE):\n")
  print(x$upr, digits = digits)
  class(x) <- "calcVEk"
  invisible(x)
}


summary.calcVEk <- function(object,...){

  VE <- object$VE
  se <- object$se
  lwr <- object$lwr
  upr <- object$upr
  ndlevel <- length(VE)
  TAB <- cbind(VE, se, lwr, upr)
  colnames(TAB) <- c("Estimate", "Std.Error", "  log(1-VE) lwr", "  log(1-VE) upr")
  rownames(TAB) <- paste("dose", 1:ndlevel, sep = "")
  res <- list("VE" = round(TAB, digits = 4))
  class(res) <- "summary.calcVEk"
  return(res)
}


print.summary.calcVEk <- function(x,...){

  cat("Per-challenge VE estimates:\n")
  printCoefmat(x$VE)
  class(x)<-"summary.calcVEk"
}
