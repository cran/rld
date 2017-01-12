
calcVEt <- function(object, nexposure, newdata, CIlevel = 0.95){

  if(class(object)!="rld") stop("'object' class should be 'rld' ")

  param <- object$coefficients
  formulaRHS <- object$VEexpr
  data <- object$augdata
  frailty <- object$frailty
  Hessian <- object$hessian
  ndlevel <- max(object$augdata$dose)

  if(length(nexposure)!=ndlevel) stop("The length of 'nexposure' should match with the number of dose levels")

  VarTheta <- ginv(-Hessian)
  critval <- abs(qnorm((1-CIlevel)/2))

  varnames <- names(data)

  indxtime <- match("time", varnames, nomatch = 0)
  indxdelta <- match("delta", varnames, nomatch = 0)
  indxid <- match("id", varnames, nomatch = 0)

  sdata <- data[,-c(indxid, indxtime, indxdelta)]

  doselevels <- c(1:ndlevel)
  expolevels <- c(1:sum(nexposure))
  repdose <- rep(doselevels, nexposure)

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
    contrastgrp <- cbind(repdose, t(replicate(sum(nexposure), contrlevel)))
  } else {
    contrastgrp <- cbind(repdose, replicate(sum(nexposure), contrlevel))
  }
  colnames(contrastgrp) <- names(sdata)

  if (length(reflevel)>1){
    refgrp <- cbind(repdose, t(replicate(sum(nexposure), reflevel)))
  } else {
    refgrp <- cbind(repdose, replicate(sum(nexposure), reflevel))
  }
  colnames(refgrp) <- names(sdata)

  tempdata <- rbind(contrastgrp, refgrp, sdata)
  tempdesignmat <- model.matrix(formulaRHS, data = tempdata)
  designmat <- tempdesignmat[1:(2*sum(nexposure)),]

  designmatcontr <- designmat[1:sum(nexposure),]
  designmatref <- designmat[(sum(nexposure)+1):(2*sum(nexposure)),]

  if(frailty){
    npar <- length(param)
    regpar <- param[1:(npar-1)]
    Nu <- param[npar]

    #A(t) function for contrast group

    Atfuncontr <- function(time){
      j <- 1
      temp1 <- exp(designmatcontr[1,]%*%regpar)
      temp2 <- exp(designmatcontr[1,]%*%regpar)*designmatcontr[1,]

      while(j<time){
        temp1 <- temp1+exp(designmatcontr[j+1,]%*%regpar)
        temp2 <- temp2+exp(designmatcontr[j+1,]%*%regpar)*designmatcontr[j+1,]
        j <- j+1
      }
      output <- list("sumexp" = temp1, "sumexpdesign" = temp2)
      return(output)
    }

    #A(t) function for reference group

    Atfunref <- function(time){
      j <- 1
      temp1 <- exp(designmatref[1,]%*%regpar)
      temp2 <- exp(designmatref[1,]%*%regpar)*designmatref[1,]

      while(j<time){
        temp1 <- temp1+exp(designmatref[j+1,]%*%regpar)
        temp2 <- temp2+exp(designmatref[j+1,]%*%regpar)*designmatref[j+1,]
        j <- j+1
      }
      output <- list("sumexp" = temp1, "sumexpdesign" = temp2)
      return(output)
    }

    #calculate VE(t)

    VEt <- c()

    for (j in expolevels){
      sumexpcontr <- Atfuncontr(time = j)$sumexp
      sumexpref <- Atfunref(time = j)$sumexp

      tempcontr <- 1-(Nu*sumexpcontr+1)^(-1/Nu)
      tempref <- 1-(Nu*sumexpref+1)^(-1/Nu)

      VEt[j] <- 1-tempcontr/tempref
    }

    #variance-covariance matrix of VEt
    VarVE <- c()
    lwr <- c()
    upr <- c()

    for (j in expolevels){
      sumexpcontr <- Atfuncontr(time = j)$sumexp
      sumexpref <- Atfunref(time = j)$sumexp

      sumexpdesigncontr <- Atfuncontr(time = j)$sumexpdesign
      sumexpdesignref <- Atfunref(time = j)$sumexpdesign

      Deparcoeff1 <- -(Nu*sumexpcontr+1)^(-1/Nu-1)*sumexpdesigncontr
      Deparcoeff0 <- -(Nu*sumexpref+1)^(-1/Nu-1)*sumexpdesignref

      DeNu1 <- ((1/Nu^2)*log(Nu*sumexpcontr+1)-(1/Nu)*sumexpcontr/(Nu*sumexpcontr+1))*(Nu*sumexpcontr+1)^(-1/Nu)
      DeNu0 <- ((1/Nu^2)*log(Nu*sumexpref+1)-(1/Nu)*sumexpref/(Nu*sumexpref+1))*(Nu*sumexpref+1)^(-1/Nu)

      Deparmat1 <- matrix(c(Deparcoeff1, DeNu1), nrow = length(param), ncol = 1)
      Deparmat0 <- matrix(c(Deparcoeff0, DeNu0), nrow = length(param), ncol = 1)

      tempDe <- -Deparmat1*(1-(Nu*c(sumexpref)+1)^(-1/Nu))+(1-(Nu*c(sumexpcontr)+1)^(-1/Nu))*Deparmat0
      DeVEk <- -tempDe/(1-(Nu*c(sumexpref)+1)^(-1/Nu))^2
      VarVE[j] <- t(DeVEk)%*%VarTheta%*%DeVEk
      lwr[j] <- VEt[j]-critval*sqrt(VarVE[j])
      upr[j] <- VEt[j]+critval*sqrt(VarVE[j])
    }
  }else{

    Atfuncontr <- function(time){
      j <- 1
      temp1 <- exp(designmatcontr[1,]%*%param)
      temp2 <- exp(designmatcontr[1,]%*%param)*designmatcontr[1,]

      while(j<time){
        temp1 <- temp1+exp(designmatcontr[j+1,]%*%param)
        temp2 <- temp2+exp(designmatcontr[j+1,]%*%param)*designmatcontr[j+1,]
        j <- j+1
      }
      output <- list("sumexp" = temp1, "sumexpdesign" = temp2)
      return(output)
    }

    #A(t) function for reference group

    Atfunref <- function(time){
      j <- 1
      temp1 <- exp(designmatref[1,]%*%param)
      temp2 <- exp(designmatref[1,]%*%param)*designmatref[1,]

      while(j<time){
        temp1 <- temp1+exp(designmatref[j+1,]%*%param)
        temp2 <- temp2+exp(designmatref[j+1,]%*%param)*designmatref[j+1,]
        j <- j+1
      }
      output <- list("sumexp" = temp1, "sumexpdesign" = temp2)
      return(output)
    }

    #calculate VE(t)

    VEt <- c()

    for (j in expolevels){
      sumexpcontr <- Atfuncontr(time = j)$sumexp
      sumexpref <- Atfunref(time = j)$sumexp

      tempcontr <- 1-(exp(-sumexpcontr))
      tempref <- 1-(exp(-sumexpref))

      VEt[j] <- 1-tempcontr/tempref
    }

    #variance-covariance matrix of VEt
    VarVE <- c()
    lwr <- c()
    upr <- c()

    for (j in expolevels){
      sumexpcontr <- Atfuncontr(time = j)$sumexp
      sumexpref <- Atfunref(time = j)$sumexp

      sumexpdesigncontr <- Atfuncontr(time = j)$sumexpdesign
      sumexpdesignref <- Atfunref(time = j)$sumexpdesign

      Deparcoeff1 <- -exp(-sumexpcontr)*sumexpdesigncontr
      Deparcoeff0 <- -exp(-sumexpref)*sumexpdesignref

      Deparmat1 <- matrix(Deparcoeff1, nrow = length(param), ncol = 1)
      Deparmat0 <- matrix(Deparcoeff0, nrow = length(param), ncol = 1)

      tempDe <- -Deparmat1*(1-exp(-c(sumexpref)))+(1-exp(-c(sumexpcontr)))*Deparmat0
      DeVEk <- -tempDe/(1-exp(-c(sumexpref)))^2
      VarVE[j] <- t(DeVEk)%*%VarTheta%*%DeVEk
      lwr[j] <- VEt[j]-critval*sqrt(VarVE[j])
      upr[j] <- VEt[j]+critval*sqrt(VarVE[j])
    }
  }
  output <- list("VE" = VEt, "se" = sqrt(VarVE), "lwr" = lwr, "upr" = upr, "time" = expolevels, "doselevel" = repdose)
  class(output) <- "calcVEt"

  return(output)
}


print.calcVEt <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("VE estimates before or at the time of challenge t:\n")
  print(x$VE, digits = digits)

  cat("Standard error estimates:\n")
  print(x$se, digits = digits)

  cat("Lower bound of confidence interval:\n")
  print(x$lwr, digits = digits)

  cat("Upper bound of confidence interval:\n")
  print(x$upr, digits = digits)

  cat("challenge time:\n")
  print(x$time)
  class(x) <- "calcVEt"
}


summary.calcVEt <- function(object,...){

  VE <- object$VE
  se <- object$se
  lwr <- object$lwr
  upr <- object$upr
  time <- object$time
  dose <- object$doselevel
  TAB <- cbind(Dose = dose, Estimate = VE, Std.Error = se, lwr = lwr, upr = upr)
  rownames(TAB) <- paste("time", 1:length(time), sep = "")
  res <- list("VE" = round(TAB, digits = 4))
  class(res) <- "summary.calcVEt"
  return(res)
}


print.summary.calcVEt <- function(x,...){

  cat("VE estimates before or at the time of challenge t:\n")
  printCoefmat(x$VE)
  class(x)<-"summary.calcVEt"
}
