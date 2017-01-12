
#Power calculation

rld.design <- function(nsim, nv, np, ndlevel, nexposure, rho, p0, RR, method = c("LRT", "log-rank"), Siglevel){

  ##data generation
    Gendata <- function(nv, np, ndlevel, nexposure, rho, p0, RR){

      if (rho>0){
      Nu <- rho/(1-rho)*(pi^2/6)
      AlphaBetap <- log((1/Nu)*(exp(-Nu*log(1-p0))-1))
      AlphaBetav <- log((1/Nu)*(exp(-Nu*log(1-RR*p0))-1))

      n <- nv+np
      W <- rgamma(n, shape = 1/Nu, rate = 1/Nu)
      R <- log(W)

      Xv <- rep(1, nv)
      Xp <- rep(0, np)
      tempX <- c(Xv, Xp)
      X <- as.matrix(sample(tempX), n, 1)

      lambda <- matrix(NA, nrow = n, ncol = ndlevel)
      for (i in 1:n){
        if (X[i]==1) {
          lambda[i,] <- 1-exp(-exp(AlphaBetav+R[i]))
        } else{
          lambda[i,] <- 1-exp(-exp(AlphaBetap+R[i]))
        }
      }
    } else if (rho==0){

        AlphaBetap <- log(-log(1-p0))
        AlphaBetav <- log(-log(1-RR*p0))

        n <- nv+np
        Xv <- rep(1, nv)
        Xp <- rep(0, np)
        tempX <- c(Xv, Xp)
        X <- as.matrix(sample(tempX), n, 1)

        lambda <- matrix(NA, nrow = n, ncol = ndlevel)
        for (i in 1:n){
          if (X[i]==1) {
            lambda[i,] <- 1-exp(-exp(AlphaBetav))
          } else{
            lambda[i,] <- 1-exp(-exp(AlphaBetap))
          }
        }
      }

      Y <- matrix(NA, nrow = n, ncol = sum(nexposure))
      for (i in 1:n){
        lambdarep <- rep(lambda[i,], nexposure)
        Y[i,] <- rbinom(n = sum(nexposure), size = 1, prob = lambdarep)
      }

      Ti <- c()
      for (i in 1:n){
        poi <- which(Y[i,]==1)
        if (length(poi)>0){
          Ti[i] <- min(poi)
        } else {
          Ti[i] <- 99
        }
      }

      delta <- c()
      delta[which(Ti<=sum(nexposure))] <- 1
      delta[which(Ti>sum(nexposure))] <- 0

      Ti[which(Ti==99)] <- sum(nexposure)

      origdata <- cbind(Ti, delta, X)
      colnames(origdata) <- c("time", "delta", "Trt")
      origdata <- data.frame(origdata)

      return(origdata)
    }

    ##do regression
    if (method=='LRT'){
      if (rho>0){
        initialval0 <- c(seq(-1, 1, length.out = ndlevel), 0.1)
        lwrb0 <- c(rep(-Inf, ndlevel), 0.01)
        uprb0 <- c(rep(Inf, ndlevel), Inf)

        initialval1 <- c(seq(-1, 1, length.out = ndlevel), 1, rep(0, ndlevel-1), 0.1)
        lwrb1 <- c(rep(-Inf, ndlevel), -Inf, rep(-Inf, ndlevel-1), 0.01)
        uprb1 <- c(rep(Inf, ndlevel), Inf, rep(Inf, ndlevel-1), Inf)

        ind <- c()

        for (i in 1:nsim){
          repeat{
            origdata <- Gendata(nv = nv, np = np, ndlevel = ndlevel, nexposure = nexposure, rho = rho,
                                p0 = p0, RR = RR)

            if (max(origdata$time)>sum(nexposure[-ndlevel])) {break}
          }

          newdata <- transdata(data = origdata, ndlevel = ndlevel, nexposure = nexposure)

          result0 <- try(rld(formula = Surv(time, delta)~factor(dose), data = newdata, initial = initialval0,
                             lower = lwrb0, upper = uprb0, frailty = TRUE))

          rldcorr1 <- function(initial){
            tempresult1 <- try(rld(formula = Surv(time, delta)~factor(dose)*factor(Trt), data = newdata,
                                   initial = initialval1, lower = lwrb1, upper = uprb1, frailty = TRUE))

            if (class(tempresult1)=='try-error'){
              cat('initial value issue, change another one.\n')
              repeat{
                newinitialval1 <- c(seq(-1, 1, length.out = ndlevel), 0.5, runif(ndlevel-1, -1, 1), 0.1)
                tempresult1 <- try(rld(formula = Surv(time, delta)~factor(dose)*factor(Trt), data = newdata,
                                       initial = newinitialval1, lower = lwrb1, upper = uprb1, frailty = TRUE))
                if (class(tempresult1)!='try-error') {break}
              }
            }
            return(tempresult1)
          }

          result1 <- rldcorr1(initialval1)

          LRTresult <- lrtest(model1 = result0, model2 = result1, TestNu = FALSE, Siglevel = Siglevel)

          if (LRTresult$pvalue<=Siglevel){
            ind[i] <- 1
          } else {
            ind[i] <- 0
          }
          print(paste("Simulation", i))
        }
      } else if (rho==0){

        initialval0 <- seq(-1, 1, length.out = ndlevel)
        lwrb0 <- rep(-Inf, ndlevel)
        uprb0 <- rep(Inf, ndlevel)

        initialval1 <- c(seq(-1, 1, length.out = ndlevel), 1, rep(0, ndlevel-1))
        lwrb1 <- c(rep(-Inf, ndlevel), -Inf, rep(-Inf, ndlevel-1))
        uprb1 <- c(rep(Inf, ndlevel), Inf, rep(Inf, ndlevel-1))

        ind <- c()

        for (i in 1:nsim){

          repeat {
            origdata <- Gendata(nv = nv, np = np, ndlevel = ndlevel, nexposure = nexposure, rho = rho,
                                p0 = p0, RR = RR)
            if (max(origdata$time)>sum(nexposure[-ndlevel])) {break}
          }

          newdata <- transdata(data = origdata, ndlevel = ndlevel, nexposure = nexposure)

          result0 <- try(rld(formula = Surv(time, delta)~factor(dose), data = newdata, initial = initialval0,
                             lower = lwrb0, upper = uprb0, frailty = FALSE))

          rldcorr1 <- function(initial){
            tempresult1 <- try(rld(formula = Surv(time, delta)~factor(dose)*factor(Trt), data = newdata,
                                   initial = initialval1, lower = lwrb1, upper = uprb1, frailty = FALSE))

            if (class(tempresult1)=='try-error'){
              cat('initial value issue, change another one.\n')
              repeat{
                newinitialval1 <- c(seq(-1, 1, length.out = ndlevel), 0.5, runif(ndlevel-1, -1, 1))
                tempresult1 <- try(rld(formula = Surv(time, delta)~factor(dose)*factor(Trt), data = newdata,
                                       initial = newinitialval1, lower = lwrb1, upper = uprb1, frailty = FALSE))
                if (class(tempresult1)!='try-error') {break}
              }
            }
            return(tempresult1)
          }

          result1 <- rldcorr1(initialval1)

          LRTresult <- lrtest(model1 = result0, model2 = result1, TestNu = FALSE, Siglevel = Siglevel)

          if (LRTresult$pvalue<=Siglevel){
            ind[i] <- 1
          } else {
            ind[i] <- 0
          }
          print(paste("Simulation", i))
        }
      }
    } else if (method=='log-rank'){

      ind <- c()
      for (i in 1:nsim){
        repeat {
          origdata <- Gendata(nv = nv, np = np, ndlevel = ndlevel, nexposure = nexposure, rho = rho,
                              p0 = p0, RR = RR)
          if (max(origdata$time)>sum(nexposure[-ndlevel])) {break}
        }

        result <- survdiff(formula = Surv(time, delta)~Trt, data = origdata, rho = 0)
        pvalue <- pchisq(q = result$chisq, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

        if (pvalue<=Siglevel){
          ind[i] <- 1
        } else {
          ind[i] <- 0
        }
        print(paste("Simulation", i))
      }
    }
    power <- sum(ind)/nsim
    output <- list(method, power)
    names(output) <- c("method", "power")

    return(output)
}


print.rld.design <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("Statistical test:\n")
  print(x$method)
  cat("Power calculation:\n")
  print(x$power, digits = digits)
  class(x) <- "rld.design"
  invisible(x)
}
