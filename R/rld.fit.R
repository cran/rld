
#data: design matrix
#if "frailty=TRUE" then use ObsLikFunNu; if "frailty=FALSE" then use ObsLikFunNoNu.

rld.fit <- function(X, C, delta, initial, lower, upper, frailty){
  
  #Observed likelihood function with frailty
  
  ObsLikFunNu <- function(par, data, C, delta){
    
    npar <- length(par)
    regpar <- par[1:(npar-1)]
    Nu <- par[npar]
    n <- length(delta)
    
    designmat <- data
    rowstart <- c()
    rowend <- c()
    
    rowstart[1] <- 1
    rowend[1] <- C[1]
    
    for(i in 2:n){
      rowstart[i] <- rowend[i-1]+1
      rowend[i] <- rowend[i-1]+C[i]
    }
    
    #create A(t) function
    Atfun <- function(regpar, i, time){
      if (time>0){
        designmati<- designmat[rowstart[i]:rowend[i],]
        j <- 1
        temp1 <- ifelse(is.null(dim(designmati)), exp(designmati%*%regpar), exp(designmati[1,]%*%regpar))
        
        while (j<time){
          temp1 <- temp1+exp(designmati[j+1,]%*%regpar)
          j <- j+1
        }
        return (temp1)
      } else if (time==0){
        return(0)
      }
    }
    
    temp2 <- c()
    
    for (i in 1:n){
      temp2[i] <- (1-delta[i])*(-1/Nu)*log(Nu*Atfun(regpar = regpar, i = i, time = C[i])+1)+
        delta[i]*log((Nu*Atfun(regpar = regpar, i = i, time = C[i]-1)+1)^(-1/Nu)-
                       ((Nu*Atfun(regpar = regpar, i = i, time = C[i])+1))^(-1/Nu))
    }
    return(sum(temp2))
  }
  
  ###################################################################################################
  
  ObsLikFunNoNu <- function(par, data, C, delta){
    
    npar <- length(par)
    n <- length(delta)
    designmat <- data
    rowstart <- c()
    rowend <- c()
    
    rowstart[1] <- 1
    rowend[1] <- C[1]
    
    for(i in 2:n){
      rowstart[i] <- rowend[i-1]+1
      rowend[i] <- rowend[i-1]+C[i]
    }
    
    #create A(t) function
    Atfun <- function(par, i, time){
      if (time>0){
        designmati<- designmat[rowstart[i]:rowend[i],]
        j <- 1
        temp1 <- ifelse(is.null(dim(designmati)), exp(designmati%*%par), exp(designmati[1,]%*%par))
        
        while (j<time){
          temp1 <- temp1+exp(designmati[j+1,]%*%par)
          j <- j+1
        }
        return (temp1)
      } else if (time==0){
        return(0)
      }
    }
    
    temp2 <- c()
    
    for (i in 1:n){
      temp2[i] <- delta[i]*log(exp(-Atfun(par = par, i = i, time = C[i]-1))-
                                 exp(-Atfun(par = par, i = i, time = C[i])))-
        (1-delta[i])*Atfun(par = par, i = i, time = C[i])
    }
    return(sum(temp2))
  }
  
  ###################################################################################################
  
  if(frailty){
    result <- optim(par = initial, fn = ObsLikFunNu, gr = NULL, method = c("L-BFGS-B"), lower = lower, 
                    upper = upper, control = list(fnscale = -1, maxit = 1000), hessian = TRUE,
                    data = X, C = C, delta = delta)
  } else{
    result <- optim(par = initial, fn = ObsLikFunNoNu, gr = NULL, method = c("L-BFGS-B"), lower = lower, 
                    upper = upper, control = list(fnscale = -1, maxit = 1000), hessian = TRUE,
                    data = X, C = C, delta = delta)
  }
  
  coefficients <- result$par
  hessian <- result$hessian
  LikFunValue <- result$value
  output <- list("coefficients" = coefficients, "hessian" = hessian, "LikFunValue" = LikFunValue)
  return(output)
}
