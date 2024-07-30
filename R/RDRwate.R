#' Rate doubly robust (RDR) estimation for weighted average treatment effect (WATE)
#' @param A treatment vector; should be 0 or 1 valued, where 0 is control, and 1 is treated
#' @param Y outcome vector
#' @param X covariate matrix
#' @param beta whether to calculate WATE based on beta weights; if so, specify the below v1 & v2 parameters
#' @param v1 model parameter in beta weights; default is NA (beta weights are not estimated); can be a vector and length must be the same as v2
#' @param v2 model parameter in beta weights; default is NA (beta weights are not estimated); can be a vector and length must be the same as v1
#' @param method method for estimation: "EIF" or "DML";
#'               "EIF" calculates the nuisance functions on the whole dataset, without using cross-fitting;
#'               "DML" uses cross-fitting on k folds (k = n.folds specified below)
#' @param n.folds number of folds used for cross fitting; default is 5
#' @param ps.library method(s) used for fitting the propensity score model;
#'                   it is chosen from the list of the SuperLearner package, and can be a vector of methods (the ensemble library)
#' @param out.library method(s) used for fitting the outcome regression models;
#'                    it is chosen from the list of the SuperLearner package, and can be a vector of methods (the ensemble library)
#' @param seed seed for generating random numbers used in set.seed() function when splitting data; default is 1
RDRwate <- function(A,
                    Y,
                    X,
                    beta=FALSE,
                    v1=NA,
                    v2=NA,
                    method="EIF",
                    n.folds=5,
                    ps.library=c("SL.glm", "SL.glm.interaction"),
                    out.library=c("SL.glm", "SL.glm.interaction"),
                    seed=1) {

  # library(SuperLearner)
  # library(caret)
  # library(tidyverse)

  #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
  #### ~~~~   The EIF-based estimators  ~~~~~~ ####
  #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
  if(method=="EIF") {
    X <- as.data.frame(X)
    nuisance <- .RDR_nuisance(A=A, Y=Y, X=X, X.pred=X,
                              ps.method=ps.library, out.method=out.library)
    e.h <- nuisance$e.h
    mu1.h <- nuisance$mu1.h
    mu0.h <- nuisance$mu0.h

    result.eif <- .Eif_Wate(A=A, Y=Y, e.h=e.h, mu1.h=mu1.h, mu0.h=mu0.h,
                            beta=beta, v1=v1, v2=v2)
    return(result.eif)
  }

  #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
  #### ~~~~   The DML-based estimators  ~~~~~~ ####
  #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
  if(method=="DML") {
    n <- length(A)
    n1 <- sum(A)
    n0 <- sum(1-A)

    set.seed(seed=seed)
    pred.folds <- createFolds(1:n, k=n.folds, list=T)
    est.dml1 <- se.dml1 <- c()
    e.h.dml2 <- mu1.h.dml2 <- mu0.h.dml2 <- A.dml2 <- Y.dml2 <- c()
    for(i in 1:n.folds) {
      pred.ind <- pred.folds[[i]]
      train.ind <- (1:n)[-pred.ind]

      nuisance <- .RDR_nuisance(A=A[train.ind],
                                Y=Y[train.ind],
                                X=X[train.ind,],
                                X.pred=X[pred.ind,],
                                ps.method=ps.library,
                                out.method=out.library)
      e.h <- nuisance$e.h
      mu1.h <- nuisance$mu1.h
      mu0.h <- nuisance$mu0.h

      result <- .Eif_Wate(A=A[pred.ind], Y=Y[pred.ind], e.h=e.h, mu1.h=mu1.h, mu0.h=mu0.h,
                          beta=beta, v1=v1, v2=v2)
      est.dml1 <- cbind(est.dml1, result$Est)
      se.dml1 <- cbind(se.dml1, result$Std.Err)

      e.h.dml2 <- c(e.h.dml2, e.h)
      mu1.h.dml2 <- c(mu1.h.dml2, mu1.h)
      mu0.h.dml2 <- c(mu0.h.dml2, mu0.h)
      A.dml2 <- c(A.dml2, A[pred.ind])
      Y.dml2 <- c(Y.dml2, Y[pred.ind])
    }
    result.dml.1 <- data.frame(weights=result$weights,
                               Est=apply(est.dml1, 1, mean),
                               Std.Err=apply(se.dml1, 1, mean)/sqrt(n.folds))

    .dml2 <- .Eif_Wate(A=A.dml2, Y=Y.dml2, e.h=e.h.dml2, mu1.h=mu1.h.dml2, mu0.h=mu0.h.dml2,
                        beta=beta, v1=v1, v2=v2)
    result.dml.2 <- data.frame(weights=.dml2$weights,
                               Est=.dml2$Est,
                               Std.Err=.dml2$Std.Err)

    return(list(result.dml.1=result.dml.1, result.dml.2=result.dml.2))
  }
}

.RDR_nuisance <- function(A,
                          Y,
                          X,
                          X.pred,
                          ps.method,
                          out.method) {

  ps.fit = SuperLearner(Y=A, X=X, family=binomial(), SL.library=ps.method)
  e.h = predict(ps.fit, X.pred)$pred

  mu1.fit = SuperLearner(Y=Y[A==1], X=X[A==1,], SL.library=out.method)
  mu1.h = predict(mu1.fit, X.pred)$pred
  mu0.fit = SuperLearner(Y=Y[A==0], X=X[A==0,], SL.library=out.method)
  mu0.h = predict(mu0.fit, X.pred)$pred

  return(list(e.h=e.h, mu1.h=mu1.h, mu0.h=mu0.h))
}

.Eif_Wate <- function(A,
                      Y,
                      e.h,
                      mu1.h,
                      mu0.h,
                      beta=FALSE,
                      v1=NA,
                      v2=NA) {

  weights <- c("overall", "treated", "control", "overlap", "matching", "entropy")
  tilt    <- data.frame(1, e.h, 1-e.h, e.h*(1-e.h), pmin(e.h, 1-e.h),      -e.h*log(e.h)-(1-e.h)*log(1-e.h) )
  d.tilt  <- data.frame(0, 1,   -1,    1-2*e.h,     I(e.h<0.5)-I(e.h>0.5), log(1-e.h)-log(e.h) )

  n <- length(A)
  if(beta) {
    if(length(v1)==length(v2)) {
      weights <- c(weights, paste0("beta(", v1, ",", v2, ")"))
      tilt.beta <- d.tilt.beta <- matrix(NA, nrow=n, ncol=length(v1))
      for(i in 1:length(v1)) {
        tilt.beta[,i]   <- e.h^v1[i]*(1-e.h)^v2[i]
        d.tilt.beta[,i] <- e.h^(v1[i]-1)*(1-e.h)^(v2[i]-1) * (v1[i]*(1-e.h) + v2[i]*e.h)
      }
      tilt    <- cbind(tilt, tilt.beta)
      d.tilt  <- cbind(d.tilt, d.tilt.beta)
    }
  }
  tau.X <- mu1.h - mu0.h
  psi <- (A*(Y-mu1.h))/e.h - ((1-A)*(Y-mu0.h))/(1-e.h) + tau.X

  k <- length(weights)
  EIF <- matrix(NA, nrow=n, ncol=k)
  EIF.est <- c()
  for(i in 1:k) {
    phi.num <- tilt[,i]*psi + d.tilt[,i]*tau.X*(A-e.h)
    phi.den <- tilt[,i] + d.tilt[,i]*(A-e.h)
    EIF.est[i] <- mean(phi.num)/mean(phi.den)
    EIF[,i] <- tilt[,i]/mean(tilt[,i]) * (psi-EIF.est[i]) + d.tilt[,i]/mean(tilt[,i]) * (tau.X-EIF.est[i])*(A-e.h)
  }
  EIF.df <- data.frame(weights=weights, Est=EIF.est, Std.Err=sqrt(colMeans(EIF^2)/n))
  return(EIF.df)
}
