#' Rate doubly robust (RDR) estimation for weighted average treatment effect (wate)
#' @param A the treatment vector, which should be valued on 0 or 1, where 0 is the control, and 1 is the treated
#' @param Y the outcome vector
#' @param X the covariate matrix or vector (if only one covariate)
#' @param beta whether to calculate WATEs based on beta weights; if so, specify the `v1` and `v2` parameters
#' @param v1 a model parameter in beta weights; default is NA (beta weights are not estimated); can be a vector and the length must be the same as v2
#' @param v2 a model parameter in beta weights; default is NA (beta weights are not estimated); can be a vector and the length must be the same as v1
#' @param method method for estimation: "EIF" or "DML";
#'               "EIF" calculates the nuisance functions on the whole dataset, without using cross-fitting;
#'               "DML" uses cross-fitting on k folds (if use DML, the `n.folds` and `n.split` parameters need to be specified)
#' @param n.folds number of folds used for cross fitting when using DML methods; default is 5
#' @param n.split number of replications for doing the sample splits for cross-fitting when using DML methods; default is 10
#' @param return.naive whether to calculate and return the two naive estimators (the propensity score weighted and outcome imputed estimators)
#' @param ps.library method(s) used for fitting the propensity score model, chosen from the list in the SuperLearner package, and can be a vector of methods ensemble; see `SuperLearner::listWrappers()`
#' @param out.library method(s) used for fitting the outcome regression models, chosen from the list in the SuperLearner package, and can be a vector of methods ensemble; see `SuperLearner::listWrappers()`
#' @param seed seed for generating random numbers used in set.seed() function when splitting data; default is 4399
RDRwate <- function(A,
                    Y,
                    X,
                    beta=FALSE,
                    v1=NA,
                    v2=NA,
                    method="EIF",
                    n.folds=5,
                    n.split=10,
                    return.naive=FALSE,
                    ps.library=c("SL.glm", "SL.glm.interaction"),
                    out.library=c("SL.glm", "SL.glm.interaction"),
                    seed=4399) {

  X <- as.matrix(X)
  set.seed(seed=seed)
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

    if(return.naive) {
      result.naive <- .naive_wate(A=A, Y=Y, e.h=e.h, mu1.h=mu1.h, mu0.h=mu0.h,
                                  beta=beta, v1=v1, v2=v2)
      return(list(result.eif=result.eif,  result.naive=result.naive))
    } else {
      return(result.eif)
    }
  }

  #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
  #### ~~~~   The DML-based estimators  ~~~~~~ ####
  #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
  if(method=="DML") {
    n <- length(A)
    n1 <- sum(A)
    n0 <- sum(1-A)
    X <- as.data.frame(X)
    colnames(X) <- paste0("X", 1:ncol(X))
    
    seeds <- round(runif(n.split, min=0, max=10^3*n.split))
    Est.dml1 <- SE.dml1 <- Est.dml2 <- SE.dml2 <- NULL

    for(k in 1:n.split) {
      set.seed(seeds[k])
      pred.folds <- createFolds(1:n, k=n.folds, list=T)
      est.dml1 <- se.dml1 <- c()
      for(i in 1:n.folds) {
        pred.ind <- pred.folds[[i]]
        train.ind <- (1:n)[-pred.ind]
        e.h.dml2 <- mu1.h.dml2 <- mu0.h.dml2 <- A.dml2 <- Y.dml2 <- c()
        nuisance <- .RDR_nuisance(A=A[train.ind],
                                  Y=Y[train.ind],
                                  X=X[train.ind,],
                                  X.pred=X[pred.ind,],
                                  ps.method=ps.library,
                                  out.method=out.library)
        e.h <- nuisance$e.h
        mu1.h <- nuisance$mu1.h
        mu0.h <- nuisance$mu0.h

        result <- .Eif_Wate(A=A[pred.ind], Y=Y[pred.ind],
                            e.h=e.h, mu1.h=mu1.h, mu0.h=mu0.h,
                            beta=beta, v1=v1, v2=v2)
        est.dml1 <- cbind(est.dml1, result$Est)
        se.dml1 <- cbind(se.dml1, result$Std.Err)

        e.h.dml2 <- c(e.h.dml2, e.h)
        mu1.h.dml2 <- c(mu1.h.dml2, mu1.h)
        mu0.h.dml2 <- c(mu0.h.dml2, mu0.h)
        A.dml2 <- c(A.dml2, A[pred.ind])
        Y.dml2 <- c(Y.dml2, Y[pred.ind])
      }
      .dml2 <- .Eif_Wate(A=A.dml2, Y=Y.dml2,
                         e.h=e.h.dml2, mu1.h=mu1.h.dml2, mu0.h=mu0.h.dml2,
                         beta=beta, v1=v1, v2=v2)

      Est.dml1 <- rbind(Est.dml1, apply(est.dml1, 1, mean))
      Est.dml2 <- rbind(Est.dml2, .dml2$Est)

      SE.dml1 <- rbind(SE.dml1, apply(se.dml1, 1, mean)/sqrt(n.folds))
      SE.dml2 <- rbind(SE.dml2, .dml2$Std.Err/sqrt(n.folds))

      print(paste0("Sample split ", k, " done"))
    }

    result.dml.1 <- data.frame(weights=result$weights,
                               Est.mean=apply(Est.dml1, 2, mean, na.rm=T),
                               SE.mean=apply(SE.dml1, 2, mean, na.rm=T),
                               Est.median=apply(Est.dml1, 2, median, na.rm=T),
                               SE.median=apply(SE.dml1, 2, median, na.rm=T))

    result.dml.2 <- data.frame(weights=result$weights,
                               Est.mean=apply(Est.dml2, 2, mean, na.rm=T),
                               SE.mean=apply(SE.dml2, 2, mean, na.rm=T),
                               Est.median=apply(Est.dml2, 2, median, na.rm=T),
                               SE.median=apply(SE.dml2, 2, median, na.rm=T))

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

.naive_wate <- function(A,
                        Y,
                        e.h,
                        mu1.h,
                        mu0.h,
                        beta=FALSE,
                        v1=NA,
                        v2=NA) {
  weights <- c("overall", "treated", "control", "overlap", "matching", "entropy")
  tilt    <- data.frame(1, e.h, 1-e.h, e.h*(1-e.h), pmin(e.h, 1-e.h), -e.h*log(e.h)-(1-e.h)*log(1-e.h) )

  n <- length(Y)
  if(beta) {
    if(length(v1)==length(v2)) {
      weights <- c(weights, paste0("beta(", v1, ",", v2, ")"))
      tilt.beta <- matrix(NA, nrow=n, ncol=length(v1))
      for(i in 1:length(v1)) {
        tilt.beta[,i]   <- e.h^v1[i]*(1-e.h)^v2[i]
      }
      tilt    <- cbind(tilt, tilt.beta)
    }
  }

  k <- length(weights)
  naive.PS <- naive.PS.sd <- naive.OR <- naive.OR.sd <- c()
  for(i in 1:k) {
    naive.OR[i] <- mean(tilt[,i]*(mu1.h-mu0.h))/mean(tilt[,i])
    naive.OR.sd[i] <- sqrt(var(tilt[,i]*(mu1.h-mu0.h)/mean(tilt[,i])) / n)
    IPW <- A*Y/e.h - (1-A)*Y/(1-e.h)
    naive.PS[i] <- mean(tilt[,i]*IPW)/mean(tilt[,i])
    naive.PS.sd[i] <- sqrt(var(tilt[,i]*IPW/mean(tilt[,i])) / n)
  }
  naive.df <- data.frame(weights=weights,
                         naive.PS=naive.PS, naive.PS.sd=naive.PS.sd,
                         naive.OR=naive.OR, naive.OR.sd=naive.OR.sd)
  return(naive.df)
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
  EIF.df <- data.frame(weights=weights,
                       Est=EIF.est,
                       Std.Err=sqrt(colMeans(EIF^2)/n))
  return(EIF.df)
}
