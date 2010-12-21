## "_helpers" don't get exported in namespace


## Pull out many coefficients from a series g(s) = \sum_{j=0}^\infty a_j s^j

power.coef.many = function(power.series, n=1024,...){
  twopii = complex(real=0,imag=2*pi)
  s.seq = exp(twopii*seq(from=0,to=n-1)/n)
  Re((1/n)*fft(power.series(s=s.seq,...)))
}

## Pull out the k-th coefficient from a series g(s) = \sum_{j=0}^\infty a_j s^j
power.coef.one = function(power.series, n=1000,k,...){
  twopii = complex(real=0,imag=2*pi)
  s.seq = exp(twopii*seq(from=0,to=n-1)/n)
  weights.seq = exp(-twopii*seq(from=0,to=n-1)*k/n)
  PSvals <- apply(as.matrix(s.seq), 1, function(s){power.series(s,...)});
                                        #mean(Re(power.series(s=s.seq,...)*weights.seq))
  mean(Re(PSvals*weights.seq))
}


timeave.laplace = function(r,s,t,lambda,mu,nu,X0) {
  myDet <- (lambda+mu+r)^2 - 4*lambda*mu
  if (myDet < 0 ){
    print("timeavelaplace: determinant<0");
    tmp <- c(r,s,t,lambda,mu,nu,X0)
    names(tmp) <- c("r","s","t","lambda","mu","nu","X0")
    print( tmp)
  }
  x1 = (lambda + mu + r - sqrt(myDet)) / (2*lambda)
  x2 = (lambda + mu + r + sqrt(myDet)) / (2*lambda)
  exp.term = exp(-lambda*(x2-x1)*t)
  ratio.term = (s - x1) / (s - x2)
  H0 = (x1 - x2*ratio.term*exp.term) / (1 - ratio.term*exp.term)
  H1 = ((x2 - x1) / ( (ratio.term*exp.term - 1)*(s-x2) ))^(nu/lambda) * exp(-nu*(1-x1)*t)
  val <- (H0)^(X0) * H1
}
##Note: H could be used in add,rem,timeave.generators but was added later.
## not sure what diff b/t a1,a2 and exparg1,exparg2 is.  room for 'r' term?
H <- function( a1, a2, exparg1, exparg2, r, s, t, i, L, nu){
  expterm <- exp(L*(exparg1-exparg2)*r*t);
  nume1 <- a1 - a2* expterm *(s-a1)/(s-a2);
  denom1 <- 1-expterm*(s-a1)/(s-a2);
  nume2 <- a1-a2;
  denom2 <- s-a2 - (s-a1)*expterm;
  (nume1/denom1)^i  *  (nume2/denom2)^(nu/L)  *  exp(nu*(r*exparg1 - 1)*t);
}



###Note: To use "hessian" function requires real arguments and real value.  
### Here we did integral over complex circle first then apply hessian.
### Still get correct result, although that is not always the case!  e.g. addhold derivs!
addrem.joint.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4){
  myAddRemGenCoef <- function(vec){  #d(addrem)/dudv evaluated at u,v both ==1. arg MUST be 's' for power.coef.one
    power.coef.one(power.series=function(s){addrem.generator(u=vec[1], v=vec[2], s=s, t=t, lambda=lambda, mu=mu, nu=nu, X0=X0)},
                   n=n, k=Xt);
  }
  return(hessianOneSided(func=myAddRemGenCoef, x=c(1,1),sides=c(-1,-1), method.args=list(d=.0001, r=r))[1,2]);
}

addhold.joint.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4){
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta);
  myGen <- function(vec,s){
    addhold.generator(u=vec[1], w=vec[2], s=s, t=t, lambda=lambda, mu=mu, nu=nu, X0=X0);
  }
  secondDerivGen <- function(s){
    tmpR <- function(vec){Re(myGen(vec,s))}
    tmpI <- function(vec){Im(myGen(vec,s))}
    if (derivTypeInfo[1] && derivTypeInfo[2]){
      require(numDeriv) ## copied genD so it accepts args;
####NOTE: delta works here as delta b/c differentiating at 1 ...
      myReal <- hessian(func=tmpR, x=c(1,0),
                        method.args=list(d=delta,eps=delta, r=r))[1,2];
      myIm <-  hessian(func=tmpI, x=c(1,0),
                       method.args=list(d=delta,eps=delta, r=r))[1,2];
    }
    else { #here min(derivTypeInfo[3:4]) <delta
      myReal <- hessianOneSided(func=tmpR, x=c(1,0), sides=c(-1,1),
                                method.args=list(d=delta, eps=delta, r=r))[1,2];
      myIm   <- hessianOneSided(func=tmpI, x=c(1,0), sides=c(-1,1),
                                method.args=list(d=delta, eps=delta, r=r))[1,2];
    }
    -complex(real=myReal,im=myIm);
  }
  joint.mean <- power.coef.one(power.series=secondDerivGen, n=n,k=Xt);
}

remhold.joint.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4){
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta);
  myGen <- function(vec,s){
    remhold.generator(v=vec[1], w=vec[2], s=s, t=t, lambda=lambda, mu=mu, nu=nu, X0=X0)
  }
  secondDerivGen <- function(s){
    tmpR <- function(vec){Re(myGen(vec,s))}
    tmpI <- function(vec){Im(myGen(vec,s))}
    if (derivTypeInfo[1] && derivTypeInfo[2]){
      require(numDeriv) ## copied genD so it accepts args;
      myReal <- hessian(func=tmpR, x=c(1,0),
                        method.args=list(d=delta,eps=delta, r=r))[1,2];
      myIm <-  hessian(func=tmpI, x=c(1,0),
                       method.args=list(d=delta,eps=delta, r=r))[1,2];
    }
    else { #here min(derivTypeInfo[3:4]) <delta ; so do onesided
      myReal <- hessianOneSided(func=tmpR, x=c(1,0), sides=c(-1,1),
                                method.args=list(d=delta,eps=delta, r=r))[1,2];
      myIm   <- hessianOneSided(func=tmpI, x=c(1,0), sides=c(-1,1),
                                method.args=list(d=delta,eps=delta,r=r))[1,2];
    }
    -complex(real=myReal,im=myIm);
  }
  joint.mean <- power.coef.one(power.series=secondDerivGen, n=n,k=Xt);
}


hold.joint.meanSq.one <- function(t,lambda,mu,nu=0,X0=1,Xt, r=4, n=1024,delta=1e-4){
  myGen <- function(w, s){hold.generator(w=w, s=s, t=t, lambda=lambda, mu=mu, nu=nu, X0=X0);}
  secondDerivGen <- function(s){
    tmp <- function(w){Re(myGen(w,s))}
    myReal <- hessianOneSided(func=tmp, x=0,sides=1,
                              method.args=list(d=delta,eps=delta, r=r));
    tmp <- function(w){Im(myGen(w,s))} 
    myIm <- hessianOneSided(func=tmp, x=0, sides=1,
                            method.args=list(d=delta,eps=delta, r=r));
    complex(re=myReal,im=myIm);
  }
  joint.mean <- power.coef.one(power.series=secondDerivGen, n=n, k=Xt);
}

## #### unused but useful to keep around as examples how to use H

## #"add.generator"
## Hplus <- function( r, s, t, i, L, m, nu){ 
##   a1 <- (L + m - sqrt((L+m)^2 - 4*L*m*r) ) / (2*L*r);
##   a2 <- (L + m + sqrt((L+m)^2 - 4*L*m*r) ) / (2*L*r);
##   H(a1,a2, a1, a2, r,s,t,i,L,nu);  
## }

## #"rem.generator"
## Hminus <- function( r, s, t, i, L, m, nu){
##   s1 <- (L + m - sqrt((L+m)^2 - 4*L*m*r) ) / (2*L);
##   s2 <- (L + m + sqrt((L+m)^2 - 4*L*m*r) ) / (2*L);
##   H(s1,s2, s1/r, s2/r, r,s,t,i,L,nu);  
## }

## Hstar <- function( r, s, t, i, L, m, nu){
##   p1 <- (L + m +r - sqrt((L+m+r)^2 - 4*L*m) ) / (2*L);
##   p2 <- (L + m +r + sqrt((L+m+r)^2 - 4*L*m) ) / (2*L);
##   H(p1,p2, p1/r, p2/r, r,s,t,i,L,nu);  
## }



                                        # Via L'Hopital's Rule
process.generator.ratesEqual <- function(s, time, lambda, mu, nu, X0){
  piece1 <- ( 1 +
             (1-s) / ( (s-1)*mu*time -1 ) ) ^ X0;
  piece2 <- (1 / (1 - mu * time *(s-1))) ^ (nu / lambda)
  piece1 * piece2;
}



##function that gets the domain of hte variables in the generating function
## and decides to do a 1sided or 2 sided deriv
### the latter two entries in the result are the values
### that can be used as deltas appropriately on both sides.
derivType <- function(L,mu, eps=1e-04){
  countsType <- TRUE;
  holdtimeType <- TRUE;
  v <- sqrt(L/mu)
  UB.counts <- (v+1/v)/2
  LB.holdtime <- - (sqrt(L) - sqrt(mu))^2;

  h.counts <- (UB.counts-1)/2;
  h.holdtime <- -LB.holdtime/2;
  if (h.counts <= eps) countsType <- FALSE;
  if (h.holdtime<=eps) holdtimeType <- FALSE;
  c(countsType,holdtimeType, h.counts,h.holdtime);
}


add.joint.mean.many = function(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024){
  return(power.coef.many(num.deriv,
                         ftn=add.generator,
                         var=1, ## differentiate add.generator at r=1
                         delta=delta,
                         n=n,
                         t=t,
                         lambda=lambda,
                         mu=mu,
                         nu=nu,
                         X0=X0)
         )
}

rem.joint.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
                               delta=0.001,n=1024,r=4){
  myMeanCoef <- function(rr){
    power.coef.one(power.series=function(s){rem.generator(r=rr,s=s,t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)},
                   n=n, k=Xt);
  }
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta)
  if (derivTypeInfo[1]){
    require(numDeriv) ## copied genD so it accepts args;
    joint.mean <- genD(myMeanCoef, x=1,
                       method.args=list(d=delta,eps=delta, r=r))$D[1]
  }
  else
    joint.mean <- genDoneSided(myMeanCoef, x=1,sides=-1,
                               method.args=list(d=delta,eps=delta, r=r))$D[1];
  return(joint.mean);
}

rem.joint.mean.many = function(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024){
  return(power.coef.many(num.deriv,
                         ftn=rem.generator,
                         var=1, ## differentiate rem.generator at r=1
                         delta=delta,
                         n=n,
                         t=t,
                         lambda=lambda,
                         mu=mu,
                         nu=nu,
                         X0=X0)
         )
}




timeave.joint.mean.many = function(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024){
  return(-1*power.coef.many(num.deriv, ## need to multiply by -1
                            ftn=timeave.laplace,
                            var=0, ## differentiate timeave.laplace at r=0
                            delta=delta,
                            n=n,
                            t=t,
                            lambda=lambda,
                            mu=mu,
                            nu=nu,
                            X0=X0)
         )
}

add.joint.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4){
  myMeanCoef <- function(r){
    power.coef.one(power.series=function(s){add.generator(r=r,s=s,t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)},
                   n=n, k=Xt);
  }
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta)
  if (derivTypeInfo[1]){
    require(numDeriv) ## copied genD so it accepts args;
    joint.mean <- genD(myMeanCoef, x=1, method.args=list(d=delta,eps=delta, r=r))$D[1]
    ## res <- try(joint.mean <- genD(myMeanCoef, x=1, method.args=list(d=delta,eps=delta, r=r))$D[1])
    ## if (class(res)=="try-error"){print(c(lambda,mu,nu,X0, Xt,delta,r))}
  }
  else
    joint.mean <- genDoneSided(myMeanCoef, x=1,sides=-1,
                               method.args=list(d=delta,eps=delta, r=r))$D[1];
  return(joint.mean);
}

#### this works the same as add.joint.meanSq.one;
#### note delta is applied to 2 different derivatives; one is richardson method, other isn't ...
add.joint.meanSq.one = function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024, r=4){
  myMeanSqCoef <- function(r){
    power.coef.one(power.series=function(s){add.generator(r=r,s=s,t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)},
                   n=n, k=Xt);
  }
  ##myFactorialMean <- num.2deriv(ftn=myMeanSqCoef, var=1, delta=delta)
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta);
  if (derivTypeInfo[1]){
    require(numDeriv) ## copied genD so it accepts args
    myFactorialMean <- hessian(func=myMeanSqCoef, x=1, method.args=list(d=delta,eps=delta,  r=r))[1,1]
  }
  else
    myFactorialMean <- hessianOneSided(func=myMeanSqCoef, x=1, sides=-1, method.args=list(d=delta, eps=delta, r=r))[1,1]
  myMean <- add.joint.mean.one(t=t, lambda=lambda, mu=mu, nu=nu, X0=X0, Xt=Xt, delta=delta, n=n);
  return(myFactorialMean+myMean);
}


#### this works the same as add.joint.meanSq.one;
#### note delta is applied to 2 different derivatives; one is richardson method, other isn't ...
rem.joint.meanSq.one = function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024, r=4){
  myMeanSqCoef <- function(r){
    power.coef.one(power.series=function(s){rem.generator(r=r,s=s,t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)},
                   n=n, k=Xt);
  }
  ##myFactorialMean <- num.2deriv(ftn=myMeanSqCoef, var=1, delta=delta)
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta);
  if (derivTypeInfo[1]){
    require(numDeriv) ## copied genD so it accepts args
    myFactorialMean <- hessian(func=myMeanSqCoef, x=1, method.args=list(d=delta,eps=delta,  r=r))[1,1]
  }
  else
    myFactorialMean <- hessianOneSided(func=myMeanSqCoef, x=1, sides=-1, method.args=list(d=delta,eps=delta,  r=r))[1,1]
  myMean <- rem.joint.mean.one(t=t, lambda=lambda, mu=mu, nu=nu, X0=X0, Xt=Xt, delta=delta, n=n);
  return(myFactorialMean+myMean);
}
rem.cond.meanSq.one = function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
  joint.mean = rem.joint.meanSq.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,delta=delta,n=n, r=r)
  trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
  prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
                     fnid="rem.cond.meanSq.one");
  ##   if(joint.mean < 10^-12){
  ##     stop("rem.cond.meanSq.one: jointmean too small");
  ##   }else{
  ##     return(joint.mean/trans.prob)
  ##   }
}


timeave.joint.mean.one = function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4){
  myMeanCoef <- function(r){
    power.coef.one(power.series=function(s){
      timeave.laplace(r=r,s=s,t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)},
                   n=n, k=Xt);
  }
  derivTypeInfo <- derivType(L=lambda,mu=mu,eps=delta)
  if (derivTypeInfo[2]){
    require(numDeriv) 
    joint.mean <- -genD(myMeanCoef, x=0,
                        method.args=list(d=0, eps=delta,r=r))$D[1]
  }
  else{
    joint.mean <- -genDoneSided(myMeanCoef, x=0,sides=1,
                                method.args=list(d=0, eps=delta,r=r))$D[1];
  }
  return(joint.mean);
}


## timeave.joint.mean.one = function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024){
##   return(-1*power.coef.one(num.deriv, ## need to multiply by -1
##                            ftn=timeave.laplace,
##                            var=0, ## differentiate timeave.laplace at r=0
##                            delta=delta,
##                            n=n,
##                            k=Xt,
##                            t=t,
##                            lambda=lambda,
##                            mu=mu,
##                            nu=nu,
##                            X0=X0)
##          )
## }


###all the *.cond.mean* functions compute joint means and divide by
### transition probabilities.  We may worry about precision errors
### if those numbers are very small.  This is the code that decides
### what to do in all those functions.
prec.error.handler <- function(joint.mean, trans.prob, prec.tol,
                               prec.fail.stop, fnid=""){
###print(fnid);
  if (trans.prob < prec.tol){
    if (prec.fail.stop==TRUE){ ##1
      ##      print("experiment");
      ##     tmp1 <- joint.mean/prec.tol;
      ##      tmp2 <- trans.prob/prec.tol;
      ##      stop( paste("first value is regular, second is multiply-by-constant.",
      ##                  joint.mean/trans.prob, tmp1/tmp2));
      stop( paste(fnid, "error: trans.prob is", trans.prob,
                  "which is too small"));
    }
    else if(prec.fail.stop==FALSE) {
      print(paste(fnid, "warning: joint.mean is", joint.mean,
                  "which is too small for comfort"));
      return(joint.mean/trans.prob)
    }
    else if (prec.fail.stop==2){
      print(paste(fnid, "error: joint.mean is", joint.mean,
                  "and trans.prob is", trans.prob,
                  ".  0 is being returned...!"))
      return(0);
    }
  }
  else return(joint.mean/trans.prob);
}

