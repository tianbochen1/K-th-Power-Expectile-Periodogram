kep = function(y, f, tau, k, n.cores=1) {
  generalized_expectile_regression <- function(X, y, alpha = 0.5, k = 2) {
    stopifnot(length(y) == nrow(X))
    stopifnot(k >= 1 && k <= 2)
    stopifnot(alpha > 0 && alpha < 1)
    X <- as.matrix(cbind(1, X))  # 第一列是截距
    n <- nrow(X)
    p <- ncol(X)
    loss_fn <- function(beta) {
      r <- y - X %*% beta
      weights <- ifelse(r >= 0, alpha, 1 - alpha)
      sum(weights * abs(r)^k)
    }
    
    beta_init <- rep(0, p)
    
    opt <- optim(par = beta_init, fn = loss_fn, method = "BFGS")
    return(list(coefficients = opt$par,
                convergence = opt$convergence,
                value = opt$value))
  }
  ep.parallel<-function(yy,ff,tau,tt) {
    # parallel computation of quantile harmonic regression
    # for a vector of tau (and single value of ff)
    # used when intercept = T
    ntau = length(tau)
    coef = matrix(0, 2, ntau)
    for(i in 1:ntau){
      taui = tau[i]
      x1 = cos(2*pi*ff*tt)
      x2 = sin(2*pi*ff*tt)
      o = rep(length(x1))
      X = cbind(x1,x2)
      fit <- generalized_expectile_regression(X, yy, alpha=taui, k = k)
      coef[,i] = fit$coefficients[-1]
    }
    return(coef)
  }
  
  ns<-length(y)
  nf<-length(f)
  tt<-c(1:ns)
  ntau<-length(tau)
  
  if(n.cores>1) {
    library(foreach) 
    library(doParallel) # libraries required for parallel computing
    classify.cl <- makeCluster(n.cores)
    registerDoParallel(classify.cl) 
  }
  
  
  yy<-y
  coef<-array(NA,dim=c(ntau,2,nf))
  if(n.cores>1) {
    tmp<-foreach(i=1:nf) %dopar% { ep.parallel(yy,f[i],tau=tau,tt) }
  } 
  else {
    library(foreach) 
    library(doParallel) # libraries required for parallel computing
    tmp<-foreach(i=1:nf) %do% { ep.parallel(yy,f[i],tau=tau,tt) }
  }
  
  
  out<-lapply(tmp,FUN=function(x) {apply(x^2,2,sum)})
  out<-matrix(unlist(out),ncol=ntau,byrow=T)
  out<-out*ns/4
  
  if(n.cores>1) stopCluster(classify.cl)
  
  out[out<0]<-0
  if(ntau==1) out<-c(out)
  out
  
}

expectile_peri_ls = function(y, f, tau, n.cores=1) {
  ep.parallel1<-function(yy,ff,tau,tt) {
    ntau = length(tau)
    coef = matrix(0, 2, ntau)
    for(i in 1:ntau){
      taui = tau[i]
      x1 = cos(2*pi*ff*tt)
      x2 = sin(2*pi*ff*tt)
      o = rep(1,length(x1))
      X = cbind(o,x1,x2)
      fit <- expectreg.ls(yy ~ X[,2]+ X[,3], mstop=100, expectiles=taui, quietly=TRUE)
      coef[,i] = c(fit$coefficients[[1]],fit$coefficients[[2]])
    }
    return(coef)
  }
  
  ns<-length(y)
  nf<-length(f)
  tt<-c(1:ns)
  ntau<-length(tau)
  
  yy=y
  tmp=list()
  for(i in 1:length(f)){
    tmp[[i]] = ep.parallel1(yy,f[i],tau,tt)
  }
  
  out<-lapply(tmp,FUN=function(x) {apply(x^2,2,sum)})
  out<-matrix(unlist(out),ncol=ntau,byrow=T)
  out<-out*ns/4
  
  
  out[out<0]<-0
  if(ntau==1) out<-c(out)
  out
  
}

lap.spec.new2 = function(y, f, tau=0.5, intercept=F, type=1, weights=NULL, method="fn", n.cores=1) {
  # type 1: squared L2 norm of coefficients
  # type 2: cost difference
  
  rsoid2 <- function(n, f, a, b) {
    p <- length(f)
    tt <- c(1:n)
    one <- rep(1, n)
    tmp <- (one %o% a) * cos(2*pi*(tt %o% f)) + (one %o% b) * sin(2*pi*(tt %o% f))
    tmp <- apply(tmp,1,sum)
    tmp
  }
  
  
  lap.cost<-function(y,tau=0.5,weights=NULL) {
    # cost function of quantile regression
    ns<-length(y)
    if(is.null(weights)) weights<-rep(1,ns)
    tmp<-tau*y
    sel<-which(y < 0)
    if(length(sel) > 0) tmp[sel]<-(tau-1)*y[sel]
    sum(tmp*weights,na.rm=T)
  }
  
  
  qh.parallel<-function(yy,ff,tau,tt,ns,weights,type) {
    # parallel computation of quantile harmonic regression
    # for a single value of tau (and single value of ff)
    # used when intercept = F
    if(ff==0.5) {
      fit<-try(rq(yy ~ 0+cos(2*pi*ff*tt),method=method,tau=tau,weights=weights),silent=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-c(0,0)
      } else {
        fit$coefficients<-c(fit$coefficients,0)
      }
    }
    if(ff==0) {
      fit<-NULL
      fit$coefficients<-c(0,0)
    }
    if(ff != 0.5 & ff>0) {
      fit<-try(rq(yy ~ 0+cos(2*pi*ff*tt)+sin(2*pi*ff*tt),method=method,tau=tau,weights=weights),silent=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-c(0,0)
      }
    }
    fit$residuals<-yy-rsoid2(ns,ff,fit$coefficients[1],fit$coefficients[2])
    if(type==1) tmp.coef<-fit$coefficients
    if(type==2) tmp.cost<-lap.cost(fit$residuals,tau=tau,weights=weights)
    rm(fit)
    if(type==1) return(tmp.coef)
    if(type==2) return(tmp.cost)
  }
  
  
  qh.parallel2<-function(yy,ff,tau,tt,ns,weights,type) {
    # parallel computation of quantile harmonic regression
    # for a vector of tau (and single value of ff)
    # used when intercept = T
    if(ff==0.5) {
      fit<-try(rq(yy ~ cos(2*pi*ff*tt),method=method,tau=tau,weights=weights),silient=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-rbind(quantile(yy,probs=tau),rep(0,length(tau)),rep(0,length(tau)))
      } else {
        fit$coefficients<-rbind(fit$coefficients,rep(0,length(tau)))
      }
    }
    if(ff==0) {
      fit<-NULL
      fit$coefficients<-rbind(quantile(yy,probs=tau),rep(0,length(tau)),rep(0,length(tau)))
    }
    if(ff != 0.5 & ff > 0) {
      fit<-try(rq(yy ~ cos(2*pi*ff*tt)+sin(2*pi*ff*tt),method=method,tau=tau,weights=weights),silent=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-rbind(quantile(yy,probs=tau),rep(0,length(tau)),rep(0,length(tau)))
      }
    }
    fit$coefficients<-matrix(fit$coefficients,ncol=length(tau))
    if(type==1) tmp.coef<-matrix(fit$coefficients[-1,],ncol=length(tau))
    if(type==2) {
      if(length(tau)==1) fit$coefficients<-matrix(fit$coefficients,ncol=1)
      tmp.cost<-rep(NA,length(tau))
      for(i.tau in c(1:length(tau))) {
        tmp.resid<-yy-fit$coefficients[1,i.tau]-rsoid2(ns,ff,fit$coefficients[2,i.tau],fit$coefficients[3,i.tau])
        tmp.cost[i.tau]<-lap.cost(tmp.resid,tau=tau[i.tau],weights=weights)
      }
      rm(tmp.resid)
    }
    rm(fit)
    if(type==1) return(tmp.coef)
    if(type==2) return(tmp.cost)
  }
  
  
  ns<-length(y)
  nf<-length(f)
  tt<-c(1:ns)
  ntau<-length(tau)
  
  if(n.cores>1) {
    library(foreach) 
    library(doParallel) # libraries required for parallel computing
    classify.cl <- makeCluster(n.cores)
    registerDoParallel(classify.cl) 
  }
  
  
  if(is.null(weights)) weights<-rep(1,ns)
  
  
  yy<-y
  if(intercept) {
    
    if(type==1) coef<-array(NA,dim=c(ntau,2,nf))
    if(type==2) { 
      cost<-matrix(NA,ntau,nf)
      fit<-rq(yy ~ 1,method=method,tau=tau,weights=weights)
      if(ntau==1) fit$coefficients<-matrix(fit$coefficients,ncol=1)
      cost0<-rep(NA,ntau)
      for(i.tau in c(1:ntau)) cost0[i.tau]<-lap.cost(yy-fit$coefficients[,i.tau],tau=tau[i.tau],weights=weights)
    }
    
    if(n.cores>1) {
      tmp<-foreach(i=1:nf,.packages="quantreg") %dopar% { qh.parallel2(yy,f[i],tau=tau,tt,ns,weights,type) }
    } else {
      library(foreach)
      tmp<-foreach(i=1:nf,.packages="quantreg") %do% { qh.parallel2(yy,f[i],tau=tau,tt,ns,weights,type) }
    }
    if(type==1) {
      out<-lapply(tmp,FUN=function(x) {apply(x^2,2,sum)})
      out<-matrix(unlist(out),ncol=ntau,byrow=T)
      out<-out*ns/4
    }
    if(type==2) {
      out<-matrix(unlist(tmp),ncol=ntau,byrow=T)
      out<-matrix(rep(cost0,nf),ncol=ntau,byrow=T)-out
    }
    
  } else {
    
    out<-NULL
    for(i.tau in c(1:ntau)) {
      if(type==1) coef<-matrix(NA,ncol=2,nrow=nf)
      if(type==2) yy<-y-quantile(y,probs=tau[i.tau],na.rm=T)
      cost0<-lap.cost(yy,tau=tau[i.tau],weights=weights)
      if(n.cores>1) {
        tmp<-foreach(i=1:nf,.packages="quantreg") %dopar% { qh.parallel(yy,f[i],tau=tau[i.tau],tt,ns,weights,type) }
      } else {
        library(foreach)
        tmp<-foreach(i=1:nf,.packages="quantreg") %do% { qh.parallel(yy,f[i],tau=tau[i.tau],tt,ns,weights,type) }
      }
      if(type==1) {
        coef<-matrix(unlist(tmp),ncol=2,byrow=T)
        out<-cbind(out,apply(coef^2,1,sum)*ns/4)
      }
      if(type==2) {
        cost<-c(unlist(tmp))
        out<-cbind(out,cost0-cost)
      }
    }
    
  }
  
  if(n.cores>1) stopCluster(classify.cl)
  
  out[out<0]<-0
  if(ntau==1) out<-c(out)
  out
  
}

spec.normalize = function(qper) {
  if(is.matrix(qper)) {
    return( apply(qper,2,FUN=function(x) { x/sum(x) }) )
  } else {
    return( qper/sum(qper) )
  }
}

generalized_expectile <- function(y, tau = 0.5, k = 1.5) {
  stopifnot(k > 1, k < 2)
  loss_fn <- function(theta) {
    u <- y - theta
    weights <- ifelse(u >= 0, tau, 1 - tau)
    sum(weights * abs(u)^k)
  }
  opt <- optim(par = mean(y), fn = loss_fn, method = "BFGS")
  return(opt$par)
}

pfish = function(x){
  q = length(x)
  g = q * max(x)/sum(x)
  js = 0:q
  s = 1-js*g/q
  for(i in 1:(q+1)){
    if (s[i] <0) s[i] = 0
  }
  pp = 1 - sum((-1)^js * choose(q,js) * s^(q-1))
  return(pp)
}

pgram = function(y){
  n =length(y)
  J = floor((n+1)/2)-1
  pg = (1/T)*(abs(fft(y)[1:J]))^2
  return(pg)
  }

mymix = function(n){
  x1 = as.numeric(arima.sim(list(order=c(1,0,0),ar=c(0.8),sd=1),n))
  x2 = as.numeric(arima.sim(list(order=c(1,0,0),ar=c(-0.75),sd=1),n))
  x3 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(0,-0.81),sd=1),n))
  w1<-0.75
  w2<-0.1
  y1<-0.8
  y2<--0.8
  y12<-0
  y22<--0.4
  w12<-0.0
  w22<-0.5
  wt<-psi.wts(x1,y1,y2,w1,w2)
  x0<-wt*x2+(1-wt)*x1
  wt2<-psi.wts(x0,y12,y22,w12,w22)
  x<-wt2*x3+(1-wt2)*x0
  return(x)
}

generalized_expectile_regression <- function(X, y, alpha = 0.5, k = 2) {
  stopifnot(length(y) == nrow(X))
  stopifnot(k >= 1 && k <= 2)
  stopifnot(alpha > 0 && alpha < 1)
  X <- as.matrix(cbind(1, X))  # 第一列是截距
  n <- nrow(X)
  p <- ncol(X)
  loss_fn <- function(beta) {
    r <- y - X %*% beta
    weights <- ifelse(r >= 0, alpha, 1 - alpha)
    sum(weights * abs(r)^k)
  }
  
  beta_init <- rep(0, p)
  
  opt <- optim(par = beta_init, fn = loss_fn, method = "BFGS",control = list(reltol=0.001,maxit=100))
  return(list(coefficients = opt$par,
              convergence = opt$convergence,
              value = opt$value))
}

kasec = function(x, alpha, k){
  n = length(x)
  mu = generalized_expectile_regression(rep(0, n), x, alpha, k)$coefficients[1]
  x2 = x - mu
  
  x3 = rep(0, n)
  for(i in 1:n){
    if(x2[i] >= 0){x3[i] = alpha*k*x2[i]^(k-1)}
    if(x2[i] < 0){x3[i] = -(1-alpha)*k*(-x2[i])^(k-1)}
  }
 return(x3)
}