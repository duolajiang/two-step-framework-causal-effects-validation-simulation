library(matrixcalc) ## for function: is.positive.semi.definite
library(mvtnorm) ## for function: pmvnorm
library(prodlim) ## for function: rowSums

CopulaTestIntegral <- function(){
  library(MASS)
  library(mvtnorm)
  m.a <- 0.2; m.b <- 0.3; m.c <- 0.4; rho <- 0
  x.a <- qnorm(m.a); x.b <- qnorm(m.b); x.c <- qnorm(m.c)
  mu <- rep(0, 3)
  Sigma <- matrix(c(1,rho,rho,
                    rho,1,rho,
                    rho,rho,1), nrow = 3, ncol = 3)
  
  p000 <- pmvnorm(lower=c(x.a,x.b,x.c), upper=c(Inf,Inf,Inf), mean = mu,corr = Sigma)
  p001 <- pmvnorm(lower=c(x.a,x.b,-Inf), upper=c(Inf,Inf,x.c), mean = mu,corr = Sigma)
  p010 <- pmvnorm(lower=c(x.a,-Inf,x.c), upper=c(Inf,x.b,Inf), mean = mu,corr = Sigma)
  p011 <- pmvnorm(lower=c(x.a,-Inf,-Inf), upper=c(Inf,x.b,x.c), mean = mu,corr = Sigma)
  
  p100 <- pmvnorm(lower=c(-Inf,x.b,x.c), upper=c(x.a,Inf,Inf), mean = mu,corr = Sigma)
  p101 <- pmvnorm(lower=c(-Inf,x.b,-Inf), upper=c(x.a,Inf,x.c), mean = mu,corr = Sigma)
  p110 <- pmvnorm(lower=c(-Inf,-Inf,x.c), upper=c(x.a,x.b,Inf), mean = mu,corr = Sigma)
  p111 <- pmvnorm(lower=c(-Inf,-Inf,-Inf), upper=c(x.a,x.b,x.c), mean = mu,corr = Sigma)
  
  p1__ <- pmvnorm(lower=c(-Inf,-Inf,-Inf), upper=c(x.a,Inf,Inf), mean = mu,corr = Sigma)
  
  
  p.a <- 0.2; p.b <- 0.3; p.c <- 0.4; rho.ab <- 0.2; rho.ac <- 0.3; rho.bc <- 0.4
  x.a <- qnorm(p.a); x.b <- qnorm(p.b); x.c <- qnorm(p.c)
  Sigma.ab <- matrix(c(1,rho.ab,
                       rho.ab,1), nrow = 2, ncol= 2)
  Sigma.bc <- matrix(c(1,rho.bc,
                       rho.bc,1), nrow = 2, ncol= 2)
  
  Sigma <- matrix(c(1,rho.ab,rho.ac,
                    rho.ab,1,rho.bc,
                    rho.ac,rho.bc,1), nrow = 3, ncol = 3)
  F.a11 <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(x.a,Inf,Inf),mean = rep(0,3),corr = Sigma)
  F.1b1 <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(Inf,x.b,Inf),mean = rep(0,3),corr = Sigma)
  F.11c <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(Inf,Inf,x.c),mean = rep(0,3),corr = Sigma)
  F.ab1 <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(x.a,x.b,Inf),mean = rep(0,3),corr = Sigma)
  F.a1c <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(x.a,Inf,x.c),mean = rep(0,3),corr = Sigma)
  F.1bc <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(Inf,x.b,x.c),mean = rep(0,3),corr = Sigma)
  F.abc <- pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(x.a,x.b,x.c),mean = rep(0,3),corr = Sigma)
  
  p000 <- F.abc; c(p000,F.abc) 
  p001 <- F.ab1 - F.abc; c(p001,pmvnorm(lower = c(-Inf,-Inf,x.c),upper = c(x.a,x.b,Inf),mean = rep(0,3),corr = Sigma))
  p010 <- F.a1c - F.abc; c(p010,pmvnorm(lower = c(-Inf,x.b,-Inf),upper = c(x.a,Inf,x.c),mean = rep(0,3),corr = Sigma))
  p100 <- F.1bc - F.abc;  c(p100,pmvnorm(lower = c(x.a,-Inf,-Inf),upper = c(Inf,x.b,x.c),mean = rep(0,3),corr = Sigma))
  p110 <- F.11c - F.a1c - F.1bc + F.abc; c(p110,pmvnorm(lower = c(x.a,x.b,-Inf),upper = c(Inf,Inf,x.c),mean = rep(0,3),corr = Sigma))
  p101 <- F.1b1 - F.ab1 - F.1bc + F.abc; c(p101,pmvnorm(lower = c(x.a,-Inf,x.c),upper = c(Inf,x.b,Inf),mean = rep(0,3),corr = Sigma))
  p011 <- F.a11 - F.ab1 - F.a1c + F.abc; c(p011,pmvnorm(lower = c(-Inf,x.b,x.c),upper = c(x.a,Inf,Inf),mean = rep(0,3),corr = Sigma))
  p111 <- 1-p000-p001-p010-p100-p110-p101-p011; c(p111,pmvnorm(lower = c(x.a,x.b,x.c),upper = c(Inf,Inf,Inf),mean = rep(0,3),corr = Sigma))
  
  
  library(MASS)
  x <- mvrnorm(mu=mu,Sigma = Sigma,n=1e6)
  x <- pnorm(x)
  x.a <- ifelse(x[,1]<p.a,0,1)
  x.b <- ifelse(x[,2]<p.b,0,1)
  x.c <- ifelse(x[,3]<p.c,0,1)
  xx <- cbind(x.a,x.b,x.c)
  mean((xx[,1]==1)&(xx[,2]==1)&(xx[,3]==0))
}


bart.m.hat <- function(data,xvars,ntree,ndpost){
  x1 <- data[data$Z==1,xvars]
  x0 <- data[data$Z==0,xvars]
  y1 <- data[data$Z==1,'Y']
  y0 <- data[data$Z==0,'Y'] 
  test <- data[,xvars]
  
  post_mean.1 <- wbart(x1,y1,test,ntree,ndpost=ndpost)
  post_mean.1$prob.test <- post_mean.1$yhat.test
  post_mean.1$prob.test.mean <- apply(post_mean.1$prob.test,2,mean)
  m1.hat <- post_mean.1$prob.test.mean
  
  post_mean.0 <- wbart(x0,y0,test,ntree,ndpost=ndpost)
  post_mean.0$prob.test <- post_mean.0$yhat.test
  post_mean.0$prob.test.mean <- apply(post_mean.0$prob.test,2,mean)
  m0.hat <- post_mean.0$prob.test.mean
  
  return(cbind(m1.hat,m0.hat))
}

potential_outcomes_direct_ipw_dr <- function(data,xvars,ntree,ndpost){
  set.seed(123)
  ntree <- ntree
  ndpost <- ndpost
  #############################################################
  direct_formula_main <- c('Z',xvars)
  direct_formula_main <- paste(direct_formula_main, collapse = '+' )
  direct_formula_interaction <- paste(xvars, collapse = ':Z +' )
  direct_formula <- paste(direct_formula_main,'+',direct_formula_interaction,':Z')
  direct_formula <- as.formula(paste('Y~(',direct_formula,')'))
  
  ps_formula <- paste(xvars,sep = '', collapse = '+')
  ps_formula <- as.formula(paste('Z~(',ps_formula,')^2'))
  
  dr_formula <- paste(xvars,sep = '', collapse = '+')
  dr_formula <- as.formula(paste('Y~(',dr_formula,')^2'))
  # =============================================================================
  # propensity score 
  # #1) LReg ps
  covariates <- model.matrix(ps_formula,data)[,-1]
  y <- as.matrix(data[,'Z'])
  cv.lasso.m <- cv.glmnet(covariates,y,alpha=1,family='gaussian',type.measure="mse")
  data$lrps <- predict(cv.lasso.m,covariates,s='lambda.min',type='response')
  
  # 2) BART ps
  x <- data[,xvars]
  y <- data[,'Z']
  post <- pbart(x,y,x,ntree = ntree, ndpost = ndpost)
  post$prob.test <- pnorm(post$yhat.test)
  post$prob.test.mean <- apply(post$prob.test,2,mean)
  data$bartps <- post$prob.test.mean
  
  
  # =============================================================================
  # caculate: 1) E(Y=1|X,Z) : joint mean model in direct method; 
  #           2) E(Y=1|X,Z=1),E(Y=1|X,Z=0) : separate mean model in doubly robust method E(Y=1|X,Z=1),E(Y=1|X,Z=0)
  #           3) E(Y=1|X,Z,PS(BART)) : joint mean model conditional on both x and ps fitted by BART
  
  # =========== 1.1) E(Y=1|X,Z) logistic regression
  x <- model.matrix(direct_formula,data)[,-1]
  y <- as.matrix(data[,'Y'])
  
  library(glmnet)
  cv.lasso.direct <- cv.glmnet(x,y,alpha=1,family='gaussian',type.measure="mse")
  
  data_lr_1 <- data_lr_0 <- data
  data_lr_1$Z <- 1
  data_lr_0$Z <- 0
  data_lr_0 <- model.matrix(direct_formula,data_lr_0)[,-1]
  data_lr_1 <- model.matrix(direct_formula,data_lr_1)[,-1]
  outcome_lr_0 <- predict(cv.lasso.direct,data_lr_0,s='lambda.min',type = 'response')
  outcome_lr_1 <- predict(cv.lasso.direct,data_lr_1,s='lambda.min',type='response')
  data$direct_lr_0 <- outcome_lr_0
  data$direct_lr_1 <- outcome_lr_1
  
  # =========== 1.2) E(Y=1|X,Z) BART
  length_test <- dim(data)[1]
  train <- data[,c('Z',xvars)]
  y <- data[,'Y']
  test <- data[,c('Z',xvars)]
  test_bart_1 <- test_bart_0 <- test
  test_bart_1$Z <- 1
  test_bart_0$Z <- 0
  test <- rbind(test_bart_1,test_bart_0)
  post_mean <- wbart(train,y,test, ndpost = ndpost)
  post_mean$prob.test <- post_mean$yhat.test
  post_mean$prob.test.mean <- apply(post_mean$prob.test,2,mean)
  bart_mean <- post_mean$prob.test.mean
  data$direct_bart_1 <- bart_mean[1:length_test]
  data$direct_bart_0 <- bart_mean[(length_test+1):(2*length_test)]
  
  # =========== 2.1) E(Y=1|X,Z=1),E(Y=1|X,Z=0) logistic regression
  covariates.all <- model.matrix(dr_formula,data)[,-1]
  test.all <- model.matrix(dr_formula,data)[,-1]
  
  data1 <- data[data[,'Z']==1,]
  covariates <- model.matrix(dr_formula,data1)[,-1]
  y <- as.matrix(data1[,'Y'])
  cv.lasso.m1 <- cv.glmnet(covariates,y,alpha=1,family='gaussian',type.measure = "mse")
  
  data0 <- data[data[,'Z']==0,]
  covariates <- model.matrix(dr_formula,data0)[,-1]
  y <- as.matrix(data0[,'Y'])
  cv.lasso.m0 <- cv.glmnet(covariates,y,alpha=1,family='gaussian',type.measure="mse")
  
  data$m1.lr.hat <- predict(cv.lasso.m1,test.all,s='lambda.min',type='response')
  data$m0.lr.hat <- predict(cv.lasso.m0,test.all,s='lambda.min',type='response')
  
  # =========== 2.2) E(Y=1|X,Z=1),E(Y=1|X,Z=0) BART
  m.bart.hat <- bart.m.hat(data,xvars,ntree=50,ndpost=ndpost)
  data$m1.bart.hat <- m.bart.hat[,1]
  data$m0.bart.hat <- m.bart.hat[,2]
  
  
  # =========== 3) E(Y=1|X,Z,PS(BART))
  length_test <- dim(data)[1]
  train <- data[,c('Z',xvars,'bartps')]
  #train$bartps <- ifelse(train$bartps<0.1,0.1,ifelse(train$bartps>0.9,0.9,train$bartps))
  y <- data[,'Y']
  test <- data[,c('Z',xvars,'bartps')]
  #test$bartps <- ifelse(test$bartps<0.1,0.1,ifelse(test$bartps>0.9,0.9,test$bartps))
  test_bart_1 <- test_bart_0 <- test
  test_bart_1$Z <- 1
  test_bart_0$Z <- 0
  test <- rbind(test_bart_1,test_bart_0)
  post_mean <- wbart(train,y,test,ndpost=ndpost)
  post_mean$prob.test <- post_mean$yhat.test
  post_mean$prob.test.mean <- apply(post_mean$prob.test,2,mean)
  bart_mean <- post_mean$prob.test.mean
  data$direct_ps_bart_1 <- bart_mean[1:length_test]
  data$direct_ps_bart_0 <- bart_mean[(length_test+1):(2*length_test)]
  
  return(data)
}


ate_ipw <- function(data,ps){
  z <- data$Z
  pi <- data[,ps]
  y <- data$Y*data$weight
  p1 <- 1/sum((z/pi)*data$weight)*sum(z*y/pi)
  p0 <- 1/sum((1-z)/(1-pi)*data$weight)*sum((1-z)*y/(1-pi))
  ate <- p1-p0
  return(ate)
}

ate_dr <- function(data,m1,m0,ps){
  Z <- data$Z
  Y <- data$Y*data$weight
  P <- data[,ps]
  IPTW1 <- Z*Y/P
  Augmentation1 <- (1-Z/P)*m1*data$weight
  IPTW0 <- (1-Z)*Y/(1-P)
  Augmentation0 <- (1-(1-Z)/(1-P))*m0*data$weight
  DR1 <- IPTW1+Augmentation1
  DR0 <- IPTW0+Augmentation0
  ate <- sum(DR1)/sum(data$weight)-sum(DR0)/sum(data$weight)
  return(ate)
}

ResultsEstimators <- function(data){
  naive <- sum(data[data$Z==1,'Y']*data[data$Z==1,'weight'])/sum(data[data$Z==1,'weight'])-sum(data[data$Z==0,'Y']*data[data$Z==0,'weight'])/sum(data[data$Z==0,'weight'])
  ipw_lr <- ate_ipw(data,'lrps')
  ipw_bart <- ate_ipw(data,'bartps')
  dr_lr_lr <- ate_dr(data,data$m1.lr.hat,data$m0.lr.hat,'lrps')
  dr_lr_bart <- ate_dr(data,data$m1.lr.hat,data$m0.lr.hat,'bartps')
  dr_bart_lr <- ate_dr(data,data$m1.bart.hat,data$m0.bart.hat,'lrps')
  dr_bart_bart <- ate_dr(data,data$m1.bart.hat,data$m0.bart.hat,'bartps')
  dr_lr <- sum(data$direct_lr_1*data$weight)/sum(data$weight)-sum(data$direct_lr_0*data$weight)/sum(data$weight)
  dr_bart <- sum(data$direct_bart_1*data$weight)/sum(data$weight)-sum(data$direct_bart_0*data$weight)/sum(data$weight)
  dr_psbart <- sum(data$direct_ps_bart_1*data$weight)/sum(data$weight)-sum(data$direct_ps_bart_0*data$weight)/sum(data$weight)
  return(c(naive,dr_lr,dr_bart,dr_psbart,ipw_lr,ipw_bart,dr_lr_lr,dr_lr_bart,dr_bart_lr,dr_bart_bart))
}

Mu.DGP <- function(xA,xB,xC,linear){
  if(linear==TRUE){
    Mu <- 0.5*xA+2*xB-xC
    #Mu <- xA+xB+xC
  } else{
    Mu <- log(xA^2*sqrt(xB)+xC+0.5)
  }
  return(Mu)
}

tau.DGP <- function(xA,xB,xC,hetero,taulinear=NULL,weight.same=NULL){
  if(hetero==TRUE){
    if(taulinear==TRUE){
      if(weight.same==TRUE){
        tau <- xA + xB + xC 
      } else{
        tau <- 3 * xA + 2*xB + xC
      }
    } else{
      tau <- xA*xB*xC
    }
  } else{
    tau <- rep(3,length(xA))
  }
  return(tau)
}


DGP <- function(Rho,n.pop,p.a.succ,p.b.succ,p.c.succ,Linear,Hetero,TauLinear,Weight.same){
  mu <- rep(0, 3)
  Sigma <- matrix(c(1,Rho,Rho,
                    Rho,1,Rho,
                    Rho,Rho,1), nrow = 3, ncol = 3)
  rawvars <- mvrnorm(n = n.pop, mu = mu, Sigma = Sigma)
  pvars <- pnorm(rawvars) #cdf when x=rawvars, range in (0,1)
  
  xA <- 1*(pvars[,1] < p.a.succ)
  xB <- 1*(pvars[,2] < p.b.succ)
  xC <- 1*(pvars[,3] < p.c.succ)
  
  Mu <- Mu.DGP(xA,xB,xC,linear=Linear)
  tau <- tau.DGP(xA,xB,xC,hetero=Hetero,taulinear=TauLinear,weight.same=Weight.same)
  population <- as.data.frame(cbind(xA,xB,xC,Mu,tau))
  return(population)
}


ProbUnitSelection <- function(population,unitselect,pR=NULL){
  if(unitselect=='xA'){
    prob <- ifelse(population$xA==1,pR,1-pR)
  } else if(unitselect=='xB'){
    prob <- ifelse(population$xB==1,pR,1-pR)
  } else if(unitselect=='xC'){
    prob <- ifelse(population$xC==1,pR,1-pR)
  } else{
    sum.x <- population$xA + population$xB + population$xC
    prob <- inv.logit((log(10)*(sum.x-mean(sum.x)))/sd(sum.x))
  }
  return(prob)
}

GenerateStrataResults <- function(dat.rct,dat.obs){
  rct.xa1 <- mean(dat.rct[(dat.rct$xA==1)&(dat.rct$Z==1),'Y']) - mean(dat.rct[(dat.rct$xA==1)&(dat.rct$Z==0),'Y'])
  rct.xa0 <- mean(dat.rct[(dat.rct$xA==0)&(dat.rct$Z==1),'Y']) - mean(dat.rct[(dat.rct$xA==0)&(dat.rct$Z==0),'Y'])
  rct.xb1 <- mean(dat.rct[(dat.rct$xB==1)&(dat.rct$Z==1),'Y']) - mean(dat.rct[(dat.rct$xB==1)&(dat.rct$Z==0),'Y'])
  rct.xb0 <- mean(dat.rct[(dat.rct$xB==0)&(dat.rct$Z==1),'Y']) - mean(dat.rct[(dat.rct$xB==0)&(dat.rct$Z==0),'Y'])
  rct.xc1 <- mean(dat.rct[(dat.rct$xC==1)&(dat.rct$Z==1),'Y']) - mean(dat.rct[(dat.rct$xC==1)&(dat.rct$Z==0),'Y'])
  rct.xc0 <- mean(dat.rct[(dat.rct$xC==0)&(dat.rct$Z==1),'Y']) - mean(dat.rct[(dat.rct$xC==0)&(dat.rct$Z==0),'Y'])
  
  dat.stra <- dat.obs[dat.obs$xA==1,] ; obs.xa1 <- c(rct.xa1,mean(dat.stra$tau),ResultsEstimators(dat.stra))
  dat.stra <- dat.obs[dat.obs$xA==0,] ; obs.xa0 <- c(rct.xa0,mean(dat.stra$tau),ResultsEstimators(dat.stra))
  dat.stra <- dat.obs[dat.obs$xB==1,] ; obs.xb1 <- c(rct.xb1,mean(dat.stra$tau),ResultsEstimators(dat.stra))
  dat.stra <- dat.obs[dat.obs$xB==0,] ; obs.xb0 <- c(rct.xb0,mean(dat.stra$tau),ResultsEstimators(dat.stra))
  dat.stra <- dat.obs[dat.obs$xC==1,] ; obs.xc1 <- c(rct.xc1,mean(dat.stra$tau),ResultsEstimators(dat.stra))
  dat.stra <- dat.obs[dat.obs$xC==0,] ; obs.xc0 <- c(rct.xc0,mean(dat.stra$tau),ResultsEstimators(dat.stra))
  
  obs.stra <- rbind(obs.xa1,obs.xa0,
                    obs.xb1,obs.xb0,
                    obs.xc1,obs.xc0)
  
  return(obs.stra)
  
}


# distance between weighted obs and true treatment effect of a RCT
GenerateOverallResults <- function(dat.rct,dat.obs){
  outcome <- c(mean(dat.rct$tau),
               mean(dat.obs$tau),
               ResultsEstimators(dat.obs))
  return(outcome)
}


# output: n x m data.frame, 
# n is the number of rows (number of possible pairwise correlations)
# m is the number of columns (tau of rct, tau of obs, estimated tau of a naive method and 9 methods, namely outcome regression, IPW, and DR )
GenerateOverallResultsEstimateRCT <- function(dat.rct,dat.obs){
  outcome <- c(mean(dat.rct[dat.rct$Z==1,'Y']) - mean(dat.rct[dat.rct$Z==0,'Y']),
               mean(dat.obs$tau),
               ResultsEstimators(dat.obs))
  return(outcome)
}


GenerateOverallResultsTrueRCT <- function(dat.rct,dat.obs){
  outcome <- c(mean(dat.rct$tau),
               mean(dat.obs$tau),
               ResultsEstimators(dat.obs))
  return(outcome)
}


library(dplyr)

PrintDistance <- function(dat,nsubvar){
  if((nsubvar==0)|(nsubvar==1)){
    dat.mean <- as.data.frame(dat)
    dat4measure <- dat.mean %>%
      mutate_if(is.numeric, funs(c(first(.), (. - first(.))[-1])) )
    dat4measure <- rowMeans(abs(dat4measure[2:12,]))
  } else {
    dat.mean <- as.data.frame(apply(dat, c(2,3), mean))
    dat4measure <- dat.mean %>%
      mutate_if(is.numeric, funs(c(first(.), (. - first(.))[-1])) )
    dat4measure <- rowMeans(abs(dat4measure[2:12,]))
  }
  return(dat4measure)
}


PrintResults4OverallStratified <- function(Outcome,nsubvar){
  #browser()
  overall <- round(PrintDistance(Outcome,nsubvar),digits = 3)
  #colnames(overall) <- c('True','Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
  #                      'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  return(overall[2:11])
  # Tauxa1 <- round(PrintDistance(Stratified.xa1,nsubvar),digits = 3)
  # Tauxa0 <- round(PrintDistance(Stratified.xa0,nsubvar),digits = 3)
  # Tauxb1 <- round(PrintDistance(Stratified.xb1,nsubvar),digits = 3)
  # Tauxb0 <- round(PrintDistance(Stratified.xb0,nsubvar),digits = 3)
  # Tauxc1 <- round(PrintDistance(Stratified.xc1,nsubvar),digits = 3)
  # Tauxc0 <- round(PrintDistance(Stratified.xc0,nsubvar),digits = 3)
  # Results4Print <- rbind('&',overall,'&',Tauxa1,'&',Tauxa0,'&',Tauxb1,'&',Tauxb0,'&',Tauxc1,'&',Tauxc0,'\\')
  # Results4Print <- t(as.matrix(Results4Print))
  # Results4Print <- as.data.frame(Results4Print)
  # methods.name <- c('multirow{11}{*}{$X_{A}$}&True','&Crude','&Rg/LReg','&Rg/BART','&Rg/ps-BART','&IPW/LReg','&IPW/BART',
  #                   '&DR/LReg-LReg','&DR/LReg-BART','&DR/BART-LReg','&DR/BART-BART')
  # Results4Print <- cbind(methods.name,Results4Print)
  # 
  # colnames(Results4Print) <- c('methods','','d(sub;rct)','','d(sub;rct;xa=1)','','d(sub;rct;xa=0)','','d(sub;rct;xb=1)','','d(sub;rct;xb=0)','','d(sub;rct;xc=1)','','d(sub;rct;xc=0)','')
  #return(Results4Print)
}


PrintTable <- function(dat.rct,dat.obs){
  n.rct <- dim(dat.rct)[1]; n.obs <- dim(dat.obs)[1]
  dif.xa <- abs(table(dat.rct$xA)[1]/n.rct - table(dat.obs$xA)[1]/n.obs)
  dif.xb <- abs(table(dat.rct$xB)[1]/n.rct - table(dat.obs$xB)[1]/n.obs)
  dif.xc <- abs(table(dat.rct$xC)[1]/n.rct - table(dat.obs$xC)[1]/n.obs)
  return(dif.xa+dif.xb+dif.xc)
}

DKL <- function(data.rct,data.obs){
  p.matrix <- matrix(c(mean((data.rct$xA==0)&(data.rct$xB==0)&(data.rct$xC==0)),
                       mean((data.obs$xA==0)&(data.obs$xB==0)&(data.obs$xC==0)),
                       mean((data.rct$xA==0)&(data.rct$xB==0)&(data.rct$xC==1)),
                       mean((data.obs$xA==0)&(data.obs$xB==0)&(data.obs$xC==1)),
                       mean((data.rct$xA==0)&(data.rct$xB==1)&(data.rct$xC==0)),
                       mean((data.obs$xA==0)&(data.obs$xB==1)&(data.obs$xC==0)),
                       mean((data.rct$xA==0)&(data.rct$xB==1)&(data.rct$xC==1)),
                       mean((data.obs$xA==0)&(data.obs$xB==1)&(data.obs$xC==1)),
                       mean((data.rct$xA==1)&(data.rct$xB==0)&(data.rct$xC==0)),
                       mean((data.obs$xA==1)&(data.obs$xB==0)&(data.obs$xC==0)),
                       mean((data.rct$xA==1)&(data.rct$xB==0)&(data.rct$xC==1)),
                       mean((data.obs$xA==1)&(data.obs$xB==0)&(data.obs$xC==1)),
                       mean((data.rct$xA==1)&(data.rct$xB==1)&(data.rct$xC==0)),
                       mean((data.obs$xA==1)&(data.obs$xB==1)&(data.obs$xC==0)),
                       mean((data.rct$xA==1)&(data.rct$xB==1)&(data.rct$xC==1)),
                       mean((data.obs$xA==1)&(data.obs$xB==1)&(data.obs$xC==1))),
                     ncol = 2, byrow = TRUE
  )
  DKL <- sum(apply(p.matrix, 1, function(x) x[1]*log(x[1]/x[2])))
  return(DKL)
}


GetWeight <- function(weight,domain,value){
  if((length(domain)==1)&&(domain=='xA')){
    p.weight <- weight[weight$xA==value[1],'weight']
  } else if ((length(domain)==1)&&(domain=='xB')){
    p.weight <- weight[weight$xB==value[2],'weight']
  } else if ((length(domain)==1)&&(domain=='xC')){
    p.weight <- weight[weight$xC==value[3],'weight']
  } else if((length(domain)==2)&&(sum(domain==c('xA','xB'))==2)){
    p.weight <- weight[(weight$xA==value[1])&(weight$xB==value[2]),'weight']
  } else if ((length(domain)==2)&&(sum(domain==c('xA','xC'))==2)){
    p.weight <- weight[(weight$xA==value[1])&(weight$xC==value[3]),'weight']
  } else if ((length(domain)==2)&&(sum(domain==c('xB','xC'))==2)){
    p.weight <- weight[(weight$xB==value[2])&(weight$xC==value[3]),'weight']
  } else{
    p.weight <- weight[(weight$xA==value[1])&(weight$xB==value[2])&(weight$xC==value[3]),'weight']
  }
  
  return(p.weight)
}

#data.rct <- dat.rct
#data.obs <- dat.obs
#domain <- c('xA','xB')
DKLWeighting <- function(data.rct,data.obs,weight,domain){
  p.matrix <- matrix(c(mean((data.rct$xA==0)&(data.rct$xB==0)&(data.rct$xC==0)),
                       mean((data.obs$xA==0)&(data.obs$xB==0)&(data.obs$xC==0)),
                       GetWeight(weight,domain,value=c(0,0,0)),
                       mean((data.rct$xA==0)&(data.rct$xB==0)&(data.rct$xC==1)),
                       mean((data.obs$xA==0)&(data.obs$xB==0)&(data.obs$xC==1)),
                       GetWeight(weight,domain,value=c(0,0,1)),
                       mean((data.rct$xA==0)&(data.rct$xB==1)&(data.rct$xC==0)),
                       mean((data.obs$xA==0)&(data.obs$xB==1)&(data.obs$xC==0)),
                       GetWeight(weight,domain,value=c(0,1,0)),
                       mean((data.rct$xA==0)&(data.rct$xB==1)&(data.rct$xC==1)),
                       mean((data.obs$xA==0)&(data.obs$xB==1)&(data.obs$xC==1)),
                       GetWeight(weight,domain,value=c(0,1,1)),
                       mean((data.rct$xA==1)&(data.rct$xB==0)&(data.rct$xC==0)),
                       mean((data.obs$xA==1)&(data.obs$xB==0)&(data.obs$xC==0)),
                       GetWeight(weight,domain,value=c(1,0,0)),
                       mean((data.rct$xA==1)&(data.rct$xB==0)&(data.rct$xC==1)),
                       mean((data.obs$xA==1)&(data.obs$xB==0)&(data.obs$xC==1)),
                       GetWeight(weight,domain,value=c(1,0,1)),
                       mean((data.rct$xA==1)&(data.rct$xB==1)&(data.rct$xC==0)),
                       mean((data.obs$xA==1)&(data.obs$xB==1)&(data.obs$xC==0)),
                       GetWeight(weight,domain,value=c(1,1,0)),
                       mean((data.rct$xA==1)&(data.rct$xB==1)&(data.rct$xC==1)),
                       mean((data.obs$xA==1)&(data.obs$xB==1)&(data.obs$xC==1)),
                       GetWeight(weight,domain,value=c(1,1,1))),
                       ncol = 3, byrow = TRUE
  )
  
  #DKL <- sum(apply(p.matrix, 1, function(x) (x[2]*x[3])*log((x[2]*x[3])/(x[1]))))
  DKL <- sum(apply(p.matrix, 1, function(x) x[1]*log(x[1]/(x[2]*x[3]))))
  return(DKL)
}


WeightCalculation <- function(source,target,domain,p.rct=NULL){
  if(is.null(p.rct)){
    n.domain <- length(domain)
    weight <- unique(source[,domain])
    weight <- as.data.frame(weight)
    colnames(weight) <- domain
    if(n.domain==1){
      for (i in seq(dim(weight)[1])){
        p.target <- mean(target[,domain]==weight[i,domain])
        p.source <- mean(source[,domain]==weight[i,domain])
        weight[i,'weight'] <- p.target/p.source
      }
    } else{
      for (i in seq(dim(weight)[1])){
        p.target <- rowSums(target[,domain]==c(weight[i,domain]))
        p.target <- mean(p.target==length(domain))
        p.source <- rowSums(source[,domain]==c(weight[i,domain]))
        p.source <- mean(p.source==length(domain))
        weight[i,'weight'] <- p.target/p.source
      }
    }
  } else{
    n.domain <- length(domain)
    weight <- unique(source[,domain])
    weight <- as.data.frame(weight)
    colnames(weight) <- domain
    if(n.domain==1){
      for (i in seq(dim(weight)[1])){
        index.rct <- (p.rct[,n.domain]== weight[i,n.domain])
        p.target <- p.rct[index.rct,'prob']
        p.source <- mean(source[,domain]==weight[i,domain])
        weight[i,'weight'] <- p.target/p.source
      }
    } else{
      for (i in seq(dim(weight)[1])){
        index.rct <- TRUE
        for (j in seq(n.domain)){
          aa <- (p.rct[,j]==weight[i,j])
          index.rct <- index.rct&aa
        }
        p.target <- p.rct[index.rct,'prob']
        p.source <- rowSums(source[,domain]==c(weight[i,domain]))
        p.source <- mean(p.source==length(domain))
        weight[i,'weight'] <- p.target/p.source
      }
    }
  }
  return(weight)
}

Weight_normalize <- function(dat.obs){
  s000 <- filter(dat.obs,(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==0))
  s001 <- filter(dat.obs,(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==1))
  s010 <- filter(dat.obs,(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==0))
  s011 <- filter(dat.obs,(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==1))
  s100 <- filter(dat.obs,(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==0))
  s101 <- filter(dat.obs,(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==1))
  s110 <- filter(dat.obs,(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==0))
  s111 <- filter(dat.obs,(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==1))
  
  N <- nrow(dat.obs)
  n000 <- nrow(s000); w000 <- s000$weight[1]; p000 <- n000/N
  n001 <- nrow(s001); w001 <- s001$weight[1]; p001 <- n001/N
  n010 <- nrow(s010); w010 <- s010$weight[1]; p010 <- n010/N
  n011 <- nrow(s011); w011 <- s011$weight[1]; p011 <- n011/N
  n100 <- nrow(s100); w100 <- s100$weight[1]; p100 <- n100/N
  n101 <- nrow(s101); w101 <- s101$weight[1]; p101 <- n101/N
  n110 <- nrow(s110); w110 <- s110$weight[1]; p110 <- n110/N
  n111 <- nrow(s111); w111 <- s111$weight[1]; p111 <- n111/N
  
  sumP <- p000*w000+p001*w001+p010*w010+p010*w010+p011*w011+
    p100*w100+p101*w101+p110*w110+p110*w110+p111*w111
  
  #dat.obs$weight <- dat.obs$weight/sumP
  return(sumP)
}



ComputeWeightedEstimatesBiasVariance <- function(Rho,p.a.succ,p.b.succ,p.c.succ,
                                                 n.pop,n.rct,n.obs,N,
                                                 lambda,
                                                 Linear,Hetero,
                                                 TauLinear,Weight.same=NULL,
                                                 unitselect,pR,
                                                 ntree,ndpost){
  #browser()
  Outcome.0 <-array(NaN,c(13,N))
  Outcome.xa <- array(NaN,c(13,N))
  Outcome.xb <- array(NaN,c(13,N))
  Outcome.xc <- array(NaN,c(13,N))
  Outcome.xab <- array(NaN,c(length(seq(-1,1,0.1)),13,N))
  Outcome.xac <- array(NaN,c(length(seq(-1,1,0.1)),13,N))
  Outcome.xbc <- array(NaN,c(length(seq(-1,1,0.1)),13,N))
  Outcome.xabc <- array(NaN,c(16,13,N))
  Outcome.xabc.Rho <- array(NaN,c(13,N))

  n <- 0
  while (n < N) {
    n <- n+1
    
    # step 1: DGP
    population <- DGP(Rho,n.pop,p.a.succ,p.b.succ,p.c.succ,Linear,Hetero,TauLinear,Weight.same)
    tau.p <- mean(population$tau)
    
    # step 2: sampling R=1/0
    dat.rct <- population[sample(1:n.pop,size = n.rct),]
    prob <- ProbUnitSelection(population,unitselect=unitselect,pR=pR)
    dat.obs <- population[sample(1:dim(population)[1],size = n.obs,prob = prob),]

    # step 3: treatment assignment Z=1/0; non-randomization induced confounding
    pi <- inv.logit((lambda*(dat.obs$Mu-mean(dat.obs$Mu)))/sd(dat.obs$Mu))
    dat.obs$Z <- rbinom(n.obs,1,pi)
    dat.obs$Y <- dat.obs$Mu+dat.obs$Z*dat.obs$tau+rnorm(n.obs)
    dat.rct$Z <- rbinom(n.rct,1,0.5)
    dat.rct$Y <- dat.rct$Mu+dat.rct$Z*dat.rct$tau+rnorm(n.rct)
    
    # the estimate is average treatment effect; return conterfactual computation for each unit
    xvars <- c('xA','xB','xC')
    dat.obs <- potential_outcomes_direct_ipw_dr(dat.obs,xvars,ntree,ndpost)
   
    pa <- (table(dat.rct[,1])/dim(dat.rct)[1])[2] #success probablity of xA in rct
    pb <- (table(dat.rct[,2])/dim(dat.rct)[1])[2] #success probablity of xB in rct
    pc <- (table(dat.rct[,3])/dim(dat.rct)[1])[2] #success probablity of xC in rct
    
    # ================== Methods comparsion: ate obtained from each methods
    dat.obs$weight <- 1

    # ================== 
    # step 1: no subsampling 
    dat.obs$weight <- 1
    Outcome.0[,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))
    # ==================
    # subsampling based on one confounding, xA
    domain <- c('xA')
    weight <- WeightCalculation(dat.obs,dat.rct,domain)
    dat.obs$weight <- weight[match(dat.obs$xA,weight$xA),'weight']
    Outcome.xa[,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))
    # ==================
    # subsampling based on one confounding, xB
    weight <- WeightCalculation(dat.obs,dat.rct,'xB')
    dat.obs$weight <- weight[match(dat.obs$xB,weight$xB),'weight']
    Outcome.xb[,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))
    # ==================
    # subsampling based on one confounding, xC
    weight <- WeightCalculation(dat.obs,dat.rct,'xC')
    dat.obs$weight <- weight[match(dat.obs$xC,weight$xC),'weight']
    Outcome.xc[,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))
    # ==================
    # subsampling based on two confounding: xA,xB
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      p00 <- cbind(0,0,pmvnorm(lower = c(-Inf,-Inf),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1)),mean=mu,corr = Sigma))
      p01 <- cbind(0,1,pmvnorm(lower = c(-Inf,qnorm(pb,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),Inf),mean=mu,corr = Sigma))
      p10 <- cbind(1,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf),upper = c(Inf,qnorm(pb,mean=0,sd=1)),mean=mu,corr = Sigma))
      p11 <- cbind(1,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1)),upper = c(Inf,Inf),mean=mu,corr = Sigma))
      rct.prob <- rbind(p00,p01,p10,p11)
      colnames(rct.prob) <- c('xA','xB','prob')
      weight <- WeightCalculation(dat.obs,dat.rct,c('xA','xB'),p.rct=rct.prob)
  
      # specify the weight for each unit in dat.obs
      domain <- c('xA','xB')
      rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
      dat.obs$weight <- weight[rowID.weight,'weight']
      Outcome.xab[i,,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))

      i <- i+1
    }

    # ==================
    # subsampling based on two confounding: xA,xC
    mu <- rep(0, 2)
    i <- 1
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      p00 <- cbind(0,0,pmvnorm(lower = c(-Inf,-Inf),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p01 <- cbind(0,1,pmvnorm(lower = c(-Inf,qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),Inf),mean=mu,corr = Sigma))
      p10 <- cbind(1,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf),upper = c(Inf,qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p11 <- cbind(1,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(Inf,Inf),mean=mu,corr = Sigma))
      rct.prob <- rbind(p00,p01,p10,p11)
      colnames(rct.prob) <- c('xA','xC','prob')
      weight <- WeightCalculation(dat.obs,dat.rct,c('xA','xC'),p.rct=rct.prob)
      
      # specify the weight for each unit in dat.obs
      domain <- c('xA','xC')
      rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
      dat.obs$weight <- weight[rowID.weight,'weight']
      Outcome.xac[i,,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))

      i <- i+1
    }

    # ==================
    # subsampling based on two confounding: xB,xC
    mu <- rep(0, 2)
    i <- 1
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      p00 <- cbind(0,0,pmvnorm(lower = c(-Inf,-Inf),upper = c(qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p01 <- cbind(0,1,pmvnorm(lower = c(-Inf,qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pb,mean=0,sd=1),Inf),mean=mu,corr = Sigma))
      p10 <- cbind(1,0,pmvnorm(lower = c(qnorm(pb,mean=0,sd=1),-Inf),upper = c(Inf,qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p11 <- cbind(1,1,pmvnorm(lower = c(qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(Inf,Inf),mean=mu,corr = Sigma))
      rct.prob <- rbind(p00,p01,p10,p11)
      colnames(rct.prob) <- c('xB','xC','prob')
      weight <- WeightCalculation(dat.obs,dat.rct,c('xB','xC'),p.rct=rct.prob)
      domain <- c('xB','xC')
      rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
      dat.obs$weight <- weight[rowID.weight,'weight']
      Outcome.xbc[i,,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))
      i <- i+1
    }

    # ==================
    # subsampling based on three confounding: xA,xB,xC
    mu <- rep(0, 3)
    i <- 1
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,rho,
                       rho, 1,rho,
                       rho,rho,1), nrow = 3, ncol = 3)
      if(is.positive.semi.definite(Sigma)){
        p000 <- cbind(0,0,0,pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,sigma = Sigma))
        p001 <- cbind(0,0,1,pmvnorm(lower = c(-Inf,-Inf,qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),Inf),mean=mu,sigma= Sigma))
        p010 <- cbind(0,1,0,pmvnorm(lower = c(-Inf,qnorm(pb,mean=0,sd=1),-Inf),upper = c(qnorm(pa,mean=0,sd=1),Inf,qnorm(pc,mean=0,sd=1)),mean=mu,sigma= Sigma))
        p011 <- cbind(0,1,1,pmvnorm(lower = c(-Inf,qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),Inf,Inf),mean=mu,sigma= Sigma))
        p100 <- cbind(1,0,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf,-Inf),upper = c(Inf,qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,sigma= Sigma))
        p101 <- cbind(1,0,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf,qnorm(pc,mean=0,sd=1)),upper = c(Inf,qnorm(pb,mean=0,sd=1),Inf),mean=mu,sigma= Sigma))
        p110 <- cbind(1,1,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),-Inf),upper = c(Inf,Inf,qnorm(pc,mean=0,sd=1)),mean=mu,sigma= Sigma))
        p111 <- cbind(1,1,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(Inf,Inf,Inf),mean=mu,sigma= Sigma))
                
        rct.prob <- rbind(p000,p001,p010,p011,p100,p101,p110,p111)
        colnames(rct.prob) <- c('xA','xB','xC','prob')
        weight <- WeightCalculation(dat.obs,dat.rct,c('xA','xB','xC'),p.rct=rct.prob)
                
        # specify the weight for each unit in dat.obs
        domain <- c('xA','xB','xC')
        rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
        dat.obs$weight <- weight[rowID.weight,'weight']
        Outcome.xabc[i,,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))

        i <- i+1
        #print(rho)
      }
    }
    
    # ==================
    # subsampling based on three confounding and correct dependency Rho: xA,xB,xC
    domain <- c('xA','xB','xC')
    weight <- WeightCalculation(dat.obs,dat.rct,domain)
    rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
    dat.obs$weight <- weight[rowID.weight,'weight']
    Outcome.xabc.Rho[,n] <- c(tau.p,GenerateOverallResultsEstimateRCT(dat.rct,dat.obs))
  }
  
  return(list(Outcome.0,
              Outcome.xa,
              Outcome.xb,
              Outcome.xc,
              Outcome.xab,
              Outcome.xac,
              Outcome.xbc,
              Outcome.xabc,
              Outcome.xabc.Rho))
}


ComputeDistance2RctAndMethodsBiasWeighting <- function(Rho,p.a.succ,p.b.succ,p.c.succ,
                                                       n.pop,n.rct,n.obs,N,
                                                       lambda,
                                                       Linear,Hetero,
                                                       TauLinear,Weight.same=NULL,
                                                       unitselect,pR,
                                                       ntree,ndpost){
  #browser()
  Outcome.0 <-array(NaN,c(12,N)); Outcome.rcttrue.0 <-array(NaN,c(12,N))
  Outcome.xa <- array(NaN,c(12,N)); Outcome.rcttrue.xa <-array(NaN,c(12,N))
  Outcome.xb <- array(NaN,c(12,N)); Outcome.rcttrue.xb <-array(NaN,c(12,N))
  Outcome.xc <- array(NaN,c(12,N)); Outcome.rcttrue.xc <-array(NaN,c(12,N))
  Outcome.xab <- array(NaN,c(length(seq(-1,1,0.1)),12,N)); Outcome.rcttrue.xab <-array(NaN,c(length(seq(-1,1,0.1)),12,N))
  Outcome.xac <- array(NaN,c(length(seq(-1,1,0.1)),12,N)); Outcome.rcttrue.xac <-array(NaN,c(length(seq(-1,1,0.1)),12,N))
  Outcome.xbc <- array(NaN,c(length(seq(-1,1,0.1)),12,N)); Outcome.rcttrue.xbc <-array(NaN,c(length(seq(-1,1,0.1)),12,N))
  Outcome.xabc <- array(NaN,c(16,12,N)); Outcome.rcttrue.xabc <- array(NaN,c(16,12,N))
  Outcome.xabc.Rho <- array(NaN,c(12,N)); Outcome.rcttrue.xabc.Rho <-array(NaN,c(12,N))
  
  Methods.Est <- array(NaN,c(N,11))
  Diff.dis <- data.frame(x0=numeric(),xa=numeric(),xb=numeric(),xc=numeric(),xab=numeric(),xac=numeric(),xbc=numeric(),xabc=numeric(),xabc.rho=numeric())
  n <- 0
  while (n < N) {
    n <- n+1
    
    # step 1: DGP
    population <- DGP(Rho,n.pop,p.a.succ,p.b.succ,p.c.succ,Linear,Hetero,TauLinear,Weight.same)
    
    # step 2: sampling R=1/0
    dat.rct <- population[sample(1:n.pop,size = n.rct),]
    prob <- ProbUnitSelection(population,unitselect=unitselect,pR=pR)
    dat.obs <- population[sample(1:dim(population)[1],size = n.obs,prob = prob),]
    
    # step 3: treatment assignment Z=1/0; non-randomization induced confounding
    pi <- inv.logit((lambda*(dat.obs$Mu-mean(dat.obs$Mu)))/sd(dat.obs$Mu))
    dat.obs$Z <- rbinom(n.obs,1,pi)
    dat.obs$Y <- dat.obs$Mu+dat.obs$Z*dat.obs$tau+rnorm(n.obs)
    dat.rct$Z <- rbinom(n.rct,1,0.5)
    dat.rct$Y <- dat.rct$Mu+dat.rct$Z*dat.rct$tau+rnorm(n.rct)
    
    # the estimate is average treatment effect; return conterfactual computation for each unit
    xvars <- c('xA','xB','xC')
    dat.obs <- potential_outcomes_direct_ipw_dr(dat.obs,xvars,ntree,ndpost)
    
    pa <- (table(dat.rct[,1])/dim(dat.rct)[1])[2] #success probablity of xA in rct
    pb <- (table(dat.rct[,2])/dim(dat.rct)[1])[2] #success probablity of xB in rct
    pc <- (table(dat.rct[,3])/dim(dat.rct)[1])[2] #success probablity of xC in rct
    
    # ================== Methods comparsion: ate obtained from each methods
    dat.obs$weight <- 1; Methods.Est[n,] <- c(mean(dat.obs$tau),ResultsEstimators(dat.obs))
    
    # ================== 
    # step 1: no subsampling 
    dat.obs$weight <- 1
    Outcome.0[,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
    Outcome.rcttrue.0[,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
    Diff.dis[n,'x0'] <- DKL(dat.rct,dat.obs)    
    # ==================
    # subsampling based on one confounding, xA
    domain <- c('xA')
    weight <- WeightCalculation(dat.obs,dat.rct,domain)
    dat.obs$weight <- weight[match(dat.obs$xA,weight$xA),'weight']
    #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
    #weight$weight <- weight$weight/sumP
    Outcome.xa[,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
    Outcome.rcttrue.xa[,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
    Diff.dis[n,'xa'] <- DKLWeighting(dat.rct,dat.obs,weight,domain=domain)
    # ==================
    # subsampling based on one confounding, xB
    weight <- WeightCalculation(dat.obs,dat.rct,'xB')
    dat.obs$weight <- weight[match(dat.obs$xB,weight$xB),'weight']
    #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
    #weight$weight <- weight$weight/sumP    
    Outcome.xb[,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
    Outcome.rcttrue.xb[,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
    Diff.dis[n,'xb'] <- DKLWeighting(dat.rct,dat.obs,weight,c('xB'))
    # ==================
    # subsampling based on one confounding, xC
    weight <- WeightCalculation(dat.obs,dat.rct,'xC')
    dat.obs$weight <- weight[match(dat.obs$xC,weight$xC),'weight']
    #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
    #weight$weight <- weight$weight/sumP
    Outcome.xc[,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
    Outcome.rcttrue.xc[,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
    Diff.dis[n,'xc'] <- DKLWeighting(dat.rct,dat.obs,weight,c('xC'))
    # ==================
    # subsampling based on two confounding: xA,xB
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      p00 <- cbind(0,0,pmvnorm(lower = c(-Inf,-Inf),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1)),mean=mu,corr = Sigma))
      p01 <- cbind(0,1,pmvnorm(lower = c(-Inf,qnorm(pb,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),Inf),mean=mu,corr = Sigma))
      p10 <- cbind(1,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf),upper = c(Inf,qnorm(pb,mean=0,sd=1)),mean=mu,corr = Sigma))
      p11 <- cbind(1,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1)),upper = c(Inf,Inf),mean=mu,corr = Sigma))
      rct.prob <- rbind(p00,p01,p10,p11)
      colnames(rct.prob) <- c('xA','xB','prob')
      weight <- WeightCalculation(dat.obs,dat.rct,c('xA','xB'),p.rct=rct.prob)
      
      # specify the weight for each unit in dat.obs
      domain <- c('xA','xB')
      rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
      dat.obs$weight <- weight[rowID.weight,'weight']
      #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
      #weight$weight <- weight$weight/sumP
      Outcome.xab[i,,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
      Outcome.rcttrue.xab[i,,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
      diff.dis.mean <- c(diff.dis.mean,DKLWeighting(dat.rct,dat.obs,weight,c('xA','xB')))
      
      i <- i+1
    }
    Diff.dis[n,'xab'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    # ==================
    # subsampling based on two confounding: xA,xC
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      p00 <- cbind(0,0,pmvnorm(lower = c(-Inf,-Inf),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p01 <- cbind(0,1,pmvnorm(lower = c(-Inf,qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),Inf),mean=mu,corr = Sigma))
      p10 <- cbind(1,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf),upper = c(Inf,qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p11 <- cbind(1,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(Inf,Inf),mean=mu,corr = Sigma))
      rct.prob <- rbind(p00,p01,p10,p11)
      colnames(rct.prob) <- c('xA','xC','prob')
      weight <- WeightCalculation(dat.obs,dat.rct,c('xA','xC'),p.rct=rct.prob)
      
      # specify the weight for each unit in dat.obs
      domain <- c('xA','xC')
      rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
      dat.obs$weight <- weight[rowID.weight,'weight']
      #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
      #weight$weight <- weight$weight/sumP
      
      Outcome.xac[i,,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
      Outcome.rcttrue.xac[i,,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
      diff.dis.mean <- c(diff.dis.mean,DKLWeighting(dat.rct,dat.obs,weight,c('xA','xC')))
      
      i <- i+1
    }
    Diff.dis[n,'xac'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    # ==================
    # subsampling based on two confounding: xB,xC
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      p00 <- cbind(0,0,pmvnorm(lower = c(-Inf,-Inf),upper = c(qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p01 <- cbind(0,1,pmvnorm(lower = c(-Inf,qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pb,mean=0,sd=1),Inf),mean=mu,corr = Sigma))
      p10 <- cbind(1,0,pmvnorm(lower = c(qnorm(pb,mean=0,sd=1),-Inf),upper = c(Inf,qnorm(pc,mean=0,sd=1)),mean=mu,corr = Sigma))
      p11 <- cbind(1,1,pmvnorm(lower = c(qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(Inf,Inf),mean=mu,corr = Sigma))
      rct.prob <- rbind(p00,p01,p10,p11)
      colnames(rct.prob) <- c('xB','xC','prob')
      weight <- WeightCalculation(dat.obs,dat.rct,c('xB','xC'),p.rct=rct.prob)
      
      # specify the weight for each unit in dat.obs
      domain <- c('xB','xC')
      rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
      dat.obs$weight <- weight[rowID.weight,'weight']
      #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
      #weight$weight <- weight$weight/sumP
      
      Outcome.xbc[i,,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
      Outcome.rcttrue.xbc[i,,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)    
      diff.dis.mean <- c(diff.dis.mean,DKLWeighting(dat.rct,dat.obs,weight,c('xB','xC')))
      
      i <- i+1
    }
    Diff.dis[n,'xbc'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    # ==================
    # subsampling based on three confounding: xA,xB,xC
    mu <- rep(0, 3)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,rho,
                        rho, 1,rho,
                        rho,rho,1), nrow = 3, ncol = 3)
      if(is.positive.semi.definite(Sigma)){
        p000 <- cbind(0,0,0,pmvnorm(lower = c(-Inf,-Inf,-Inf),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,sigma = Sigma))
        p001 <- cbind(0,0,1,pmvnorm(lower = c(-Inf,-Inf,qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),Inf),mean=mu,sigma= Sigma))
        p010 <- cbind(0,1,0,pmvnorm(lower = c(-Inf,qnorm(pb,mean=0,sd=1),-Inf),upper = c(qnorm(pa,mean=0,sd=1),Inf,qnorm(pc,mean=0,sd=1)),mean=mu,sigma= Sigma))
        p011 <- cbind(0,1,1,pmvnorm(lower = c(-Inf,qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(qnorm(pa,mean=0,sd=1),Inf,Inf),mean=mu,sigma= Sigma))
        p100 <- cbind(1,0,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf,-Inf),upper = c(Inf,qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),mean=mu,sigma= Sigma))
        p101 <- cbind(1,0,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),-Inf,qnorm(pc,mean=0,sd=1)),upper = c(Inf,qnorm(pb,mean=0,sd=1),Inf),mean=mu,sigma= Sigma))
        p110 <- cbind(1,1,0,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),-Inf),upper = c(Inf,Inf,qnorm(pc,mean=0,sd=1)),mean=mu,sigma= Sigma))
        p111 <- cbind(1,1,1,pmvnorm(lower = c(qnorm(pa,mean=0,sd=1),qnorm(pb,mean=0,sd=1),qnorm(pc,mean=0,sd=1)),upper = c(Inf,Inf,Inf),mean=mu,sigma= Sigma))
        
        rct.prob <- rbind(p000,p001,p010,p011,p100,p101,p110,p111)
        colnames(rct.prob) <- c('xA','xB','xC','prob')
        weight <- WeightCalculation(dat.obs,dat.rct,c('xA','xB','xC'),p.rct=rct.prob)
        
        # specify the weight for each unit in dat.obs
        domain <- c('xA','xB','xC')
        rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
        dat.obs$weight <- weight[rowID.weight,'weight']
        #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
        #weight$weight <- weight$weight/sumP
        Outcome.xabc[i,,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
        Outcome.rcttrue.xabc[i,,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
        diff.dis.mean <- c(diff.dis.mean,DKLWeighting(dat.rct,dat.obs,weight,c('xA','xB','xC')))
        
        i <- i+1
        #print(rho)
      }
    }
    Diff.dis[n,'xabc'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    
    # ==================
    # subsampling based on three confounding and correct dependency Rho: xA,xB,xC
    diff.dis.mean <- NULL
    domain <- c('xA','xB','xC')
    weight <- WeightCalculation(dat.obs,dat.rct,domain)
    # specify the weight for each unit in dat.obs
    rowID.weight <- row.match(dat.obs[,domain],weight[,domain])
    dat.obs$weight <- weight[rowID.weight,'weight']
    #sumP <- Weight_normalize(dat.obs); dat.obs$weight <- dat.obs$weight/sumP
    #weight$weight <- weight$weight/sumP
    Outcome.xabc.Rho[,n] <- GenerateOverallResultsEstimateRCT(dat.rct,dat.obs)
    Outcome.rcttrue.xabc.Rho[,n] <- GenerateOverallResultsTrueRCT(dat.rct,dat.obs)
    Diff.dis[n,'xabc.rho'] <- DKLWeighting(dat.rct,dat.obs,weight,c('xA','xB','xC'))
  }
  
  results.0 <- PrintResults4OverallStratified(Outcome.0,nsubvar=0)
  results.rcttrue.0 <- PrintResults4OverallStratified(Outcome.rcttrue.0,nsubvar=0)
  
  results.xa  <- PrintResults4OverallStratified(Outcome.xa, nsubvar=1)
  results.rcttrue.xa  <- PrintResults4OverallStratified(Outcome.rcttrue.xa, nsubvar=1)
  
  results.xb  <- PrintResults4OverallStratified(Outcome.xb,nsubvar=1)
  results.rcttrue.xb  <- PrintResults4OverallStratified(Outcome.rcttrue.xb,nsubvar=1)
  
  results.xc  <- PrintResults4OverallStratified(Outcome.xc,nsubvar=1)
  results.rcttrue.xc  <- PrintResults4OverallStratified(Outcome.rcttrue.xc,nsubvar=1)
  
  results.xab <- PrintResults4OverallStratified(Outcome.xab,nsubvar=2)
  results.rcttrue.xab <- PrintResults4OverallStratified(Outcome.rcttrue.xab,nsubvar=2)
  
  results.xac <- PrintResults4OverallStratified(Outcome.xac,nsubvar=2)
  results.rcttrue.xac <- PrintResults4OverallStratified(Outcome.rcttrue.xac,nsubvar=2)
  
  results.xbc <- PrintResults4OverallStratified(Outcome.xbc,nsubvar=2)
  results.rcttrue.xbc <- PrintResults4OverallStratified(Outcome.rcttrue.xbc,nsubvar=2)
  
  results.xabc <- PrintResults4OverallStratified(Outcome.xabc,nsubvar=3)
  results.rcttrue.xabc <- PrintResults4OverallStratified(Outcome.rcttrue.xabc,nsubvar=3)
  
  results.xabc.Rho <- PrintResults4OverallStratified(Outcome.xabc.Rho,nsubvar=1)
  results.rcttrue.xabc.Rho <- PrintResults4OverallStratified(Outcome.rcttrue.xabc.Rho,nsubvar=1)
  
  ResultsEstimate <- cbind(results.0,results.xa,results.xb,results.xc,
                           results.xab,results.xac,results.xbc,
                           results.xabc,results.xabc.Rho)
  rownames(ResultsEstimate) <- c('Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                                 'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  
  
  ResultsTrue <- cbind(results.rcttrue.0,results.rcttrue.xa,results.rcttrue.xb,results.rcttrue.xc,
                       results.rcttrue.xab,results.rcttrue.xac,results.rcttrue.xbc,
                       results.rcttrue.xabc,results.rcttrue.xabc.Rho)
  rownames(ResultsTrue) <- c('Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                             'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  
  return(list(Methods.Est,ResultsEstimate,ResultsTrue,Diff.dis))
}



ComputeDistance2RctAndMethodsBias <- function(Rho,p.a.succ,p.b.succ,p.c.succ,
                                              n.pop,n.rct,n.obs,N,
                                              lambda,
                                              Linear,Hetero,
                                              TauLinear,Weight.same=NULL,
                                              unitselect,pR,
                                              ntree,ndpost){
  #browser()
  Outcome.0 <- Stratified.0.xa1 <- Stratified.0.xa0 <- Stratified.0.xb1 <- Stratified.0.xb0 <- Stratified.0.xc1 <- Stratified.0.xc0 <- array(NaN,c(12,N))
  Outcome.xa <- Stratified.xa.xa1 <- Stratified.xa.xa0 <- Stratified.xa.xb1 <- Stratified.xa.xb0 <- Stratified.xa.xc1 <- Stratified.xa.xc0 <- array(NaN,c(12,N))
  Outcome.xb <- Stratified.xb.xa1 <- Stratified.xb.xa0 <- Stratified.xb.xb1 <- Stratified.xb.xb0 <- Stratified.xb.xc1 <- Stratified.xb.xc0 <- array(NaN,c(12,N))
  Outcome.xc <- Stratified.xc.xa1 <- Stratified.xc.xa0 <- Stratified.xc.xb1 <- Stratified.xc.xb0 <- Stratified.xc.xc1 <- Stratified.xc.xc0 <- array(NaN,c(12,N))
  Outcome.xab <- Stratified.xab.xa1 <- Stratified.xab.xa0 <- Stratified.xab.xb1 <- Stratified.xab.xb0 <- Stratified.xab.xc1 <- Stratified.xab.xc0 <- array(NaN,c(length(seq(-1,1,0.1)),12,N))
  Outcome.xac <- Stratified.xac.xa1 <- Stratified.xac.xa0 <- Stratified.xac.xb1 <- Stratified.xac.xb0 <- Stratified.xac.xc1 <- Stratified.xac.xc0 <- array(NaN,c(length(seq(-1,1,0.1)),12,N))
  Outcome.xbc <- Stratified.xbc.xa1 <- Stratified.xbc.xa0 <- Stratified.xbc.xb1 <- Stratified.xbc.xb0 <- Stratified.xbc.xc1 <- Stratified.xbc.xc0 <- array(NaN,c(length(seq(-1,1,0.1)),12,N))
  Outcome.xabc <-Stratified.xabc.xa1 <- Stratified.xabc.xa0 <- Stratified.xabc.xb1 <- Stratified.xabc.xb0 <- Stratified.xabc.xc1 <- Stratified.xabc.xc0 <- array(NaN,c(16,12,N))
  Outcome.xabc.Rho <-Stratified.xabc.Rho.xa1 <- Stratified.xabc.Rho.xa0 <- Stratified.xabc.Rho.xb1 <- Stratified.xabc.Rho.xb0 <- Stratified.xabc.Rho.xc1 <- Stratified.xabc.Rho.xc0 <- array(NaN,c(12,N))
  ATE.true <- data.frame(population=numeric(),rct.tau=numeric(),rct.obs=numeric(),obs.tau=numeric(),obs.obs=numeric())
  Diff.dis <- data.frame(x0=numeric(),xa=numeric(),xb=numeric(),xc=numeric(),xab=numeric(),xac=numeric(),xbc=numeric(),xabc=numeric(),xabc.rho=numeric())
  
  Methods.Est <- array(NaN,c(N,11))
  
  n <- 0
  while (n < N) {
    n <- n+1
    
    # step 1: DGP
    population <- DGP(Rho,n.pop,p.a.succ,p.b.succ,p.c.succ,Linear,Hetero,TauLinear,Weight.same)
    
    # step 2: sampling R=1/0
    dat.rct <- population[sample(1:n.pop,size = n.rct),]
    prob <- ProbUnitSelection(population,unitselect=unitselect,pR=pR)
    dat.obs <- population[sample(1:dim(population)[1],size = n.obs,prob = prob),]
    # A valid prob vector will be normalised to sum to 1 and used as sampling weights.
    #dat.rct <- population[sample(1:n.pop,size = n.rct,prob=prob),]
    #dat.obs <- population[!(rownames(population) %in% rownames(dat.rct)),]
    
    # step 3: treatment assignment Z=1/0; non-randomization induced confounding
    pi <- inv.logit((lambda*(dat.obs$Mu-mean(dat.obs$Mu)))/sd(dat.obs$Mu))
    dat.obs$Z <- rbinom(n.obs,1,pi)
    dat.obs$Y <- dat.obs$Mu+dat.obs$Z*dat.obs$tau
    dat.rct$Z <- rbinom(n.rct,1,0.5)
    dat.rct$Y <- dat.rct$Mu+dat.rct$Z*dat.rct$tau
    
    # the estimate is average treatment effect; return conterfactual computation for each unit
    xvars <- c('xA','xB','xC')
    dat.obs <- potential_outcomes_direct_ipw_dr(dat.obs,xvars,ntree,ndpost)
    
    
    ATE.true[n,'population'] <- mean(population$tau); ATE.true[n,'rct.tau'] <- mean(dat.rct$tau); 
    ATE.true[n,'rct.obs'] <- mean(dat.rct[dat.rct$Z==1,'Y']) - mean(dat.rct[dat.rct$Z==0,'Y']); 
    ATE.true[n,'obs.tau'] <- mean(dat.obs$tau); ATE.true[n,'obs.obs'] <- mean(dat.obs[dat.obs$Z==1,'Y']) - mean(dat.obs[dat.obs$Z==0,'Y'])
    
    pa <- (table(dat.rct[,1])/dim(dat.rct)[1])[2] #success probablity of xA in rct
    pb <- (table(dat.rct[,2])/dim(dat.rct)[1])[2] #success probablity of xB in rct
    pc <- (table(dat.rct[,3])/dim(dat.rct)[1])[2] #success probablity of xC in rct
    
    # ================== Methods comparsion: ate obtained from each methods
    Methods.Est[n,] <- c(mean(dat.obs$tau),ResultsEstimators(dat.obs))
    # ================== 
    # step 1: no subsampling 
    Diff.dis[n,'x0'] <- DKL(dat.rct,dat.obs)
    results.stra <- GenerateStrataResults(dat.rct,dat.obs)
    Stratified.0.xa1[,n] <- results.stra[1,]; Stratified.0.xa0[,n] <- results.stra[2,]
    Stratified.0.xb1[,n] <- results.stra[3,]; Stratified.0.xb0[,n] <- results.stra[4,]
    Stratified.0.xc1[,n] <- results.stra[5,]; Stratified.0.xc0[,n] <- results.stra[6,]
    
    Outcome.0[,n] <- GenerateOverallResults(dat.rct,dat.obs)
    # ==================
    # subsampling based on one confounding, xA
    n.success <- dim(dat.obs)[1]*pa
    dat.0 <- dat.obs[dat.obs$xA==0,]; dat.0 <- dat.0[sample(1:dim(dat.0)[1],size = dim(dat.obs)[1] - n.success ,replace = TRUE),]
    dat.1 <- dat.obs[dat.obs$xA==1,]; dat.1 <- dat.1[sample(1:dim(dat.1)[1],size = n.success ,replace = TRUE),]
    dat.obs.sub <- rbind(dat.0,dat.1)
    Diff.dis[n,'xa'] <-  DKL(dat.rct,dat.obs.sub)
    
    results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
    Stratified.xa.xa1[,n] <- results.stra[1,]; Stratified.xa.xa0[,n] <- results.stra[2,]
    Stratified.xa.xb1[,n] <- results.stra[3,]; Stratified.xa.xb0[,n] <- results.stra[4,]
    Stratified.xa.xc1[,n] <- results.stra[5,]; Stratified.xa.xc0[,n] <- results.stra[6,]
    
    Outcome.xa[,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
    # ==================
    # subsampling based on one confounding, xB
    n.success <- dim(dat.obs)[1]*pb
    dat.0 <- dat.obs[dat.obs$xB==0,]; dat.0 <- dat.0[sample(1:dim(dat.0)[1],size = dim(dat.obs)[1] - n.success ,replace = TRUE),]
    dat.1 <- dat.obs[dat.obs$xB==1,]; dat.1 <- dat.1[sample(1:dim(dat.1)[1],size = n.success ,replace = TRUE),]
    dat.obs.sub <- rbind(dat.0,dat.1)
    Diff.dis[n,'xb'] <-  DKL(dat.rct,dat.obs.sub)
    
    results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
    Stratified.xb.xa1[,n] <- results.stra[1,]; Stratified.xb.xa0[,n] <- results.stra[2,]
    Stratified.xb.xb1[,n] <- results.stra[3,]; Stratified.xb.xb0[,n] <- results.stra[4,]
    Stratified.xb.xc1[,n] <- results.stra[5,]; Stratified.xb.xc0[,n] <- results.stra[6,]
    
    Outcome.xb[,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
    # ==================
    # subsampling based on one confounding, xC
    n.success <- dim(dat.obs)[1]*pc
    dat.0 <- dat.obs[dat.obs$xC==0,]; dat.0 <- dat.0[sample(1:dim(dat.0)[1],size = dim(dat.obs)[1] - n.success ,replace = TRUE),]
    dat.1 <- dat.obs[dat.obs$xC==1,]; dat.1 <- dat.1[sample(1:dim(dat.1)[1],size = n.success ,replace = TRUE),]
    dat.obs.sub <- rbind(dat.0,dat.1)
    Diff.dis[n,'xc'] <-  DKL(dat.rct,dat.obs.sub)
    
    results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
    Stratified.xc.xa1[,n] <- results.stra[1,]; Stratified.xc.xa0[,n] <- results.stra[2,]
    Stratified.xc.xb1[,n] <- results.stra[3,]; Stratified.xc.xb0[,n] <- results.stra[4,]
    Stratified.xc.xc1[,n] <- results.stra[5,]; Stratified.xc.xc0[,n] <- results.stra[6,]
    
    Outcome.xc[,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
    
    # ==================
    # subsampling based on two confounding: xA,xB
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      
      rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
      pvars <- pnorm(rawvars)
      xA <- qbinom(pvars[,1],1,pa)
      xB <- qbinom(pvars[,2],1,pb)
      dat.00 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0),]; dat.00 <- dat.00[sample(1:dim(dat.00)[1],size = sum((xA==0)&(xB==0)),replace = TRUE),]
      dat.10 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0),]; dat.10 <- dat.10[sample(1:dim(dat.10)[1],size = sum((xA==1)&(xB==0)),replace = TRUE),]
      dat.01 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1),]; dat.01 <- dat.01[sample(1:dim(dat.01)[1],size = sum((xA==0)&(xB==1)),replace = TRUE),]
      dat.11 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1),]; dat.11 <- dat.11[sample(1:dim(dat.11)[1],size = sum((xA==1)&(xB==1)),replace = TRUE),]
      dat.obs.sub <- rbind(dat.00,dat.01,dat.10,dat.11)
      diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
      # table(dat.obs.sub[,1])/dim(dat.obs.sub)[1]
      # table(dat.rct[,1])/dim(dat.rct)[1]
      
      results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
      Stratified.xab.xa1[i,,n] <- results.stra[1,]; Stratified.xab.xa0[i,,n] <- results.stra[2,]
      Stratified.xab.xb1[i,,n] <- results.stra[3,]; Stratified.xab.xb0[i,,n] <- results.stra[4,]
      Stratified.xab.xc1[i,,n] <- results.stra[5,]; Stratified.xab.xc0[i,,n] <- results.stra[6,]
      
      Outcome.xab[i,,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
      i <- i+1
    }
    Diff.dis[n,'xab'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    # ==================
    # subsampling based on two confounding: xA,xC
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      
      rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
      pvars <- pnorm(rawvars)
      xA <- qbinom(pvars[,1],1,pa)
      xC <- qbinom(pvars[,2],1,pc)
      dat.00 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xC==0),]; dat.00 <- dat.00[sample(1:dim(dat.00)[1],size = sum((xA==0)&(xC==0)),replace = TRUE),]
      dat.10 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xC==0),]; dat.10 <- dat.10[sample(1:dim(dat.10)[1],size = sum((xA==1)&(xC==0)),replace = TRUE),]
      dat.01 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xC==1),]; dat.01 <- dat.01[sample(1:dim(dat.01)[1],size = sum((xA==0)&(xC==1)),replace = TRUE),]
      dat.11 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xC==1),]; dat.11 <- dat.11[sample(1:dim(dat.11)[1],size = sum((xA==1)&(xC==1)),replace = TRUE),]
      dat.obs.sub <- rbind(dat.00,dat.01,dat.10,dat.11)
      diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
      
      results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
      Stratified.xac.xa1[i,,n] <- results.stra[1,]; Stratified.xac.xa0[i,,n] <- results.stra[2,]
      Stratified.xac.xb1[i,,n] <- results.stra[3,]; Stratified.xac.xb0[i,,n] <- results.stra[4,]
      Stratified.xac.xc1[i,,n] <- results.stra[5,]; Stratified.xac.xc0[i,,n] <- results.stra[6,]
      
      Outcome.xac[i,,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
      i <- i+1
    }
    Diff.dis[n,'xac'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    # ==================
    # subsampling based on two confounding: xB,xC
    mu <- rep(0, 2)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      Sigma <- matrix(c(1,rho,
                        rho, 1), nrow = 2, ncol = 2)
      
      rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
      pvars <- pnorm(rawvars)
      xB <- qbinom(pvars[,1],1,pb)
      xC <- qbinom(pvars[,2],1,pc)
      dat.00 <- dat.obs[(dat.obs$xB==0)&(dat.obs$xC==0),]; dat.00 <- dat.00[sample(1:dim(dat.00)[1],size = sum((xB==0)&(xC==0)),replace = TRUE),]
      dat.10 <- dat.obs[(dat.obs$xB==1)&(dat.obs$xC==0),]; dat.10 <- dat.10[sample(1:dim(dat.10)[1],size = sum((xB==1)&(xC==0)),replace = TRUE),]
      dat.01 <- dat.obs[(dat.obs$xB==0)&(dat.obs$xC==1),]; dat.01 <- dat.01[sample(1:dim(dat.01)[1],size = sum((xB==0)&(xC==1)),replace = TRUE),]
      dat.11 <- dat.obs[(dat.obs$xB==1)&(dat.obs$xC==1),]; dat.11 <- dat.11[sample(1:dim(dat.11)[1],size = sum((xB==1)&(xC==1)),replace = TRUE),]
      dat.obs.sub <- rbind(dat.00,dat.01,dat.10,dat.11)
      diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
      # table(dat.obs.sub[,1])/dim(dat.obs.sub)[1]
      # table(dat.rct[,1])/dim(dat.rct)[1]
      
      results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
      Stratified.xbc.xa1[i,,n] <- results.stra[1,]; Stratified.xbc.xa0[i,,n] <- results.stra[2,]
      Stratified.xbc.xb1[i,,n] <- results.stra[3,]; Stratified.xbc.xb0[i,,n] <- results.stra[4,]
      Stratified.xbc.xc1[i,,n] <- results.stra[5,]; Stratified.xbc.xc0[i,,n] <- results.stra[6,]
      
      Outcome.xbc[i,,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
      i <- i+1
    }
    Diff.dis[n,'xbc'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    # ==================
    # subsampling based on three confounding: xA,xB,xC
    mu <- rep(0, 3)
    i <- 1
    diff.dis.mean <- NULL
    for (rho in seq(-1,1,0.1)){
      try(
        {
          Sigma <- matrix(c(1,rho,rho,
                            rho, 1,rho,
                            rho,rho,1), nrow = 3, ncol = 3)
          
          rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
          pvars <- pnorm(rawvars)
          xA <- qbinom(pvars[,1],1,pa)
          xB <- qbinom(pvars[,2],1,pb)
          xC <- qbinom(pvars[,3],1,pc)
          dat.000 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
          dat.000 <- dat.000[sample(1:dim(dat.000)[1],size = sum((xA==0)&(xB==0)&(xC==0)),replace = TRUE),]
          dat.001 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
          dat.001 <- dat.001[sample(1:dim(dat.001)[1],size = sum((xA==0)&(xB==0)&(xC==1)),replace = TRUE),]
          dat.010 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
          dat.010 <- dat.010[sample(1:dim(dat.010)[1],size = sum((xA==0)&(xB==1)&(xC==0)),replace = TRUE),]
          dat.011 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
          dat.011 <- dat.011[sample(1:dim(dat.011)[1],size = sum((xA==0)&(xB==1)&(xC==1)),replace = TRUE),]
          dat.100 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
          dat.100 <- dat.100[sample(1:dim(dat.100)[1],size = sum((xA==1)&(xB==0)&(xC==0)),replace = TRUE),]
          dat.101 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
          dat.101 <- dat.101[sample(1:dim(dat.101)[1],size = sum((xA==1)&(xB==0)&(xC==1)),replace = TRUE),]
          dat.110 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
          dat.110 <- dat.110[sample(1:dim(dat.110)[1],size = sum((xA==1)&(xB==1)&(xC==0)),replace = TRUE),]
          dat.111 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
          dat.111 <- dat.111[sample(1:dim(dat.111)[1],size = sum((xA==1)&(xB==1)&(xC==1)),replace = TRUE),]
          
          dat.obs.sub <- rbind(dat.000,dat.001,dat.010,dat.011,dat.100,dat.101,dat.110,dat.111)
          diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
          
          results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
          Stratified.xabc.xa1[i,,n] <- results.stra[1,]; Stratified.xabc.xa0[i,,n] <- results.stra[2,]
          Stratified.xabc.xb1[i,,n] <- results.stra[3,]; Stratified.xabc.xb0[i,,n] <- results.stra[4,]
          Stratified.xabc.xc1[i,,n] <- results.stra[5,]; Stratified.xabc.xc0[i,,n] <- results.stra[6,]
          
          Outcome.xabc[i,,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
          i <- i+1
        },
        silent = TRUE
      )
      
    }
    Diff.dis[n,'xabc'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
    
    
    # ==================
    # subsampling based on three confounding and correct dependency Rho: xA,xB,xC
    rhoAB <- simstudy:::.findRhoBin(p1 = pa, 
                                    p2 = pb, d = cor(dat.rct$xA,dat.rct$xB))
    rhoAC <- simstudy:::.findRhoBin(p1 = pa, 
                                    p2 = pc, d = cor(dat.rct$xA,dat.rct$xC))
    rhoBC <- simstudy:::.findRhoBin(p1 = pb, 
                                    p2 = pc, d = cor(dat.rct$xB,dat.rct$xC))
    
    mu <- rep(0, 3)
    
    Sigma <- matrix(c(1,rhoAB,rhoAC,
                      rhoAB,1,rhoBC,
                      rhoAC,rhoBC,1), nrow = 3, ncol = 3)
    
    rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
    pvars <- pnorm(rawvars)
    xA <- 1*(pvars[,1] < pa)
    xB <- 1*(pvars[,2] < pb)
    xC <- 1*(pvars[,3] < pc)
    dat.000 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
    dat.000 <- dat.000[sample(1:dim(dat.000)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==0)&(dat.rct$xC==0)),replace = TRUE),]
    dat.001 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
    dat.001 <- dat.001[sample(1:dim(dat.001)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==0)&(dat.rct$xC==1)),replace = TRUE),]
    dat.010 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
    dat.010 <- dat.010[sample(1:dim(dat.010)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==1)&(dat.rct$xC==0)),replace = TRUE),]
    dat.011 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
    dat.011 <- dat.011[sample(1:dim(dat.011)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==1)&(dat.rct$xC==1)),replace = TRUE),]
    dat.100 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
    dat.100 <- dat.100[sample(1:dim(dat.100)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==0)&(dat.rct$xC==0)),replace = TRUE),]
    dat.101 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
    dat.101 <- dat.101[sample(1:dim(dat.101)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==0)&(dat.rct$xC==1)),replace = TRUE),]
    dat.110 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
    dat.110 <- dat.110[sample(1:dim(dat.110)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==1)&(dat.rct$xC==0)),replace = TRUE),]
    dat.111 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
    dat.111 <- dat.111[sample(1:dim(dat.111)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==1)&(dat.rct$xC==1)),replace = TRUE),]
    
    dat.obs.sub <- rbind(dat.000,dat.001,dat.010,dat.011,dat.100,dat.101,dat.110,dat.111)
    Diff.dis[n,'xabc.rho'] <- DKL(dat.rct,dat.obs.sub)
    
    results.stra <- GenerateStrataResults(dat.rct,dat.obs.sub)
    Stratified.xabc.Rho.xa1[,n] <- results.stra[1,]; Stratified.xabc.Rho.xa0[,n] <- results.stra[2,]
    Stratified.xabc.Rho.xb1[,n] <- results.stra[3,]; Stratified.xabc.Rho.xb0[,n] <- results.stra[4,]
    Stratified.xabc.Rho.xc1[,n] <- results.stra[5,]; Stratified.xabc.Rho.xc0[,n] <- results.stra[6,]
    
    Outcome.xabc.Rho[,n] <- GenerateOverallResults(dat.rct,dat.obs.sub)
    
  }
  
  results.0 <- PrintResults4OverallStratified(  Outcome.0,
                                                Stratified.0.xa1,Stratified.0.xa0,
                                                Stratified.0.xb1,Stratified.0.xb0,
                                                Stratified.0.xc1,Stratified.0.xc0,
                                                nsubvar=0)
  
  results.xa  <- PrintResults4OverallStratified(Outcome.xa,
                                                Stratified.xa.xa1,Stratified.xa.xa0,
                                                Stratified.xa.xb1,Stratified.xa.xb0,
                                                Stratified.xa.xc1,Stratified.xa.xc0,
                                                nsubvar=1)
  
  results.xb  <- PrintResults4OverallStratified(Outcome.xb,
                                                Stratified.xb.xa1,Stratified.xb.xa0,
                                                Stratified.xb.xb1,Stratified.xb.xb0,
                                                Stratified.xb.xc1,Stratified.xb.xc0,
                                                nsubvar=1)
  
  results.xc  <- PrintResults4OverallStratified(Outcome.xc,
                                                Stratified.xc.xa1,Stratified.xc.xa0,
                                                Stratified.xc.xb1,Stratified.xc.xb0,
                                                Stratified.xc.xc1,Stratified.xc.xc0,
                                                nsubvar=1)
  
  results.xab <- PrintResults4OverallStratified(Outcome.xab,
                                                Stratified.xab.xa1,Stratified.xab.xa0,
                                                Stratified.xab.xb1,Stratified.xab.xb0,
                                                Stratified.xab.xc1,Stratified.xab.xc0,
                                                nsubvar=2)
  
  results.xac <- PrintResults4OverallStratified(Outcome.xac,
                                                Stratified.xac.xa1,Stratified.xac.xa0,
                                                Stratified.xac.xb1,Stratified.xac.xb0,
                                                Stratified.xac.xc1,Stratified.xac.xc0,
                                                nsubvar=2)
  
  results.xbc <- PrintResults4OverallStratified(Outcome.xbc,
                                                Stratified.xbc.xa1,Stratified.xbc.xa0,
                                                Stratified.xbc.xb1,Stratified.xbc.xb0,
                                                Stratified.xbc.xc1,Stratified.xbc.xc0,
                                                nsubvar=2)
  
  results.xabc <- PrintResults4OverallStratified(Outcome.xabc,
                                                 Stratified.xabc.xa1,Stratified.xabc.xa0,
                                                 Stratified.xabc.xb1,Stratified.xabc.xb0,
                                                 Stratified.xabc.xc1,Stratified.xabc.xc0,
                                                 nsubvar=3)
  
  results.xabc.Rho <- PrintResults4OverallStratified(Outcome.xabc.Rho,
                                                     Stratified.xabc.Rho.xa1,Stratified.xabc.xa0,
                                                     Stratified.xabc.Rho.xb1,Stratified.xabc.xb0,
                                                     Stratified.xabc.Rho.xc1,Stratified.xabc.xc0,
                                                     nsubvar=1)
  
  #Outcome <- GenerateOverallResults(dat.rct,dat.obs.sub)
  #round(PrintDistance(Outcome,nsubvar),digits = 3)
  
  Results <- cbind(results.0,results.xa,results.xb,results.xc,
                   results.xab,results.xac,results.xbc,
                   results.xabc,results.xabc.Rho)
  rownames(Results) <- c('Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                         'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  #Results <- rbind(results.0,0,results.xa,0,results.xb,0,results.xc,0,
  #                 results.xab,0,results.xac,0,results.xbc,0,results.xabc,0,results.xabc.Rho,0)
  
  #Methods.Est <- as.data.frame(Methods.Est)
  #colnames(Methods.Est) <- c('true','naive','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART','DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  
  return(list(Methods.Est,Results,ATE.true,Diff.dis))
  
}





## ===============
BiasEmpSeMSE <- function(MethodsEst){
  nsim <- dim(MethodsEst)[1]
  true <- MethodsEst[,1]
  theta.hat <- MethodsEst[,2:ncol(MethodsEst)]
  bias <- colMeans(sweep(theta.hat,1,true))
  theta.hat.mean <- colMeans(theta.hat)
  bias.mcse <- sqrt(colSums(sweep(theta.hat,2,theta.hat.mean)^2)/(nsim*(nsim-1)))
  EmpSe <- sqrt(colSums(sweep(theta.hat,2,theta.hat.mean)^2)/(nsim-1))
  Empse.mcse <- EmpSe/sqrt(2*(nsim-1))
  MSE <- colSums(sweep(theta.hat,1,true)^2)/nsim
  MSE.mcse <- sqrt(colSums(sweep((sweep(theta.hat,1,true)^2),2,MSE)^2)/(nsim*(nsim-1)))
  Performance <- cbind(round(bias,digits = 5),'&',
                       round(EmpSe,digits = 5),'&',
                       round(MSE,digits = 5),'\\') 
  colnames(Performance) <- c('Bias','',
                             'EmpSE','',
                             'MSE','')
  rownames(Performance) <- c('multirow{11}{*}{$X_{A}$}&Crude&','&Rg/LReg&','&Rg/BART&','&Rg/ps-BART&','&IPW/LReg&','&IPW/BART&','&DR/LReg-LReg&','&DR/LReg-BART&','&DR/BART-LReg&','&DR/BART-BART&')
  return(Performance)
}

SelectEstimator <- function(dat,method){
  estimate <- dat[rownames(dat)==method,]
  #estimate <- estimate[!is.na(estimate)]
  estimate <- c(dat[1,1],estimate) #dat[1,1] is crude without sampling; Z not control, R not control
  return(estimate)
}

PlotEstimatorSubsampling <- function(ResultsDistanceSubsampling,title){
  RgLReg <- SelectEstimator(ResultsDistanceSubsampling,'Rg/LReg')
  RgBART <- SelectEstimator(ResultsDistanceSubsampling,'Rg/BART')
  RgpsBART <- SelectEstimator(ResultsDistanceSubsampling,'Rg/ps-BART')
  IPWLReg <- SelectEstimator(ResultsDistanceSubsampling,'IPW/LReg')
  IPWBART <- SelectEstimator(ResultsDistanceSubsampling,'IPW/BART')
  DRLRegLReg <- SelectEstimator(ResultsDistanceSubsampling,'DR/LReg-LReg')
  DRLRegBART <- SelectEstimator(ResultsDistanceSubsampling,'DR/LReg-BART')
  DRBARTLReg <- SelectEstimator(ResultsDistanceSubsampling,'DR/BART-LReg')
  DRBARTBART <- SelectEstimator(ResultsDistanceSubsampling,'DR/BART-BART')
  
  #x <- seq(10)
  x <- factor(c('T-S-','T+S-','T+S1','T+S2','T+S3','T+S12','T+S13','T+S23','T+S123','T+Sx'), levels = c('T-S-','T+S-','T+S1','T+S2','T+S3','T+S12','T+S13','T+S23','T+S123','T+Sx'))
  data <- data.frame(x=x,RgLReg=RgLReg,RgBART=RgBART,RgpsBART=RgpsBART,IPWLReg=IPWLReg,IPWBART=IPWBART,DRLRegLReg=DRLRegLReg,DRLRegBART=DRLRegBART,DRBARTLReg=DRBARTLReg,DRBARTBART=DRBARTBART)
  
  
  library(ggplot2)
  colors <- c("Rg/LReg" = "yellow", "Rg/BART" = "yellow3", "Rg/ps-BART" = "yellow4", 
              "IPW/LReg" = "tomato","IPW/BART" = "tomato4",
              "DR/LReg-LReg" = "springgreen4", "DR/LReg-BART" = "springgreen3",
              "DR/BART-LReg" = "turquoise4", "DR/BART-BART" = "turquoise"
  )
  
  p <-ggplot(data=data) +
    geom_line(aes(x=x, y=RgLReg,color="Rg/LReg",group=1),size=1) +
    geom_line(aes(x=x, y=RgBART,color="Rg/BART",group=1),size=1) +
    geom_line(aes(x=x, y=RgpsBART,color="Rg/ps-BART",group=1),size=1) +
    geom_line(aes(x=x, y=IPWLReg,color="IPW/LReg",group=1),size=1) +
    geom_line(aes(x=x, y=IPWBART,color="IPW/BART",group=1),size=1) +
    geom_line(aes(x=x, y=DRLRegLReg,color="DR/LReg-LReg",group=1),size=1) +
    geom_line(aes(x=x, y=DRLRegBART,color="DR/LReg-BART",group=1),size=1) +
    geom_line(aes(x=x, y=DRBARTLReg,color="DR/BART-LReg",group=1),size=1) +
    geom_line(aes(x=x, y=DRBARTBART,color="DR/BART-BART",group=1),size=1) +
    labs(x = title,
         y = bquote(D[1](hat(tau)^{rct}*","*hat(tau)^{obs})),
         color = "Legend") +
    scale_color_manual(values = colors)
  
  return(p)
}

FourSamplingPlot <- function(Estimation.1,Estimation.2,Estimation.3,Estimation.4){
  p1 <- PlotEstimatorSubsampling(Estimation.1,title=bquote(pi[s](x[1])))
  p2 <- PlotEstimatorSubsampling(Estimation.2,title=bquote(pi[s](x[2])))
  p3 <- PlotEstimatorSubsampling(Estimation.3,title=bquote(pi[s](x[3])))
  p4 <- PlotEstimatorSubsampling(Estimation.4,title=bquote(pi[s](bold(x))))
  ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend = TRUE,legend = "right")
}

PrintDistanceLatex <- function(dat){
  dat.latex <- NULL
  for (i in seq(dim(dat)[2])){
    aa <- cbind('&',dat[,i])
    dat.latex <- cbind(dat.latex,aa)
  }
  dat.latex <- cbind(dat.latex,'\\')
  noquote(dat.latex)
}



FourSamplingBoxPlot <- function(Estimation.S1,
                                Estimation.S2,
                                Estimation.S3,
                                Estimation.S123){
  
  p1 <- PlotOneDGM(Estimation.S1,title=bquote(pi[s](x[1])))
  p2 <- PlotOneDGM(Estimation.S2,title=bquote(pi[s](x[2])))
  p3 <- PlotOneDGM(Estimation.S3,title=bquote(pi[s](x[3])))
  p4 <- PlotOneDGM(Estimation.S123,title=bquote(pi[s](bold(x))))
  p <- ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend = TRUE,legend = "right")
  return(p)
}

FourSamplingBoxPlot_no_color <- function(Estimation.S1,
                                Estimation.S2,
                                Estimation.S3,
                                Estimation.S123){
  
  p1 <- PlotOneDGM_no_color(Estimation.S1,title=bquote(pi[s](x[1])))
  p2 <- PlotOneDGM_no_color(Estimation.S2,title=bquote(pi[s](x[2])))
  p3 <- PlotOneDGM_no_color(Estimation.S3,title=bquote(pi[s](x[3])))
  p4 <- PlotOneDGM_no_color(Estimation.S123,title=bquote(pi[s](bold(x))))
  p <- ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend = TRUE,legend = "right")
  return(p)
}

FourSamplingBoxPlot_simple <- function(Estimation.S1,
                                Estimation.S2,
                                Estimation.S3,
                                Estimation.S123){
  
  p1 <- PlotOneDGM_simple(Estimation.S1,title=bquote(beta[11]))
  p2 <- PlotOneDGM_simple(Estimation.S2,title=bquote(beta[12]))
  p3 <- PlotOneDGM_simple(Estimation.S3,title=bquote(beta[21]))
  p4 <- PlotOneDGM_simple(Estimation.S123,title=bquote(beta[22]))
  p <- ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend = TRUE,legend = "none")
  return(p)
}

get_mean_var <- function(df){
  n.dim <- length(dim(df))
  if(n.dim==2){
    df.mean <- apply(df, 1, mean)
    df.var <- apply(df, 1, var)
  } else {
    df <- apply(df, c(2,3), mean)
    df.mean <- apply(df, 1, mean)
    df.var <- apply(df,1,var)
  }
  return(list(df.mean, df.var))
}

PlotOneDGM <- function(Estimation.S1,title){
  Outcome.S1.w0 <- get_mean_var(Estimation.S1[[1]])
  Outcome.S1.w1 <- get_mean_var(Estimation.S1[[2]])
  Outcome.S1.w2 <- get_mean_var(Estimation.S1[[3]])
  Outcome.S1.w3 <- get_mean_var(Estimation.S1[[4]])
  Outcome.S1.w12 <- get_mean_var(Estimation.S1[[5]])
  Outcome.S1.w13 <- get_mean_var(Estimation.S1[[6]])
  Outcome.S1.w23 <- get_mean_var(Estimation.S1[[7]])
  Outcome.S1.w123 <- get_mean_var(Estimation.S1[[8]])
  Outcome.S1.w123true <- get_mean_var(Estimation.S1[[9]])
  
  Outcome.S1.mean <- rbind(Outcome.S1.w0[[1]],
                           Outcome.S1.w1[[1]],
                           Outcome.S1.w2[[1]],
                           Outcome.S1.w3[[1]],
                           Outcome.S1.w12[[1]],
                           Outcome.S1.w13[[1]],
                           Outcome.S1.w23[[1]],
                           Outcome.S1.w123[[1]],
                           Outcome.S1.w123true[[1]])
  
  rownames(Outcome.S1.mean) <- c("1","X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3","Xi")
  colnames(Outcome.S1.mean) <- c("pop.true","rct.est","obs.true",'Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                                 'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  outcome.S1.mean <- Outcome.S1.mean[,-c(1,2,3)]
  outcome.S1.mean <- melt(outcome.S1.mean)
  colnames(outcome.S1.mean) <- c("weight","tau","ATE")
  estimator_model <- outcome.S1.mean$tau
  estimator <- sapply(str_split(estimator_model, "/"),"[[",1)
  outcome.S1.mean[,"estimator"] <- estimator
  levels.order <-  c(1,"X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3","Xi")
  outcome.S1.mean$weight <- factor(outcome.S1.mean$weight, levels = levels.order)
  rct.est <- mean(Outcome.S1.mean[,2])
  
  p <- ggplot(data = outcome.S1.mean %>% filter(!weight %in% c(1,"Xi")), aes(x = ATE, y = weight,colour=estimator)) +
    scale_y_discrete(limits=c(1,"X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3","Xi")) +
    geom_point(data=outcome.S1.mean%>%filter(weight %in% c(1,"Xi")), aes(x=ATE, y=weight,fill=estimator),colour="black", size=4,shape = 21, stroke=2, show.legend = FALSE) +
    geom_point(position = position_dodge(0.5), size=4) +
    geom_vline(xintercept=rct.est,linetype="dashed",size=1) +
    theme(axis.text.y = element_text(size = 10,face="bold"))  +
    #theme(panel.grid = element_line(linetype = "dotted",colour="black")) +
    ggtitle(title) 
  
  return(p)
}

PlotOneDGM_no_color <- function(Estimation.S1,title){
  #browser()
  Outcome.S1.w1 <- get_mean_var(Estimation.S1[[2]])
  Outcome.S1.w2 <- get_mean_var(Estimation.S1[[3]])
  Outcome.S1.w3 <- get_mean_var(Estimation.S1[[4]])
  Outcome.S1.w12 <- get_mean_var(Estimation.S1[[5]])
  Outcome.S1.w13 <- get_mean_var(Estimation.S1[[6]])
  Outcome.S1.w23 <- get_mean_var(Estimation.S1[[7]])
  Outcome.S1.w123 <- get_mean_var(Estimation.S1[[8]])

  Outcome.S1.mean <- rbind(Outcome.S1.w1[[1]],
                           Outcome.S1.w2[[1]],
                           Outcome.S1.w3[[1]],
                           Outcome.S1.w12[[1]],
                           Outcome.S1.w13[[1]],
                           Outcome.S1.w23[[1]],
                           Outcome.S1.w123[[1]])
  
  Outcome.S1.var <- rbind(Outcome.S1.w1[[2]],
                          Outcome.S1.w2[[2]],
                          Outcome.S1.w3[[2]],
                          Outcome.S1.w12[[2]],
                          Outcome.S1.w13[[2]],
                          Outcome.S1.w23[[2]],
                          Outcome.S1.w123[[2]])
  
  rownames(Outcome.S1.mean) <- c("X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3")
  colnames(Outcome.S1.mean) <- c("pop.true","rct.est","obs.true",'Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                                 'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  outcome.S1.mean <- Outcome.S1.mean[,c("Rg/ps-BART")]
  outcome.S1.mean <- data.frame(weight = c("X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3"),
                                ATE=outcome.S1.mean)
  levels.order <-  c("X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3")
  outcome.S1.mean$weight <- factor(outcome.S1.mean$weight, levels = levels.order)
  rct.est <- mean(Outcome.S1.mean[,2])

  rownames(Outcome.S1.var) <- c("X1","X2","X3","X1,X2","X1,X3","X2,X3","X1,X2,X3")
  colnames(Outcome.S1.var) <- c("pop.true","rct.est","obs.true",'Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                                'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  outcome.S1.var <- Outcome.S1.var[,c("Rg/ps-BART")]
  outcome.S1.mean$lb <- outcome.S1.mean$ATE - sqrt(outcome.S1.var)
  outcome.S1.mean$ub <- outcome.S1.mean$ATE + sqrt(outcome.S1.var)
  
  p <- ggplot(data = outcome.S1.mean, aes(x = ATE, y= weight)) +
    geom_errorbar(aes(xmin=lb, xmax=ub), width=0.2) +
    geom_point(position = position_dodge2(width = 1)) +
    geom_vline(xintercept=rct.est,linetype="dashed",size=1) +
    theme(axis.text.y = element_text(size = 12,face="bold"),
          axis.text.x = element_text(size = 10,face="bold"),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))  +
    ggtitle(title)  
  
  return(p)
}



PlotOneDGM_simple <- function(Estimation.S1,title){
  #browser()
  Outcome.S1.w0 <- get_mean_var(Estimation.S1[[1]])
  Outcome.S1.w123true <- get_mean_var(Estimation.S1[[9]])
  
  Outcome.S1.mean <- rbind(Outcome.S1.w0[[1]],
                           Outcome.S1.w123true[[1]])
  
  Outcome.S1.var <- rbind(Outcome.S1.w0[[2]],
                           Outcome.S1.w123true[[2]])
  
  rownames(Outcome.S1.mean) <- c("1","Xi")
  colnames(Outcome.S1.mean) <- c("pop.true","rct.est","obs.true",'Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                                 'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  outcome.S1.mean <- Outcome.S1.mean[,-c(1,2,3)]
  outcome.S1.mean <- melt(outcome.S1.mean)
  colnames(outcome.S1.mean) <- c("weight","tau","ATE")
  estimator_model <- outcome.S1.mean$tau
  estimator <- sapply(str_split(estimator_model, "/"),"[[",1)
  outcome.S1.mean[,"estimator"] <- estimator
  levels.order <-  c(1,"Xi")
  outcome.S1.mean$weight <- factor(outcome.S1.mean$weight, levels = levels.order)
  rct.est <- mean(Outcome.S1.mean[,2])
  outcome.S1.mean$weight <- ifelse(outcome.S1.mean$weight %in% c(1),"unweighted","weighted")
  outcome.S1.mean$estimator <- ifelse(outcome.S1.mean$estimator %in% c("Crude"),"unadjusted","adjusted")
  outcome.S1.mean$estimator <- apply(outcome.S1.mean[ , c("weight","estimator")] , 1 , paste , collapse = "+" ) 
  
  
  rownames(Outcome.S1.var) <- c("1","Xi")
  colnames(Outcome.S1.var) <- c("pop.true","rct.est","obs.true",'Crude','Rg/LReg','Rg/BART','Rg/ps-BART','IPW/LReg','IPW/BART',
                                 'DR/LReg-LReg','DR/LReg-BART','DR/BART-LReg','DR/BART-BART')
  outcome.S1.var <- Outcome.S1.var[,-c(1,2,3)]
  outcome.S1.var <- melt(outcome.S1.var)
  colnames(outcome.S1.var) <- c("weight","tau","var")
  estimator_model <- outcome.S1.var$tau
  estimator <- sapply(str_split(estimator_model, "/"),"[[",1)
  outcome.S1.var[,"estimator"] <- estimator
  levels.order <-  c(1,"Xi")
  outcome.S1.var$weight <- factor(outcome.S1.var$weight, levels = levels.order)
  rct.est <- mean(outcome.S1.var[,2])
  outcome.S1.var$weight <- ifelse(outcome.S1.var$weight %in% c(1),"S-","S+")
  outcome.S1.var$estimator <- ifelse(outcome.S1.var$estimator %in% c("Crude"),"T-","T+")
  outcome.S1.var$estimator <- apply(outcome.S1.var[ , c("weight","estimator")] , 1 , paste , collapse = "" ) 
  outcome.S1.mean$estimator <- outcome.S1.var$estimator
  outcome.S1.mean$lb <- outcome.S1.mean$ATE - sqrt(outcome.S1.var$var)
  outcome.S1.mean$ub <- outcome.S1.mean$ATE + sqrt(outcome.S1.var$var)
  rct.est <- mean(Outcome.S1.mean[,2])
  
  p <- ggplot(data = outcome.S1.mean, aes(x = ATE, y= estimator)) +
    geom_errorbar(data = outcome.S1.mean %>% filter(estimator %in% c("S+T-","S-T-")),aes(xmin=lb, xmax=ub), position = position_dodge2(width = 1), width=0.2) +
    geom_errorbar(data = outcome.S1.mean %>% filter(!(estimator %in% c("S+T-","S-T-"))),aes(xmin=lb, xmax=ub), position = position_dodge2(width = 1.5), width=1) +
    geom_point(position = position_dodge2(width = 1)) +
    geom_vline(xintercept=rct.est,linetype="dashed",size=1) +
    theme(axis.text.y = element_text(size = 12,face="bold"),
          axis.text.x = element_text(size = 10,face="bold"),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))  +
    ggtitle(title) 
  
  return(p)
}

# p.a.succ <- p.b.succ <- p.c.succ <- 0.5
# n.pop <- 1e3
# n.rct <- 250
# n.obs <- 750
# lambda <- -log(3)
# N <- 250 # the number of times for DGP
# ntree <- 50
# ndpost <- 100L
# Rho <- 0.5; Linear <- TRUE; Hetero <- TRUE
# 
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# ComputeKLD <- function(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                               n.pop,n.rct,n.obs,N,
#                                               Linear,Hetero,
#                                               TauLinear,Weight.same=NULL,
#                                               unitselect,pR){
#   #browser()
#   Outcome.0 <- Stratified.0.xa1 <- Stratified.0.xa0 <- Stratified.0.xb1 <- Stratified.0.xb0 <- Stratified.0.xc1 <- Stratified.0.xc0 <- array(NaN,c(12,N))
#   Outcome.xa <- Stratified.xa.xa1 <- Stratified.xa.xa0 <- Stratified.xa.xb1 <- Stratified.xa.xb0 <- Stratified.xa.xc1 <- Stratified.xa.xc0 <- array(NaN,c(12,N))
#   Outcome.xb <- Stratified.xb.xa1 <- Stratified.xb.xa0 <- Stratified.xb.xb1 <- Stratified.xb.xb0 <- Stratified.xb.xc1 <- Stratified.xb.xc0 <- array(NaN,c(12,N))
#   Outcome.xc <- Stratified.xc.xa1 <- Stratified.xc.xa0 <- Stratified.xc.xb1 <- Stratified.xc.xb0 <- Stratified.xc.xc1 <- Stratified.xc.xc0 <- array(NaN,c(12,N))
#   Outcome.xab <- Stratified.xab.xa1 <- Stratified.xab.xa0 <- Stratified.xab.xb1 <- Stratified.xab.xb0 <- Stratified.xab.xc1 <- Stratified.xab.xc0 <- array(NaN,c(length(seq(-1,1,0.1)),12,N))
#   Outcome.xac <- Stratified.xac.xa1 <- Stratified.xac.xa0 <- Stratified.xac.xb1 <- Stratified.xac.xb0 <- Stratified.xac.xc1 <- Stratified.xac.xc0 <- array(NaN,c(length(seq(-1,1,0.1)),12,N))
#   Outcome.xbc <- Stratified.xbc.xa1 <- Stratified.xbc.xa0 <- Stratified.xbc.xb1 <- Stratified.xbc.xb0 <- Stratified.xbc.xc1 <- Stratified.xbc.xc0 <- array(NaN,c(length(seq(-1,1,0.1)),12,N))
#   Outcome.xabc <-Stratified.xabc.xa1 <- Stratified.xabc.xa0 <- Stratified.xabc.xb1 <- Stratified.xabc.xb0 <- Stratified.xabc.xc1 <- Stratified.xabc.xc0 <- array(NaN,c(16,12,N))
#   Outcome.xabc.Rho <-Stratified.xabc.Rho.xa1 <- Stratified.xabc.Rho.xa0 <- Stratified.xabc.Rho.xb1 <- Stratified.xabc.Rho.xb0 <- Stratified.xabc.Rho.xc1 <- Stratified.xabc.Rho.xc0 <- array(NaN,c(12,N))
#   ATE.true <- data.frame(population=numeric(),rct.tau=numeric(),rct.obs=numeric(),obs.tau=numeric(),obs.obs=numeric())
#   Diff.dis <- data.frame(x0=numeric(),xa=numeric(),xb=numeric(),xc=numeric(),xab=numeric(),xac=numeric(),xbc=numeric(),xabc=numeric(),xabc.rho=numeric())
#   
#   Methods.Est <- array(NaN,c(N,11))
#   
#   n <- 0
#   while (n < N) {
#     n <- n+1
#     
#     # step 1: DGP
#     population <- DGP(Rho,n.pop,p.a.succ,p.b.succ,p.c.succ,Linear,Hetero,TauLinear,Weight.same)
#     
#     # step 2: sampling R=1/0
#     dat.rct <- population[sample(1:n.pop,size = n.rct),]
#     prob <- ProbUnitSelection(population,unitselect=unitselect,pR=pR)
#     dat.obs <- population[sample(1:dim(population)[1],size = n.obs,prob = prob),]
# 
#     pa <- (table(dat.rct[,1])/dim(dat.rct)[1])[2] #success probablity of xA in rct
#     pb <- (table(dat.rct[,2])/dim(dat.rct)[1])[2] #success probablity of xB in rct
#     pc <- (table(dat.rct[,3])/dim(dat.rct)[1])[2] #success probablity of xC in rct
#     
#     # ================== 
#     # step 1: no subsampling 
#     Diff.dis[n,'x0'] <- DKL(dat.rct,dat.obs)
# 
#     # ==================
#     # subsampling based on one confounding, xA
#     n.success <- dim(dat.obs)[1]*pa
#     dat.0 <- dat.obs[dat.obs$xA==0,]; dat.0 <- dat.0[sample(1:dim(dat.0)[1],size = dim(dat.obs)[1] - n.success ,replace = TRUE),]
#     dat.1 <- dat.obs[dat.obs$xA==1,]; dat.1 <- dat.1[sample(1:dim(dat.1)[1],size = n.success ,replace = TRUE),]
#     dat.obs.sub <- rbind(dat.0,dat.1)
#     Diff.dis[n,'xa'] <-  DKL(dat.rct,dat.obs.sub)
# 
#     # ==================
#     # subsampling based on one confounding, xB
#     n.success <- dim(dat.obs)[1]*pb
#     dat.0 <- dat.obs[dat.obs$xB==0,]; dat.0 <- dat.0[sample(1:dim(dat.0)[1],size = dim(dat.obs)[1] - n.success ,replace = TRUE),]
#     dat.1 <- dat.obs[dat.obs$xB==1,]; dat.1 <- dat.1[sample(1:dim(dat.1)[1],size = n.success ,replace = TRUE),]
#     dat.obs.sub <- rbind(dat.0,dat.1)
#     Diff.dis[n,'xb'] <-  DKL(dat.rct,dat.obs.sub)
# 
#     # ==================
#     # subsampling based on one confounding, xC
#     n.success <- dim(dat.obs)[1]*pc
#     dat.0 <- dat.obs[dat.obs$xC==0,]; dat.0 <- dat.0[sample(1:dim(dat.0)[1],size = dim(dat.obs)[1] - n.success ,replace = TRUE),]
#     dat.1 <- dat.obs[dat.obs$xC==1,]; dat.1 <- dat.1[sample(1:dim(dat.1)[1],size = n.success ,replace = TRUE),]
#     dat.obs.sub <- rbind(dat.0,dat.1)
#     Diff.dis[n,'xc'] <-  DKL(dat.rct,dat.obs.sub)
#     
#     # ==================
#     # subsampling based on two confounding: xA,xB
#     mu <- rep(0, 2)
#     i <- 1
#     diff.dis.mean <- NULL
#     for (rho in seq(-1,1,0.1)){
#       Sigma <- matrix(c(1,rho,
#                         rho, 1), nrow = 2, ncol = 2)
#       
#       rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
#       pvars <- pnorm(rawvars)
#       xA <- qbinom(pvars[,1],1,pa)
#       xB <- qbinom(pvars[,2],1,pb)
#       dat.00 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0),]; dat.00 <- dat.00[sample(1:dim(dat.00)[1],size = sum((xA==0)&(xB==0)),replace = TRUE),]
#       dat.10 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0),]; dat.10 <- dat.10[sample(1:dim(dat.10)[1],size = sum((xA==1)&(xB==0)),replace = TRUE),]
#       dat.01 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1),]; dat.01 <- dat.01[sample(1:dim(dat.01)[1],size = sum((xA==0)&(xB==1)),replace = TRUE),]
#       dat.11 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1),]; dat.11 <- dat.11[sample(1:dim(dat.11)[1],size = sum((xA==1)&(xB==1)),replace = TRUE),]
#       dat.obs.sub <- rbind(dat.00,dat.01,dat.10,dat.11)
#       diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
#       i <- i+1
#     }
#     Diff.dis[n,'xab'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
#     # ==================
#     # subsampling based on two confounding: xA,xC
#     mu <- rep(0, 2)
#     i <- 1
#     diff.dis.mean <- NULL
#     for (rho in seq(-1,1,0.1)){
#       Sigma <- matrix(c(1,rho,
#                         rho, 1), nrow = 2, ncol = 2)
#       
#       rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
#       pvars <- pnorm(rawvars)
#       xA <- qbinom(pvars[,1],1,pa)
#       xC <- qbinom(pvars[,2],1,pc)
#       dat.00 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xC==0),]; dat.00 <- dat.00[sample(1:dim(dat.00)[1],size = sum((xA==0)&(xC==0)),replace = TRUE),]
#       dat.10 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xC==0),]; dat.10 <- dat.10[sample(1:dim(dat.10)[1],size = sum((xA==1)&(xC==0)),replace = TRUE),]
#       dat.01 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xC==1),]; dat.01 <- dat.01[sample(1:dim(dat.01)[1],size = sum((xA==0)&(xC==1)),replace = TRUE),]
#       dat.11 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xC==1),]; dat.11 <- dat.11[sample(1:dim(dat.11)[1],size = sum((xA==1)&(xC==1)),replace = TRUE),]
#       dat.obs.sub <- rbind(dat.00,dat.01,dat.10,dat.11)
#       diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
#       i <- i+1
#     }
# 
#     # ==================
#     # subsampling based on two confounding: xB,xC
#     mu <- rep(0, 2)
#     i <- 1
#     diff.dis.mean <- NULL
#     for (rho in seq(-1,1,0.1)){
#       Sigma <- matrix(c(1,rho,
#                         rho, 1), nrow = 2, ncol = 2)
#       
#       rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
#       pvars <- pnorm(rawvars)
#       xB <- qbinom(pvars[,1],1,pb)
#       xC <- qbinom(pvars[,2],1,pc)
#       dat.00 <- dat.obs[(dat.obs$xB==0)&(dat.obs$xC==0),]; dat.00 <- dat.00[sample(1:dim(dat.00)[1],size = sum((xB==0)&(xC==0)),replace = TRUE),]
#       dat.10 <- dat.obs[(dat.obs$xB==1)&(dat.obs$xC==0),]; dat.10 <- dat.10[sample(1:dim(dat.10)[1],size = sum((xB==1)&(xC==0)),replace = TRUE),]
#       dat.01 <- dat.obs[(dat.obs$xB==0)&(dat.obs$xC==1),]; dat.01 <- dat.01[sample(1:dim(dat.01)[1],size = sum((xB==0)&(xC==1)),replace = TRUE),]
#       dat.11 <- dat.obs[(dat.obs$xB==1)&(dat.obs$xC==1),]; dat.11 <- dat.11[sample(1:dim(dat.11)[1],size = sum((xB==1)&(xC==1)),replace = TRUE),]
#       dat.obs.sub <- rbind(dat.00,dat.01,dat.10,dat.11)
#       diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
#       i <- i+1
#     }
#     Diff.dis[n,'xbc'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
#     
#     # ==================
#     # subsampling based on three confounding: xA,xB,xC
#     mu <- rep(0, 3)
#     i <- 1
#     diff.dis.mean <- NULL
#     for (rho in seq(-1,1,0.1)){
#       try(
#         {
#           Sigma <- matrix(c(1,rho,rho,
#                             rho, 1,rho,
#                             rho,rho,1), nrow = 3, ncol = 3)
#           
#           rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
#           pvars <- pnorm(rawvars)
#           xA <- qbinom(pvars[,1],1,pa)
#           xB <- qbinom(pvars[,2],1,pb)
#           xC <- qbinom(pvars[,3],1,pc)
#           dat.000 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
#           dat.000 <- dat.000[sample(1:dim(dat.000)[1],size = sum((xA==0)&(xB==0)&(xC==0)),replace = TRUE),]
#           dat.001 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
#           dat.001 <- dat.001[sample(1:dim(dat.001)[1],size = sum((xA==0)&(xB==0)&(xC==1)),replace = TRUE),]
#           dat.010 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
#           dat.010 <- dat.010[sample(1:dim(dat.010)[1],size = sum((xA==0)&(xB==1)&(xC==0)),replace = TRUE),]
#           dat.011 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
#           dat.011 <- dat.011[sample(1:dim(dat.011)[1],size = sum((xA==0)&(xB==1)&(xC==1)),replace = TRUE),]
#           dat.100 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
#           dat.100 <- dat.100[sample(1:dim(dat.100)[1],size = sum((xA==1)&(xB==0)&(xC==0)),replace = TRUE),]
#           dat.101 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
#           dat.101 <- dat.101[sample(1:dim(dat.101)[1],size = sum((xA==1)&(xB==0)&(xC==1)),replace = TRUE),]
#           dat.110 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
#           dat.110 <- dat.110[sample(1:dim(dat.110)[1],size = sum((xA==1)&(xB==1)&(xC==0)),replace = TRUE),]
#           dat.111 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
#           dat.111 <- dat.111[sample(1:dim(dat.111)[1],size = sum((xA==1)&(xB==1)&(xC==1)),replace = TRUE),]
#           
#           dat.obs.sub <- rbind(dat.000,dat.001,dat.010,dat.011,dat.100,dat.101,dat.110,dat.111)
#           diff.dis.mean <- c(diff.dis.mean,DKL(dat.rct,dat.obs.sub))
#           i <- i+1
#         },
#         silent = TRUE
#       )
#       
#     }
#     Diff.dis[n,'xabc'] <- mean(diff.dis.mean[!is.infinite(diff.dis.mean)])
#     
#     
#     # ==================
#     # subsampling based on three confounding and correct dependency Rho: xA,xB,xC
#     rhoAB <- simstudy:::.findRhoBin(p1 = pa, 
#                                     p2 = pb, d = cor(dat.rct$xA,dat.rct$xB))
#     rhoAC <- simstudy:::.findRhoBin(p1 = pa, 
#                                     p2 = pc, d = cor(dat.rct$xA,dat.rct$xC))
#     rhoBC <- simstudy:::.findRhoBin(p1 = pb, 
#                                     p2 = pc, d = cor(dat.rct$xB,dat.rct$xC))
#     
#     mu <- rep(0, 3)
#     
#     Sigma <- matrix(c(1,rhoAB,rhoAC,
#                       rhoAB,1,rhoBC,
#                       rhoAC,rhoBC,1), nrow = 3, ncol = 3)
#     
#     rawvars <- mvrnorm(n = n.obs, mu = mu, Sigma = Sigma)
#     pvars <- pnorm(rawvars)
#     xA <- 1*(pvars[,1] < pa)
#     xB <- 1*(pvars[,2] < pb)
#     xC <- 1*(pvars[,3] < pc)
#     dat.000 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
#     dat.000 <- dat.000[sample(1:dim(dat.000)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==0)&(dat.rct$xC==0)),replace = TRUE),]
#     dat.001 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
#     dat.001 <- dat.001[sample(1:dim(dat.001)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==0)&(dat.rct$xC==1)),replace = TRUE),]
#     dat.010 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
#     dat.010 <- dat.010[sample(1:dim(dat.010)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==1)&(dat.rct$xC==0)),replace = TRUE),]
#     dat.011 <- dat.obs[(dat.obs$xA==0)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
#     dat.011 <- dat.011[sample(1:dim(dat.011)[1],size = sum((dat.rct$xA==0)&(dat.rct$xB==1)&(dat.rct$xC==1)),replace = TRUE),]
#     dat.100 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==0),]; 
#     dat.100 <- dat.100[sample(1:dim(dat.100)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==0)&(dat.rct$xC==0)),replace = TRUE),]
#     dat.101 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==0)&(dat.obs$xC==1),]; 
#     dat.101 <- dat.101[sample(1:dim(dat.101)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==0)&(dat.rct$xC==1)),replace = TRUE),]
#     dat.110 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==0),]; 
#     dat.110 <- dat.110[sample(1:dim(dat.110)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==1)&(dat.rct$xC==0)),replace = TRUE),]
#     dat.111 <- dat.obs[(dat.obs$xA==1)&(dat.obs$xB==1)&(dat.obs$xC==1),]; 
#     dat.111 <- dat.111[sample(1:dim(dat.111)[1],size = sum((dat.rct$xA==1)&(dat.rct$xB==1)&(dat.rct$xC==1)),replace = TRUE),]
#     
#     dat.obs.sub <- rbind(dat.000,dat.001,dat.010,dat.011,dat.100,dat.101,dat.110,dat.111)
#     Diff.dis[n,'xabc.rho'] <- DKL(dat.rct,dat.obs.sub)
#   }
#  
#   return(Diff.dis)
# }



# # check only 16 rho can generate positive semi-definite covariance matrix
# i <- 0
# for (rho12 in seq(-1,1,0.1)) {
#   for (rho13 in seq(-1,1,0.1)) {
#     for (rho23 in seq(-1,1,0.1)) {
#       try({
#         mu <- rep(0, 3)
#         Sigma <- matrix(c(1,rho12,rho13,
#                           rho12,1,rho23,
#                           rho13,rho23,1), nrow = 3, ncol = 3)
#         rawvars <- mvrnorm(n = 1e3, mu = mu, Sigma = Sigma)
#         print(c(rho12,rho13,rho23))
#         i <- i+1
#       },
#       silent = TRUE)
#     }
#   }
# }
# 

#So the distance of estimate to rct lies in two aspects:
#  1) if variables are used to increase or decrese the probability of sampling RCT units (prob=ifelse(xA==1,prob,1-prob)), here xA vs. no xA and prob informative vs. non informative matters.
#  2) if variables are correlated with tau, and the weight of variables. 
#     prob=ifelse(xA==1,prob,1-prob); tau <- xB+xC vs. tau<- xA+xB+xC
#     prob=ifelse(xAxBxC==1,prob,1-prob); tau <- xA+xB+xC vs. tau<- 3*xA+xB+xC

