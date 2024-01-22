setwd("~/Documents/iknl/code/simulation_1")
#rm(list = ls())
library(glmnet)
library(BART)
library(MASS)
library(boot)
library(prodlim) # for function ComputeDistance2RctAndMethodsBiasWeighting()
source('function.R')
#So the distance of estimate to rct lies in two aspects:
#  1) if variables are used to increase or decrese the probability of sampling RCT units (prob=ifelse(xA==1,prob,1-prob)), here xA vs. no xA and prob informative vs. non informative matters.
#  2) if variables are correlated with tau, and the weight of variables. 
#     prob=ifelse(xA==1,prob,1-prob); tau <- xB+xC vs. tau<- xA+xB+xC
#     prob=ifelse(xAxBxC==1,prob,1-prob); tau <- xA+xB+xC vs. tau<- 3*xA+xB+xC

# ==================== general parameters
# p.a.succ <- p.b.succ <- p.c.succ <- 0.5
# n.pop <- 1e4
# n.rct <- 3000
# n.obs <- 7000
# lambda <- -log(3)
# N <- 50 # the number of times for DGP
# ntree <- 50
# ndpost <- 100L
# 
# # ====================
# # ====================
# # scenarios 1.1: mu linear + tau heterogenous + Rho 0.5
# Rho <- 0.5; Linear <- TRUE; Hetero <- TRUE
# 
# # DGM: 1.1.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.1 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 1.1.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.2 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 1.1.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.3 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 1.1.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.4 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# 
# 
# 
# # ====================
# # scenarios 1.2: mu linear + tau homogeneous + Rho 0.5
# Rho <- 0.5; Linear <- TRUE; Hetero <- FALSE
# 
# # DGM: 1.2.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.1 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 1.2.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.2 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 1.2.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.3 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 1.2.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.4 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# 
# 
# # ====================
# # scenarios 2.1: mu non-linear + tau heterogenous + Rho 0.5
# Rho <- 0.5; Linear <- FALSE; Hetero <- TRUE
# 
# # DGM: 2.1.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.1 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 2.1.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.2 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 2.1.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.3 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 2.1.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.4 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# 
# 
# # ====================
# # scenarios 2.2:  mu non-linear + tau homogenous + Rho 0.5
# Rho <- 0.5; Linear <- FALSE; Hetero <- FALSE
# # DGM: 2.2.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.1 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# 
# # DGM: 2.2.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.2 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 2.2.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.3 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# # DGM: 2.2.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.4 <- ComputeDistance2RctAndMethodsBiasWeighting(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                     n.pop,n.rct,n.obs,N,
#                                                     lambda,
#                                                     Linear,Hetero,
#                                                     TauLinear,Weight.same,
#                                                     unitselect,pR,
#                                                     ntree,ndpost)
# 
# save(Estimation.1.1,Estimation.1.2,Estimation.1.3,Estimation.1.4,
#      Estimation.2.1,Estimation.2.2,Estimation.2.3,Estimation.2.4,
#      Estimation.3.1,Estimation.3.2,Estimation.3.3,Estimation.3.4,
#      Estimation.4.1,Estimation.4.2,Estimation.4.3,Estimation.4.4,
#      file = 'results/simulation_weighting.RData')




# ======================================================
# ======================================================
# PRINT RESULTS 
# ======================================================
# ======================================================
load('results/simulation_weighting.RData')


# results of distance d to estimated RCT: DGM rbind(DGM: 1.1,1.2,2.1,2.2)
# rbind(Estimation.1.1[[2]],Estimation.1.2[[2]],Estimation.1.3[[2]],Estimation.1.4[[2]]) 
# rbind(Estimation.2.1[[2]],Estimation.2.2[[2]],Estimation.2.3[[2]],Estimation.2.4[[2]])
# rbind(Estimation.3.1[[2]],Estimation.3.2[[2]],Estimation.3.3[[2]],Estimation.3.4[[2]])
# rbind(Estimation.4.1[[2]],Estimation.4.2[[2]],Estimation.4.3[[2]],Estimation.4.4[[2]])
# 
# # results of distance d to the true RCT: DGM rbind(DGM: 1.1,1.2,2.1,2.2)
# rbind(Estimation.1.1[[3]],Estimation.1.2[[3]],Estimation.1.3[[3]],Estimation.1.4[[3]])
# rbind(Estimation.2.1[[3]],Estimation.2.2[[3]],Estimation.2.3[[3]],Estimation.2.4[[3]])
# rbind(Estimation.3.1[[3]],Estimation.3.2[[3]],Estimation.3.3[[3]],Estimation.3.4[[3]])
# rbind(Estimation.4.1[[3]],Estimation.4.2[[3]],Estimation.4.3[[3]],Estimation.4.4[[3]])
# 
# # results of KL divergence
# KL.div <- data.frame(rbind(t(cbind(round(apply(Estimation.1.1[[4]],2,mean),digits = 3),round(apply(Estimation.2.1[[4]],2,mean),digits = 3),
#                                    round(apply(Estimation.3.1[[4]],2,mean),digits = 3),round(apply(Estimation.4.1[[4]],2,mean),digits = 3))),
#                            t(cbind(round(apply(Estimation.1.2[[4]],2,mean),digits = 3),round(apply(Estimation.2.2[[4]],2,mean),digits = 3),
#                                    round(apply(Estimation.3.2[[4]],2,mean),digits = 3),round(apply(Estimation.4.2[[4]],2,mean),digits = 3))),
#                            t(cbind(round(apply(Estimation.1.3[[4]],2,mean),digits = 3),round(apply(Estimation.2.3[[4]],2,mean),digits = 3),
#                                    round(apply(Estimation.3.3[[4]],2,mean),digits = 3),round(apply(Estimation.4.3[[4]],2,mean),digits = 3))),
#                            t(cbind(round(apply(Estimation.1.4[[4]],2,mean),digits = 3),round(apply(Estimation.2.4[[4]],2,mean),digits = 3),
#                                    round(apply(Estimation.3.4[[4]],2,mean),digits = 3),round(apply(Estimation.4.4[[4]],2,mean),digits = 3)))
#                            )
#                       )
# KL.DIV <- NULL
# Colname <- NULL
# for(i in seq(dim(KL.div)[2])){
#   kl.div <- KL.div[,i]
#   KL.DIV <- cbind(KL.DIV,kl.div,'&')
#   colname.i <- colnames(KL.div)[i]
#   Colname <- cbind(Colname,colname.i,' ')
# }  
# KL.DIV <- as.data.frame(KL.DIV)
# colnames(KL.DIV) <- Colname
# print(KL.DIV,quote=FALSE,row.names=FALSE)
# 
# 
# # plot the simulation results of d between weighted estimated of obs and estimated of rct
# library(ggpubr)
# # height=8, width=12
# FourSamplingPlot(Estimation.1.1[[2]],Estimation.1.2[[2]],Estimation.1.3[[2]],Estimation.1.4[[2]]) #scenario 1: Rho <- 0.5; Linear <- TRUE; Hetero <- TRUE
# FourSamplingPlot(Estimation.2.1[[2]],Estimation.2.2[[2]],Estimation.2.3[[2]],Estimation.2.4[[2]]) #scenario 2: Rho <- 0.5; Linear <- TRUE; Hetero <- FALSE
# FourSamplingPlot(Estimation.3.1[[2]],Estimation.3.2[[2]],Estimation.3.3[[2]],Estimation.3.4[[2]]) #scenario 3: Rho <- 0.5; Linear <- FALSE; Hetero <- TRUE
# FourSamplingPlot(Estimation.4.1[[2]],Estimation.4.2[[2]],Estimation.4.3[[2]],Estimation.4.4[[2]]) #scenario 4: Rho <- 0.5; Linear <- FALSE; Hetero <- FALSE

# performance measurement for methods of 4 scenarios
# DGM 1.1
MethodsEst11 <- Estimation.1.1[[1]];MethodsEst12 <- Estimation.1.2[[1]];
MethodsEst13 <- Estimation.1.3[[1]];MethodsEst14 <- Estimation.1.4[[1]];
print(rbind(BiasEmpSeMSE(MethodsEst11),
            BiasEmpSeMSE(MethodsEst12),
            BiasEmpSeMSE(MethodsEst13),
            BiasEmpSeMSE(MethodsEst14)),
      quote=FALSE)

# DGM 1.2
MethodsEst21 <- Estimation.2.1[[1]];MethodsEst22 <- Estimation.2.2[[1]];
MethodsEst23 <- Estimation.2.3[[1]];MethodsEst24 <- Estimation.2.4[[1]];
print(rbind(BiasEmpSeMSE(MethodsEst21),
            BiasEmpSeMSE(MethodsEst22),
            BiasEmpSeMSE(MethodsEst23),
            BiasEmpSeMSE(MethodsEst24)),
      quote=FALSE)

# DGM 2.1
MethodsEst31 <- Estimation.3.1[[1]];MethodsEst32 <- Estimation.3.2[[1]];
MethodsEst33 <- Estimation.3.3[[1]];MethodsEst34 <- Estimation.3.4[[1]];
print(rbind(BiasEmpSeMSE(MethodsEst31),
            BiasEmpSeMSE(MethodsEst32),
            BiasEmpSeMSE(MethodsEst33),
            BiasEmpSeMSE(MethodsEst34)),
      quote=FALSE)

# DGM 2.2
MethodsEst41 <- Estimation.4.1[[1]];MethodsEst42 <- Estimation.4.2[[1]];
MethodsEst43 <- Estimation.4.3[[1]];MethodsEst44 <- Estimation.4.4[[1]];
print(rbind(BiasEmpSeMSE(MethodsEst41),
            BiasEmpSeMSE(MethodsEst42),
            BiasEmpSeMSE(MethodsEst43),
            BiasEmpSeMSE(MethodsEst44)),
      quote=FALSE
      )
