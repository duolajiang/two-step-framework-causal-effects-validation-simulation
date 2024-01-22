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
# Estimation.1.1 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                        n.pop,n.rct,n.obs,N,
#                                                        lambda,
#                                                        Linear,Hetero,
#                                                        TauLinear,Weight.same,
#                                                        unitselect,pR,
#                                                        ntree,ndpost)
# # DGM: 1.1.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.2 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                        n.pop,n.rct,n.obs,N,
#                                                        lambda,
#                                                        Linear,Hetero,
#                                                        TauLinear,Weight.same,
#                                                        unitselect,pR,
#                                                        ntree,ndpost)
# # DGM: 1.1.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.3 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 1.1.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.1.4 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# 
# 
# 
# # ====================
# # scenarios 1.2: mu linear + tau homogeneous + Rho 0.5
# Rho <- 0.5; Linear <- TRUE; Hetero <- FALSE
# 
# # DGM: 1.2.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.1 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 1.2.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.2 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 1.2.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.3 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 1.2.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.2.4 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# 
# 
# # ====================
# # scenarios 2.1: mu non-linear + tau heterogenous + Rho 0.5
# Rho <- 0.5; Linear <- FALSE; Hetero <- TRUE
# 
# # DGM: 2.1.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.1 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 2.1.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.2 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 2.1.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.3 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 2.1.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.3.4 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# 
# 
# # ====================
# # scenarios 2.2:  mu non-linear + tau homogenous + Rho 0.5
# Rho <- 0.5; Linear <- FALSE; Hetero <- FALSE
# # DGM: 2.2.1
# set.seed(961); unitselect <- 'xA'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.1 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# 
# # DGM: 2.2.2
# set.seed(961); unitselect <- 'xB'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.2 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 2.2.3
# set.seed(961); unitselect <- 'xC'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.3 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# # DGM: 2.2.4
# set.seed(961); unitselect <- 'all'; pR <- 0.8; TauLinear <- TRUE; Weight.same <- TRUE
# Estimation.4.4 <- ComputeWeightedEstimatesBiasVariance(Rho,p.a.succ,p.b.succ,p.c.succ,
#                                                              n.pop,n.rct,n.obs,N,
#                                                              lambda,
#                                                              Linear,Hetero,
#                                                              TauLinear,Weight.same,
#                                                              unitselect,pR,
#                                                              ntree,ndpost)
# 
# save(Estimation.1.1,file = 'results/simulation_weighting_11.RData')
# save(Estimation.1.2,file = 'results/simulation_weighting_12.RData')
# save(Estimation.1.3,file = 'results/simulation_weighting_13.RData')
# save(Estimation.1.4,file = 'results/simulation_weighting_14.RData')
# save(Estimation.2.1,file = 'results/simulation_weighting_21.RData')
# save(Estimation.2.2,file = 'results/simulation_weighting_22.RData')
# save(Estimation.2.3,file = 'results/simulation_weighting_23.RData')
# save(Estimation.2.4,file = 'results/simulation_weighting_24.RData')
# save(Estimation.3.1,file = 'results/simulation_weighting_31.RData')
# save(Estimation.3.2,file = 'results/simulation_weighting_32.RData')
# save(Estimation.3.3,file = 'results/simulation_weighting_33.RData')
# save(Estimation.3.4,file = 'results/simulation_weighting_34.RData')
# save(Estimation.4.1,file = 'results/simulation_weighting_41.RData')
# save(Estimation.4.2,file = 'results/simulation_weighting_42.RData')
# save(Estimation.4.3,file = 'results/simulation_weighting_43.RData')
# save(Estimation.4.4,file = 'results/simulation_weighting_44.RData')

# =============================================
# generate plots 
# =============================================
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# load("results/simulation_weighting_11.RData")
# load("results/simulation_weighting_12.RData")
# load("results/simulation_weighting_13.RData")
# load("results/simulation_weighting_14.RData")
# ## note only p1 can be in between pdf() and dev.off(), otherwise the first page is blank
# p1 <- FourSamplingBoxPlot(Estimation.1.1,Estimation.1.2,Estimation.1.3,Estimation.1.4)
# pdf(file = "DGM_linear_hetero.pdf",width = 8,height = 6.5)
# p1
# dev.off()
# 
# load("results/simulation_weighting_21.RData")
# load("results/simulation_weighting_22.RData")
# load("results/simulation_weighting_23.RData")
# load("results/simulation_weighting_24.RData")
# p2 <- FourSamplingBoxPlot (Estimation.2.1,Estimation.2.2,Estimation.2.3,Estimation.2.4)
# pdf(file = "DGM_linear_home.pdf",width = 8,height = 6.5)
# p2
# dev.off()
# 
# load("results/simulation_weighting_31.RData")
# load("results/simulation_weighting_32.RData")
# load("results/simulation_weighting_33.RData")
# load("results/simulation_weighting_34.RData")
# p3 <- FourSamplingBoxPlot (Estimation.3.1,Estimation.3.2,Estimation.3.3,Estimation.3.4)
# pdf(file = "DGM_nonlinear_hetero.pdf",width = 8,height = 6.5)
# p3
# dev.off()
# 
# load("results/simulation_weighting_41.RData")
# load("results/simulation_weighting_42.RData")
# load("results/simulation_weighting_43.RData")
# load("results/simulation_weighting_44.RData")
# p4 <- FourSamplingBoxPlot (Estimation.4.1,Estimation.4.2,Estimation.4.3,Estimation.4.4)
# pdf(file = "DGM_nonlinear_homo.pdf",width = 8,height = 6.5)
# p4
# dev.off()

# =============================================
# plot main figure for the paper
# =============================================
load("results/simulation_weighting_14.RData")
load("results/simulation_weighting_24.RData")
load("results/simulation_weighting_34.RData")
load("results/simulation_weighting_44.RData")

p <- FourSamplingBoxPlot_simple(Estimation.1.4,Estimation.2.4,Estimation.3.4,Estimation.4.4)
pdf(file = "plots/4DGMs_pi_X.pdf",width = 7,height = 6.5)
p
dev.off()

# =============================================
# plot secondary plot for the paper
# =============================================
load("results/simulation_weighting_11.RData")
load("results/simulation_weighting_12.RData")
load("results/simulation_weighting_13.RData")
load("results/simulation_weighting_14.RData")
## note only p1 can be in between pdf() and dev.off(), otherwise the first page is blank
p <- FourSamplingBoxPlot_no_color(Estimation.1.1,Estimation.1.2,Estimation.1.3,Estimation.1.4)
pdf(file = "plots/DGM_linear_hetero_psBART_4S.pdf",width = 8,height = 6.5)
p
dev.off()

