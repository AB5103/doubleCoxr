library(doubleCox)
library(ucminf)
library(MASS)
library(xtable)
library(survival)
data("kidney",package="survival")
kidney$sex <- kidney$sex - 1
formula.scale=as.formula('Surv(time, status) ~ age + sex')
formula.shape=as.formula('Surv(time, status) ~ 1')
cluster='id'
dist='Weibull'
ResultTest<-ParNPHCox(formula.scale,formula.shape,cluster,dist,data=kidney,NULL,NULL)
