#' @title GrGenNPH
#'
#' @description   This function calculates the neg-gradient of the loglikelihood
#' for the parametric gamma-frailty model with non-proportional hazard functions
#'
#' @param y Vector of parameters in the form
#' \deqn{y = (ln(a), ln(b), \beta _{shape}, \beta _{scale}, ln(\sigma ^2))}
#' for Weibull hazard function and
#' \deqn{y = (ln(10^3a), ln(10^2b), \beta _{shape}, \beta _{scale}, ln(\sigma ^2))}
#' for Gompertz hazard function, where a and b are slope and shape parameters,
#' \eqn{\beta _{shape}} and \eqn{\beta _{scale}} are the Cox-regression parameters
#'  for shape and scale, respectively, and \eqn{\sigma ^2} is the variance of frailty.
#'  This vector must include at least two parameters, \eqn{ln(a)} and \eqn{ln(b)}.
#' @param D  A data.frame in which to interpret the variables named in the formula.
#' The data set includes the following fields:
#' \enumerate{
#' \item time-to-failure and censoring in the case without left truncation
#' or time-of-start, time-of-failure, and censoring in the case with left truncation at the time of begin
#' (censoring must be either 0 for no event or 1 for event);
#' \item Covariates (continuous or categorical) used in a study (can be empty set).
#' }
#' @param nf The number of continuous and binary factors in the data set D corresponding to the covariates
#' used in the Cox-regression for proportional hazard term.
#' @param nk  The number of continuous and binary factors in the data set D corresponding
#' to the covariates used in the Cox-regression for shape b.
#' @param ncl The number of clusters in the data set D corresponding to the cluster covariate.
#' Is equal to 0 for the fixed-effect model.
#' @param dist Baseline hazard function ('Weibull' or 'Gompertz').
#'
#' @return Neg-gradient of the loglikelihood
#'
#' @examples
#' \dontrun{
#' GrGenNPH(y, D, nf, nk, ncl, dist)
#' }
#'
#' @export
#'
GrGenNPH=function(y,D,nf,nk,ncl,dist){
  b=y
  if (dist=='Weibull'){
    lambda0=exp(b[1])
    k0=exp(b[2])} else if (dist=='Gompertz') {
      lambda0=1e-3*exp(b[1])
      k0=1e-2*exp(b[2])
    }
  Coxk=c(rep(0,nrow(D)))
  Cox=c(rep(0,nrow(D)))
  if (nk>0){
    for (i in 1:nk){
      Coxk=Coxk+D[,i]*b[2+i]
    }
  }
  if (nf>0){
    for (i in 1:nf){
      Cox=Cox+D[,(nk+i)]*b[2+nk+i]
    }
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
  if (ncl>0) {
    G2=exp(b[3+nk+nf])
  }
  list=unique(D$cluster)
  nl=length(list)

  k0=k0*Coxk
  Cens=D$event
  x1=D$time
  x0=D$trunc
  Gr=c(rep(0,length(b)))
  if (dist=='Weibull'){
    for (i in 1:nl){
      ID=list[i]
      ind=which(D$cluster==ID)
      nn=length(ind)
      ic1=1*(Cens[ind]==1)
      nc1=sum(ic1)
      Hfull1= Cox[ind]*(x1[ind]/lambda0)^k0[ind]
      Hfull0= Cox[ind]*(x0[ind]/lambda0)^k0[ind]
      H_lambda_1=-Hfull1*k0[ind]
      lmu_lambda_1=-k0[ind]
      H_lambda_0=-Hfull0*k0[ind]
      H_k0_1=Hfull1*(log(x1[ind])-log(lambda0))*k0[ind]
      lmu_k0_1=k0[ind]*(log(x1[ind])-log(lambda0))+1
      H_k0_0=c(rep(0,nn))
      H_k0_0[x0[ind]>0]=Hfull0[x0[ind]>0]*(log(x0[ind][x0[ind]>0])-log(lambda0))*k0[ind][x0[ind]>0]
      H_beta_1={}
      if (nf>0)  H_beta_1=matrix(rep(Hfull1,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
      H_beta_0={}
      if (nf>0)  H_beta_0=matrix(rep(Hfull0,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
      H_betak_1={}
      if (nk>0)  {
        H_betak_1=matrix(rep(Hfull1*(log(x1[ind])-log(lambda0))*k0[ind],nk),nn,nk)*D[ind,1:nk]
      }
      H_betak_0={}
      if (nk>0)  {
        H_betak_0=matrix(0,nn,nk)
        if (sum(x0[ind]>0)>0) H_betak_0[x0[ind]>0,]=matrix(rep(Hfull0*(log(x1[ind])-log(lambda0))*k0[ind],nk),nn,nk)[x0[ind]>0,]*D[ind,1:nk][x0[ind]>0,]
      }
      lmu_beta_1=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
      lmu_beta_0=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
      lmu_betak_1={}
      if (nk>0)  {
        lmu_betak_1=matrix(rep(ic1*(k0[ind]*(log(x1[ind])-log(lambda0))+1),nk),nn,nk)*D[ind,1:nk]
      }
      if (ncl>0){
        Gr[1]=Gr[1]-(1+G2*nc1)*sum(H_lambda_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)/(1+G2*sum(Hfull0))
        Gr[2]=Gr[2]-(1+G2*nc1)*sum(H_k0_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_k0_1)+sum(H_k0_0)/(1+G2*sum(Hfull0))
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-(1+G2*nc1)*c(apply(H_betak_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))/(1+G2*sum(Hfull0))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-(1+G2*nc1)*c(apply(H_beta_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(H_beta_0,2,sum))/(1+G2*sum(Hfull0))+c(apply(lmu_beta_1,2,sum))
        Gr[length(b)]=Gr[length(b)]+log(1+G2*sum(Hfull1))/G2-(1+G2*nc1)*sum(Hfull1)/(1+G2*sum(Hfull1))+(digamma(1/G2)-digamma(nc1+1/G2))/G2-log(1+G2*sum(Hfull0))/G2-sum(Hfull0)/(1+G2*sum(Hfull0))+nc1
      } else {
        Gr[1]=Gr[1]-sum(H_lambda_1)+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)
        Gr[2]=Gr[2]-sum(H_k0_1)+sum(ic1*lmu_k0_1)+sum(H_k0_0)
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-c(apply(H_betak_1,2,sum))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-c(apply(H_beta_1,2,sum))+c(apply(lmu_beta_1,2,sum))+c(apply(H_beta_0,2,sum))
      }
    }
  }

  if (dist=='Gompertz'){
    for (i in 1:nl){
      ID=list[i]
      ind=which(D$cluster==ID)
      nn=length(ind)
      ic1=1*(Cens[ind]==1)
      nc1=sum(ic1)
      Hfull1= Cox[ind]*(lambda0/k0[ind])*(exp(k0[ind]*x1[ind])-1)
      Hfull0= Cox[ind]*(lambda0/k0[ind])*(exp(k0[ind]*x0[ind])-1)
      mufull1=Cox[ind]*lambda0*exp(k0[ind]*x1[ind])
      mufull0=Cox[ind]*lambda0*exp(k0[ind]*x0[ind])
      Lmufull1=LCox[ind]+log(lambda0)+k0[ind]*x1[ind]
      H_lambda_1=Hfull1
      H_lambda_0=Hfull0
      lmu_lambda_1=1
      H_k0_1=-Hfull1+mufull1*x1[ind]
      H_k0_0=-Hfull0+mufull0*x0[ind]
      lmu_k0_1=k0[ind]*x1[ind]
      H_beta_1={}
      if (nf>0)  H_beta_1=matrix(rep(Hfull1,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
      H_beta_0={}
      if (nf>0)  H_beta_0=matrix(rep(Hfull0,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
      H_betak_1={}
      if (nk>0)  {
        H_betak_1=matrix(rep(-Hfull1+mufull1*x1[ind],nk),nn,nk)*D[ind,1:nk]
      }
      H_betak_0={}
      if (nk>0)  {
        H_betak_0=matrix(rep(-Hfull0+mufull0*x0[ind],nk),nn,nk)*D[ind,1:nk]
      }
      lmu_beta_1=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
      lmu_beta_0=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
      lmu_betak_1={}
      if (nk>0)  {
        lmu_betak_1=matrix(rep(ic1*k0[ind]*x1[ind],nk),nn,nk)*D[ind,1:nk]
      }
      if (ncl>0){
        Gr[1]=Gr[1]-(1+G2*nc1)*sum(H_lambda_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)/(1+G2*sum(Hfull0))
        Gr[2]=Gr[2]-(1+G2*nc1)*sum(H_k0_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_k0_1)+sum(H_k0_0)/(1+G2*sum(Hfull0))
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-(1+G2*nc1)*c(apply(H_betak_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))/(1+G2*sum(Hfull0))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-(1+G2*nc1)*c(apply(H_beta_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(H_beta_0,2,sum))/(1+G2*sum(Hfull0))+c(apply(lmu_beta_1,2,sum))
        Gr[length(b)]=Gr[length(b)]+log(1+G2*sum(Hfull1))/G2-(1+G2*nc1)*sum(Hfull1)/(1+G2*sum(Hfull1))+(digamma(1/G2)-digamma(nc1+1/G2))/G2-log(1+G2*sum(Hfull0))/G2-sum(Hfull0)/(1+G2*sum(Hfull0))+nc1
      } else {
        Gr[1]=Gr[1]-sum(H_lambda_1)+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)
        Gr[2]=Gr[2]-sum(H_k0_1)+sum(ic1*lmu_k0_1)+sum(H_k0_0)
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-c(apply(H_betak_1,2,sum))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-c(apply(H_beta_1,2,sum))+c(apply(lmu_beta_1,2,sum))+c(apply(H_beta_0,2,sum))
      }
    }
  }
  return(-Gr)
}
