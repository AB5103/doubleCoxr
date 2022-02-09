#' @title cumhazard
#'
#' @description This function calculates cumulative hazard
#'
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
#' to the covariates used in the Cox-regression for shape \eqn{b}.
#' @param ncl The number of clusters in the data set D corresponding to the cluster covariate.
#' Is equal to 0 for the fixed-effect model).
#' @param dist Baseline hazard function ('Weibull' or 'Gompertz').
#' @param vpar The sample of the vector-row-parameters in matrix form.
#' @param ID The order number of an object in the data set.
#'
#' @return Cumulative hazard for the object with order number ID
#'
#' @examples
#' \dontrun{
#' cumhazard(D,nf,nk,ncl,dist,vpar,ID)
#' }
#'
#' @export
#'
cumhazard=function(D,nf,nk,ncl,dist,vpar,ID){
    if (dist=='Weibull'){
    lambda0=exp(vpar[,1])
    k0=exp(vpar[,2])} else if (dist=='Gompertz') {
      lambda0=1e-3*exp(vpar[,1])
      k0=1e-2*exp(vpar[,2])
    }
  Coxk=0
  Cox=0
  if (nk>0){
    for (i in 1:nk){
      Coxk=Coxk+D[ID,i]*vpar[,2+i]}
  }
  if (nf>0){
    for (i in 1:nf){
      Cox=Cox+D[ID,(nk+i)]*vpar[,2+nk+i]}
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
  if (ncl>0) {
    G2=exp(vpar[,3+nk+nf])
  }
  k0=k0*Coxk
  x1=D$time[ID]
  x0=D$trunc[ID]
  if (dist=='Weibull'){
    Hfull1= Cox*(x1/lambda0)^k0
    Hfull0= Cox*(x0/lambda0)^k0
    mufull1=Cox*k0*x1^(k0-1)/lambda0^k0
    cumh=Hfull1-Hfull0
  }

  if (dist=='Gompertz'){
    Hfull1= (Cox*lambda0/k0)*(exp(k0*x1)-1)
    Hfull0= (Cox*lambda0/k0)*(exp(k0*x0)-1)
    mufull1= (Cox*lambda0)*exp(k0*x1)
    cumh=Hfull1-Hfull0
  }
  return(cumh)
}
