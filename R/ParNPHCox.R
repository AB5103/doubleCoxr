#' @title ParNPHCox
#'
#' @description   Provides parameter estimates, their standard errors, loglikelihood and other attributes
#' for the parametric gamma-frailty model with non-proportional hazard functions
#'
#' @param formula.scale A formula object, with the response on the left of a ~ operator, and the terms on
#' the right. The response must be a survival object as returned by the Surv() function.
#' The status indicator 'event' in the Surv object must be 0=censored, 1=non-censored.
#' This object describes the effect of several factors in the Cox proportional hazards model.
#' @param formula.shape A formula object, with the response on the left of a ~ operator, and the terms on
#' the right. The response must be a survival object as returned by the Surv() function.
#' The status indicator 'event' in the Surv object must be 0=censored, 1=non-censored.
#' This object describes the effect of several factors on the shape parameter.
#' @param cluster The name of a cluster variable in data (is equal to NULL for the fixed-effect model)
#' @param dist Baseline hazard function ('Weibull' or 'Gompertz').
#' @param data  A data.frame in which to interpret the variables named in the formula.
#' The data set includes the following fields:
#' \enumerate{
#' \item time-to-failure and censoring in the case without left truncation
#' or time-of-start, time-of-failure, and censoring in the case with left truncation at the time of begin
#' (censoring must be either 0 for no event or 1 for event);
#' \item  covariates (continuous or categorical) used in a study (can be empty set).}
#' @param expr The vector of expressions for contrasts. NULL, otherwise.
#' @param strata  Full list of the factors used in the study and their levels in the form:
#' strata=list(factor1=level of the factor1,...,last factor=level of the last factor),
#' NULL, otherwise.
#' @details Two kinds of the baseline hazards are used in this function:
#' \enumerate{
#' \item Weibull with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=(t/a)^b;}
#' \item Gompertz with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=(a/b)(exp(bt)-1)).}
#' }
#' The cumulative hazard function for the vector of covariates \emph{\strong{u}} is defined by
#' \deqn{H_{cum}(t;a,b,\beta _{scale},\beta _{shape},u)=exp(\beta _{scale}u)H_{0}(t;a,exp(\beta _{shape}u)b).}
#' The marginal survival function is defined by
#' \deqn{S(t;a,b,\beta _{scale},\beta _{shape},\sigma ^2,u)=Eexp(-ZH_{cum}(t;a,b,\beta _{scale},\beta _{shape},u))=(1+\sigma ^2H_{cum}(t;a,b,\beta _{scale},\beta _{shape},u))^{-1/\sigma ^2}}
#' for gamma-distributed frailty Z with mean 0 and variance \eqn{\sigma ^2}.
#' The 'formula.scale' and 'formula.shape' are formula objects used in the R-package 'survival'
#' and have the form \deqn{Surv(time, Cens) ~ factor_1+...+factor_k).}
#' Interactions between factors are allowed.
#' If no factors are used in the Cox regression the value of 1 stands on
#' the right side of the formula. The records with NA value for at least one factor in
#' both formulas are excluded from the analysis.
#' \deqn{}
#' Remark 1. The concordance is evaluated using the function 'survConcordance' from
#' the R-package 'survival'.
#' \deqn{}
#' Remark 2. The mean contrasts, CIs, and p-values are calculated on the
#' basis of the empirical distribution for contrasts constructed using \eqn{10^6} of the generations
#' of the vector of parameters.
#' \deqn{}
#' Remark 3. In some cases (flat likelihood function, multimodality, etc.) the hessian
#' cannot be correctly estimated and the function issues the following error message:
#'   \deqn{Error in ParNPHCox(formula.scale, formula.shape, cluster, dist, data = ...) :
#'   hessian cannot be correctly calculated.}
#' Change the model and try again.
#' It can occur if, for example, the data set is not large enough the estimates of the shape and scale
#' parameters are poorly calculated. In this case it
#' is recommended to simplify the model excluding some factors from the formulas for
#' shape or scale and try to calculate the estimates again.
#' \deqn{}
#' Remark 4. The number of factors can increase after conversion of any categorical factor in binary ones. For example, if factor "type" has three levels "A", "B",
#' and "C" we will get after conversion two binary factors - "typeB" and "typeC" (factor "typeA" is the baseline one and does not appear in the list of binary factors).
#' The command "strata=list(...)" must include at most one non-zero level for each
#' non-baseline level of the converted categorical factor. For example, it is correct to
#' write "strata=list(typeB=1,typeC=0,...)" or briefly "strata=list(typeB=1,...)" and
#' it is not correct to write "strata=list(typeB=1,typeC=1,...)".
#' \deqn{}
#' Remark 5. If the covariate is a numerical one but we want to consider it as a categorical variable we can convert it to categorical one using the function "as.factor(variable)".
#'
#' @return List containing the following components:
#' \itemize{
#'  \item par: parameter estimates
#' \item  se: standard errors for parameter estimates
#' \item  LogLik: the value of the loglikelihood
#' \item  Tabs: the table of parameter estimates, their standard errors, and p-values in the Latex format
#' \item  Names: the names of estimated parameters
#' \item  Conc: concordance and its standard error
#' \item  pval: the vector of p-values for parameter estimates.
#' For null hypothesis the values of parameters are equal to zero
#' \item  p.contrast: the data frame for means, CIs, and p-values for contrasts
#' \item  pstrata: data frame for times, means, and 95% CIs of the marginal survivals,
#' cumulative hazards and instant hazards. The column names are:
#' "Time", "Mean.survival", "Low.survival", "Upper.survival",
#' "Mean.cumulative.hazard", "Low.cumulative.hazard", "Upper.cumulative.hazard",
#' "Mean.hazard", "Low.hazard", "Upper.hazard". The duplicates are excluded
#' }
#'
#' @import MASS
#' @import survival
#' @import numDeriv
#' @import ucminf
#' @import stats
#' @import xtable
#'
#' @examples
#' ########## Example with 'Weibull baseline hazard function ############
#' require(survival)
#' require(numDeriv)
#' data("lung",package="survival")      #download data set "lung" from the r-package "survival"
#' lung$sex=lung$sex-1                  #converting numerical variable 'sex' with values 1 and 2
#' #to a binary one with values 0 and 1
#' lung$status=lung$status-1            #status variable must be 0 or 1
#' lung$ph.ecog=as.factor(lung$ph.ecog) #convering categorical variable in binary ones
#' formula.scale=as.formula('Surv(time, status)  ~ age + sex*ph.ecog')
#' formula.shape=as.formula('Surv(time, status) ~ 1') # no dependency
#' cluster='ph.karno'
#' dist='Weibull'
#' expr=expression(a -500* b, log(msurv(19)/msurv(14))-1, age.scale,sex.scale)
#' NamFact(lung,formula.scale,formula.shape)
#' #[1] "age" "ph.ecog1" "ph.ecog2" "ph.ecog3" "sex" "sex:ph.ecog1" "sex:ph.ecog2" "sex:ph.ecog3"
#' strata=list(age=0,sex=1,ph.ecog1=1)
#' dcWeibull=ParNPHCox(formula.scale,formula.shape,cluster,dist,data=lung,expr,strata)
#'
#' ########## Example with 'Gompertz baseline hazard function ############
#' formula.scale=as.formula('Surv(time, status) ~ sex')
#' formula.shape=as.formula('Surv(time, status) ~ sex')
#' cluster='pat.karno'
#' dist='Gompertz'
#' expr=expression(a - b, msurv(19) - msurv(14), sex.scale + sex.shape)
#' NamFact(lung,formula.scale,formula.shape)
#' #[1] "sex"
#' strata=list(sex=1)
#' dcGompertz=ParNPHCox(formula.scale,formula.shape,cluster,dist,data=lung,expr,strata)
#'
#' @export
#'
ParNPHCox=function(formula.scale, formula.shape,cluster,dist,data,expr,strata)
{
  nr=nrow(data)
  obsdata <- NULL
  if (length(formula.scale[[2]]) == 3) {
    obsdata$trunc <- c(rep(0,nrow(data)))
    obsdata$time <- eval(formula.scale[[2]][[2]], envir = data)
    obsdata$event <- eval(formula.scale[[2]][[3]], envir = data)
  }   else if (length(formula.scale[[2]]) == 4) {
    obsdata$trunc <- eval(formula.scale[[2]][[2]], envir = data)
    obsdata$time <- eval(formula.scale[[2]][[3]], envir = data)
    obsdata$event <- eval(formula.scale[[2]][[4]], envir = data)
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
    stop(paste("The status indicator 'event' in the Surv object",
               "in the left-hand side of the formula object", "must be either 0 (no event) or 1 (event)."))
  }
  obsdata$x <- as.data.frame(model.matrix.lm(formula.scale, data = data,na.action='na.pass'))
  obsdata$xs <- as.data.frame(model.matrix.lm(formula.shape, data = data,na.action='na.pass')) #factors for shape
  names(obsdata$x)=paste(names(obsdata$x),"scale",sep=".")
  names(obsdata$xs)=paste(names(obsdata$xs),"shape",sep=".")
  ind.x=which(is.na(c(apply(obsdata$x,1,sum))))
  ind.xs=which(is.na(c(apply(obsdata$xs,1,sum))))

  if (is.null(cluster)) {
    obsdata$ncl <- 0
    obsdata$di <- sum(obsdata$event)
    obsdata$cluster <- c(rep(0,nrow(data)))
    ind.cl <- as.numeric({})
  }   else {
    if (!cluster %in% names(data)) {
      stop(paste0("object '", cluster, "' not found"))
    }
    obsdata$cluster <- eval(as.name(cluster), envir = data)
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    obsdata$di <- aggregate(obsdata$event, by = list(obsdata$cluster),
                            FUN = sum)[, , drop = FALSE]
    cnames <- obsdata$di[, 1]
    obsdata$di <- as.vector(obsdata$di[, 2])
    names(obsdata$di) <- cnames
    ind.cl=which(is.na(obsdata$cl))
  }
  ind=sort(unique(c(ind.x,ind.xs,ind.cl)))
  nx=rep(0,nr)
  nx[ind]=1
  if (is.factor(obsdata$cluster)) obsdata$cluster=as.character(obsdata$cluster)
  obs=data.frame(obsdata$xs[-1],obsdata$x[-1],obsdata$event,obsdata$trunc,obsdata$time,obsdata$cluster)[nx!=1,]
  nr=nrow(obs)
  namesk=names(obsdata$xs)[-1]
  namesf=names(obsdata$x)[-1]
  names(obs)=c(names(obsdata$xs)[-1],names(obsdata$x)[-1],names(obsdata$event),"event","trunc","time","cluster")
  ncl=obsdata$ncl
  nk=length(namesk)
  nf=length(namesf)
  D=obs
  par0=c(0,0)
  Result=ucminf(par=par0,LikGenNPH,gr=NULL,D=obs,nf=0,nk=0,ncl=0,dist=dist,hessian=1)
  if (ncl>0){
    par0=c(Result$par,rep(0,(1+nk+nf)))} else {
      par0=c(Result$par,rep(0,(nk+nf)))
    }
  Result=ucminf(par=par0,LikGenNPH,gr=GrGenNPH,D=obs,nf=nf,nk=nk,ncl=ncl,dist=dist,hessian=1)
  par=Result$par
  if (any(!is.finite(as.matrix(Result$hessian))))
    stop("infinite or missing values in hessian. It is not possible to calculate the matrix of covariates. \n  Change the model and try again.")
  if (any(suppressWarnings(diag(ginv(Result$hessian)))<0))
    stop("hessian cannot be correctly calculated. \n  Change the model and try again.")
  invHes=sqrt(diag(ginv(Result$hessian)))
  Lik=-Result$value
  negH=-cumhazard(obs,nf,nk,ncl,dist,matrix(Result$par,1,length(Result$par)),c(1:nrow(obs)))
  HCT=data.frame(obs$time,negH,obs$event)
  colnames(HCT)=c("Time","H","Cens")
  CSE=survConcordance(Surv(Time,Cens) ~ H,HCT)
  Conc=c(as.numeric(CSE$concordance),as.numeric(CSE$std.err))
  Vnames1=c("a","b",namesk,namesf)
  Vnames={}
  if ((nk+nf)>0)    Vnames=paste("exp(",c(namesk,namesf),")",sep="")
  if (ncl>0) Vnames1=c(Vnames1,"Sigma2")

  if (ncl==0 & dist=='Weibull'){
    Names=c("Sample size","Number of non-censored","a","b",Vnames,"Concordance (se)","Loglik","AIC")} else if (ncl>0 & dist=='Weibull'){
      Names=c("Sample size","Number of non-censored","a  ","b",Vnames,"Sigma2","Concordance (se)","Loglik","AIC")} else if (ncl==0 & dist=='Gompertz'){
        Names=c("Sample size","Number of non-censored","1000a ","100b",Vnames,"Concordance (se)","Loglik","AIC")} else if (ncl>0 & dist=='Gompertz'){
          Names=c("Sample size","Number of non-censored","1000a ","100b",Vnames,"Sigma2","Concordance (se)","Loglik","AIC")
        }
  set.seed(123)
  vpar<- mvrnorm(n = 1000000, Result$par, ginv(Result$hessian), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  p_h={}
  p_hl={}
  p_hu={}
  for (lika in 1:length(par)){
    hh0=ecdf(vpar[,lika])
    if (lika<3 | (lika==length(par) & ncl>0)) {
      p_h=c(p_h,round(hh0(-Inf),4))} else {
        p_h=c(p_h,round(2*min(hh0(0),1-hh0(0)),4))
      }
    p_hl=c(p_hl,round(exp(quantile(hh0,probs=0.025)),3))
    p_hu=c(p_hu,round(exp(quantile(hh0,probs=0.975)),3))
  }
  parval=c(as.character(nr),as.character(sum(obs$event)),round(exp(Result$par),3),paste0(round(Conc[1],3)," (",round(Conc[2],3),")"),round(-Result$value,2),round(2*(Result$value+length(par0)),2))
  parvalm=c("","",p_hl,"","","")
  parvalp=c("","",p_hu,"","","")
  CI=c("","",paste(parvalm[3:(2+length(par))],parvalp[3:(2+length(par))],sep=" - "),"","","")
  pval=c("","",p_h,"","","")
  pval0=as.numeric(pval[3:(2+length(par0))])
  Tab=data.frame(parval,CI,pval)
  colnames(Tab)=c("Estimates","CI","p-value")
  rownames(Tab)=Names
  capt=paste("Parameter estimates.",dist,"model.",sep=" ")
  print(xtable(Tab,caption=capt))
  pcon=NULL
  ru=NULL
  if (!is.null(expr)){
    if(is.expression(expr)){
      pcon_={}
      for (iex in 1:length(expr)){
        for (ln in 1:length(par)){
          eval(parse(text=paste0("ru","=vpar[,ln]")))
          if (dist=='Weibull'){
            assign(Vnames1[ln],ru)
            if (ln<3 | (ln==length(par) & ncl>0))    assign(Vnames1[ln],exp(ru))
          }
          if (dist=='Gompertz'){
            assign(Vnames1[ln],ru)
            if (ln==length(par) & ncl>0)    assign(Vnames1[ln],exp(ru))
            if (ln==1)    assign(Vnames1[ln],1e-3*exp(ru))
            if (ln==2)    assign(Vnames1[ln],1e-2*exp(ru))
          }
        }
        expr1=gsub("msurv(","msurv(obs,nf,nk,ncl,dist,vpar,",deparse(expr[iex]),fixed = TRUE)
        expr1=gsub("hazard(","hazard(obs,nf,nk,ncl,dist,vpar,",expr1,fixed = TRUE)
        expr1=gsub("cumhazard(","cumhazard(obs,nf,nk,ncl,dist,vpar,",expr1,fixed = TRUE)

        hi=ecdf(eval(eval(parse(text=expr1))))
        expriex=eval(eval(parse(text=expr1)))
        mea=mean(expriex)
        CIL=as.numeric(quantile(expriex,probs=0.025))
        CIU=as.numeric(quantile(expriex,probs=0.975))
        p_expr=2*min(hi(0),1-hi(0))
        pcon_=rbind(pcon_,c(round(mea,4),paste0(as.character(round(CIL,4)),"-",as.character(round(CIU,4))),as.character(round(p_expr,4)),deparse(expr[iex][[1]])))
      }
      pcon=data.frame(pcon_[,1:3])
      colnames(pcon)=c("contrast","CI","p-value")
      rownames(pcon)=pcon_[,4]
      capt=paste("Table of contrasts.",dist,"model.",sep=" ")
      print(xtable(pcon,caption=capt))
    } else {
      print(c("Variable 'expr' must have type 'expression'"))
    }
  }
  vpar<- mvrnorm(n = 10000, Result$par, ginv(Result$hessian), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  pstrata=NULL
  if (!is.null(strata)){
    cn=colnames(obs)
    cs=names(strata)
    psurv_={}
    pcumhaz_={}
    phaz_= {}
    tic=unique(obs$time)
    D1<-obs[1:length(tic),]
    D1$time=tic
    if ((nk+nf)>0) {
      D1[,1:(nk+nf)]=0
      for (ik in 1:(nk+nf)){
        ix=which(cs==substr(cn[ik],1,(nchar(cn[ik])-6)))
        if (length(ix)>0) D1[,ik]=as.numeric(strata[ix])
      }
      for (ik in 1:(nk+nf)){
        ss=1
        iv=0
        for (ix in 1:length(cs)){
          if (regexpr(cs[ix],cn[ik])[1]>0) {
            ss=ss*as.numeric(strata[ix])
            iv=iv+1
          }
        }
        if (iv>0) D1[,ik]=ss
      }
    }
    D1<-D1[order(D1$time),]
    Time=D1$time
    for (ln in 1:length(par)){
      eval(parse(text=paste0("ru","=vpar[,ln]")))
      if (dist=='Weibull'){
        assign(Vnames1[ln],ru)
        if (ln<3 | (ln==length(par) & ncl>0))    assign(Vnames1[ln],exp(ru))
      }
      if (dist=='Gompertz'){
        assign(Vnames1[ln],ru)
        if (ln==length(par) & ncl>0)    assign(Vnames1[ln],exp(ru))
        if (ln==1)    assign(Vnames1[ln],1e-3*exp(ru))
        if (ln==2)    assign(Vnames1[ln],1e-2*exp(ru))
      }
    }
    for (iy in 1:nrow(D1)){
      expriex=eval(expression(msurv(D1,nf,nk,ncl,dist,vpar,iy)))
      mea=mean(expriex)
      med=median(expriex)
      CIL=as.numeric(quantile(expriex,probs=0.025))
      CIU=as.numeric(quantile(expriex,probs=0.975))
      psurv_=rbind(psurv_,c(round(mea,4),round(med,4),round(CIL,4),round(CIU,4)))
      colnames(psurv_)=c("Mean survival","Median survival","Low survival","Upper survival")
      expriex=eval(expression(cumhazard(D1,nf,nk,ncl,dist,vpar,iy)))
      mea=mean(expriex)
      med=median(expriex)
      CIL=as.numeric(quantile(expriex,probs=0.025))
      CIU=as.numeric(quantile(expriex,probs=0.975))
      pcumhaz_=rbind(pcumhaz_,c(round(mea,4),round(med,4),round(CIL,4),round(CIU,4)))
      colnames(pcumhaz_)=c("Mean cumulative hazard","Median cumulative hazard","Low cumulative hazard","Upper cumulative hazard")
      expriex=eval(expression(hazard(D1,nf,nk,ncl,dist,vpar,iy)))
      mea=mean(expriex)
      med=median(expriex)
      CIL=as.numeric(quantile(expriex,probs=0.025))
      CIU=as.numeric(quantile(expriex,probs=0.975))
      phaz_=rbind(phaz_,c(round(mea,4),round(med,4),round(CIL,4),round(CIU,4)))
      colnames(phaz_)=c("Mean hazard","Median hazard","Low hazard","Upper hazard")
    }
    Time=sort(Time)
    pstrata=data.frame(Time,psurv_,pcumhaz_,phaz_)
    rownames(pstrata)<-NULL
  }

  list(par=par,se=invHes,LogLik=Lik,Tab=Tab,Names=Vnames1,Conc=Conc,pval=pval0,p.contrast=pcon,pstrata=pstrata)
}
