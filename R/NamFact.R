#' @title NamFact
#'
#' @description This function returns the list of factor names (after possible converting categorical variables to multiple binary ones as well)
#'
#' @param formula.scale A formula object, with the response on the left of a ~ operator, and the terms on
#' the right. The response must be a survival object as returned by the Surv() function.
#' The status indicator 'event' in the Surv object must be 0=censored, 1=non-censored.
#' This object describes the effect of several factors in the Cox proportional hazards model.
#' @param formula.shape A formula object, with the response on the left of a ~ operator, and the terms on
#' the right. The response must be a survival object as returned by the Surv() function.
#' The status indicator 'event' in the Surv object must be 0=censored, 1=non-censored.
#' This object describes the effect of several factors on the shape parameter.
#' @param data A data.frame in which to interpret the variables named in the formula.
#' The data set includes the following fields:
#' \enumerate{
#' \item Time-to-failure and censoring in the case without left truncation
#' or time-of-start, time-of-failure, and censoring in the case with left truncation at the time of begin;
#' Censoring must be either 0 (no event) or 1 (event);
#' \item Covariates (continuous or categorical) used in a study (can be empty set).
#' }
#'
#' @return List of the factor names
#'
#' @examples
#' \dontrun{
#' NamFact(data, formula.scale, formula.shape)
#' }
#'
#' @export
#'
NamFact=function(data,formula.scale,formula.shape){
  nsc=names(as.data.frame(model.matrix.lm(formula.scale, data = data,na.action='na.pass')))[-1]
  nsh=names(as.data.frame(model.matrix.lm(formula.shape, data = data,na.action='na.pass')))[-1]
  if ((length(nsc)+length(nsh))>0) {
    return(sort(unique(c(nsc,nsh))))} else {
      return(NULL)}
}
