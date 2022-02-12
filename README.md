# doubleCoxr, Parametric Gamma-Frailty Models with 
# Non-Proportional Hazard Functions in ```R``` 
Alexander Begun and Felix Begun
## Description
Fits Parametric Gamma-Frailty Models with Non-Proportional Hazard Functions by maximum marginal likelihood.
Provides parameter estimates, their standard errors and *p*-values, log-likelihood, concordance, contrasts (optional), 
and other characteristics (optional). Possible baseline hazards: Weibull, Gompertz.
## Details
The proportional hazard assumption is violated by using the
second Cox-regression term multiplying the shape parameter. The first Cox-regression term is one used in 
the Cox proportional hazards model as hazard multiplier.

It is assumed that all subjects in a cluster share the same unobserved risk of failure (frailty) that is 
a gamma distributed random variable. The Weibull and the Gompertz baseline hazard functions are implemented.
