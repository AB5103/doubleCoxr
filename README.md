# doubleCoxr, Parametric Gamma-Frailty Models with Non-Proportional Hazard Functions in ```R``` 
Alexander Begun and Felix Begun
## Description
Fits Parametric Gamma-Frailty Models with Non-Proportional Hazard Functions by maximum marginal likelihood. 
Possible baseline hazards: Weibull, Gompertz. The proportional hazard assumption is violated by using the
second Cox-regression term relating to the shape parameter. The first Cox-regression term is the Cox-regression term 
used in proportional hazards model.

It is assumed that all subjects in a cluster share the same unobserved risk of failure (frailty) that is 
a gamma distributed random variable. The Weibull and the Gompertz baseline hazard functions are implemented.
