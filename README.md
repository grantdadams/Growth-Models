# Growth-Models
Some code examples of simulating length-at-age data and fitting von Bertallanfy growth models in TMB, JAGS, and Stan that I wrote while studying spatially-varying growth of Sheepshead. Stan models used in the final publication are included, but see [here](https://github.com/grantdadams/Spatial-Growth-Models) for full analysis. The following models are included in the R folder:

1. Maximum likelihood three parameter von Bertallanfy model fit in TMB [link](https://github.com/grantdadams/Growth-Models/blob/master/R/von_Bertalanffy_growth_TMB.R)
2. Maximum likelihood three parameter von Bertallanfy model fit in TMB where the likelihood function is a truncated normal distribution (upper and lower truncation points [link](https://github.com/grantdadams/Growth-Models/blob/master/R/Truncated_von_Bertalanffy_growth_TMB.R)
3. Bayesian likelihood three parameter von Bertallanfy model fit in JAGS with normally distributed population level parameters and global hyperparameters [link](https://github.com/grantdadams/Growth-Models/blob/master/R/Hierarchical_growth_model_JAGS.R)
4. Bayesian likelihood three parameter von Bertallanfy model fit in JAGS with multivariate normal population and individual level parameters and global hyperparameters and environmental covariates [link](https://github.com/grantdadams/Growth-Models/blob/master/R/Hierarchical_growth_model_JAGS.R). The model is sructured similar to [Helser and Lai (2004)](https://www.sciencedirect.com/science/article/pii/S0304380004001577), but accounts for repeated sampling and sex/location effects.

Citation:
Adams, G.D., Leaf, R.T., Ballenger, J.C., Arnott, S.A., Mcdonough, C.J., 2018. Spatial variability in the growth of Sheepshead (Archosargus probatocephalus) in the Southeast US : Implications for assessment and management. Fish. Res. 206, 35–43. [doi:10.1016/j.fishres.2018.04.023](https://www.sciencedirect.com/science/article/abs/pii/S0165783618301279)
