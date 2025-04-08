# Growth-Models
Some code examples of simulating length-at-age data and fitting von Bertallanfy growth models in TMB, JAGS, and Stan that I wrote while studying spatially-varying growth of Sheepshead. Stan models used in the final publication are included, but see [here](https://github.com/grantdadams/Spatial-Growth-Models) for full analysis. The following models are included in the R folder:

1. Maximum likelihood three parameter von Bertallanfy model fit in TMB [link](https://github.com/grantdadams/Growth-Models/blob/master/R/1_von_Bertalanffy_growth_TMB.R)
2. Maximum likelihood three parameter von Bertallanfy model fit in TMB where the likelihood function is a truncated normal distribution (upper and lower truncation points [link](https://github.com/grantdadams/Growth-Models/blob/master/R/2_truncated_von_Bertalanffy_growth_TMB.R)
3. Bayesian likelihood three parameter von Bertallanfy model fit in JAGS with normally distributed group level parameters and global hyperparameters [link](https://github.com/grantdadams/Growth-Models/blob/master/R/3_hierarchical_growth_model_JAGS.R)
4. Bayesian likelihood three parameter von Bertallanfy model fit in JAGS with multivariate normal group and individual level parameters and global hyperparameters and environmental covariates [link](https://github.com/grantdadams/Growth-Models/blob/master/R/4_hierarchical_growth_model_w_individual_and_group_effects_JAGS.R). The model is sructured similar to [Helser and Lai (2004)](https://www.sciencedirect.com/science/article/pii/S0304380004001577), but accounts for repeated sampling and sex/location effects.
5. Bayesian likelihood three parameter von Bertallanfy model fit in Stan with normally distributed group level parameters and global hyperparameters [link](https://github.com/grantdadams/Growth-Models/blob/master/R/5_hierarchical_growth_model_w_group_effects_Stan.R)
6. Base R code to create a prior predictive distribution of a horshoe prior. Horseshoe priors are commonly recommended as shrinkage priors in Bayesian models.

TMB (Template Model Builder) is an R package for fitting statistical latent variable models to data. It is strongly inspired by ADMB. Unlike most other R packages the model is formulated in C++. This provides great flexibility, but requires some familiarity with the C/C++ programming language. See their [GitHub](https://github.com/kaskr/adcomp/wiki) for more information. You can also fit TMB models in [Stan](https://mc-stan.org/) using the [adnuts](https://github.com/Cole-Monnahan-NOAA/adnuts) package. Alternatively, [JAGS](https://mcmc-jags.sourceforge.io/) is Just Another Gibbs Sampler. [Stan](https://mc-stan.org/) is a more efficient sampler and my recommended program for Bayesian analysis.

Citation:
Adams, G.D., Leaf, R.T., Ballenger, J.C., Arnott, S.A., Mcdonough, C.J., 2018. Spatial variability in the growth of Sheepshead (Archosargus probatocephalus) in the Southeast US : Implications for assessment and management. Fish. Res. 206, 35–43. [doi:10.1016/j.fishres.2018.04.023](https://www.sciencedirect.com/science/article/abs/pii/S0165783618301279)
