## nntcalc: The Number Needed to Treat (NNT) calculator 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17068632.svg)](https://doi.org/10.5281/zenodo.17068632)

[Valentin Vancak](https://www.linkedin.com/in/valentin-vancak-0a56227a/?originalSubdomain=il)

The `nntcalc` R-package provides functions to calculate the unadjusted, the conditional (adjusted), and the harmonic mean (marginal) Laupacis NNT with the corresponding 95% confidence intervals. Available regression models include ANOVA, linear regression and logistic regression, and Cox models. In addition, the package provides a function to calculate the estimators of the Kraemer & Kupfer's NNT. 

For the installation of `nntcalc` you need to load first the `devtools` package by `library(devtools)`, and then type

`install_github("https://github.com/vancak/NNTcalculator")`

Tha package is called `nntcalc`. To load it after installation type

`library(nntcalc)`

For further details and examples please see [package manual](https://github.com/vancak/NNTcalculator/blob/main/manual.pdf). 

*NOTE*: The source code for the simulations presented in the manuscript "The Number Needed to Treat Adjusted for Explanatory Variables in Regression and Survival Analysis: Theory and Application" is available [here](https://github.com/vancak/nntx_simulations).

**Cite:** Vancak V. *nntcalc: The Number Needed to Treat (NNT) calculator* (v1.1). Zenodo.  
Version DOI: https://doi.org/10.5281/zenodo.17068632  
