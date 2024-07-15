# RDRwate
Rate doubly robust (RDR) estimation for weighted average treatment effects (WATEs). We developed an R package (`RDRwate`) implementing three RDR estimators for WATE, a general class of causal estimands. We allow users implement the estimation and inference for a number of popular WATEs in the class, including ATE, ATT, ATC, ATO (overlap weights), ATM (matching weights), ATEN (entropy weights), and ATB (beta family weights). The proposed estimators include: (i) an efficient influence function (EIF)-based estimator; (ii) two double/debiased machine learning (DML)-based estimators (DML-1 and DML-2), which actually applies the cross-fitting algorithm on the EIF-based estimator. We also implement the variance estimations of these estimators in our package.  

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("yiliu1998/Rate-DR-WATE")
```

## Demonstration
You can download and run this Rmd file ([click here]) in your R Studio after downloading the package, which gives an illustrative example of our package.  

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat[at]gmail[dot]com, if you have any questions).

## Reference
Please cite the following paper:

Wang, Y., Liu, Y., & Yang, S. (2024). Rate doubly robust estimation for weighted average treatment effects. Under review. 
