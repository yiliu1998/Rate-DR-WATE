# RDRwate
Rate doubly robust (RDR) estimation for weighted average treatment effects (WATEs). We developed an R package (`RDRwate`) implementing three RDR estimators for WATE, a general class of causal estimands. We allow users implement the estimation and inference for a number of popular WATEs in the class, including ATE, ATT, ATC, ATO (overlap weights), ATM (matching weights), ATEN (entropy weights), and ATB (beta family weights). The proposed estimators include: 

* an efficient influence function (EIF)-based estimator; 
* two double/debiased machine learning (DML)-based estimators (DML-1 and DML-2), which actually applies the cross-fitting algorithm on the EIF-based estimator.

We also implement the variance estimations of these estimators in our package.  

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("yiliu1998/Rate-DR-WATE")
```

## Demonstration
You can download and run this Rmd file ([click here](https://github.com/yiliu1998/Rate-DR-WATE/tree/main/vignettes)) in your R Studio after downloading the package, which gives an illustrative example of our package. We also demonstrate the R code in it as follows. 

Load our package: 

```{r}
library(RDRwate)
```

The following packages are also required:

```{r}
library(dplyr)
library(ggplot2)
library(SuperLearner)
library(glmnet)
library(caret)
```

Running the following code, two sets of pre-generated observational data (`df.homo` and `df.hete`) will show up in the R environment. 

```{r}
load(system.file("data", "example.Rdata", package="RDRwate"))
```

We also allow users to specify different methods for nuisance function estimation and prediction in our RDR estimators of WATE. The complete list of these methods can be found using the following code. 

```{r}
SuperLearner::listWrappers()
```

Next, we can implement our package using code as follows. We consider evaluating the following estimands: ATE (overall), ATT (treated), ATC (controls), ATO (overlap weights), ATM (matching weights), ATEN (entropy), and two ATB's (beta weights, with $\nu_1=\nu_2\in\{3,5\}$). As an illustration, we apply our method on the dataset `df.hete`. We include `SL.glm`, `SL.glm.interaction` and `SL.glmnet` in the ensemble library of `SuperLearner` for running our methods. 

As an important remark, we do not currently allow missing data (in all covariates, treatment and outcome) in the dataset. Therefore, if your data have missing values, we suggest either a complete-case analysis or imputing missing data first. Once the data is cleaned and prepared, it is ready to run our methods. The first step is to extract the binary treatment vector, outcome vector and covariates matrix (can also be a data frame). 

```{r}
A <- df.hete$A
Y <- df.hete$Y
X <- df.hete %>% select(-A, -Y)
```

Then, plug-in these arguments to our main function `RDRwate()`: 

```{r}
v1 <- c(3,5)
v2 <- c(3,5)
beta=TRUE
result.eif <- RDRwate(A=A, Y=Y, X=X, beta=beta, v1=v1, v2=v2,
                      method="EIF", 
                      ps.library=c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
                      out.library=c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
                      seed=1)
result.eif
```

```{r}
v1 <- c(3,5)
v2 <- c(3,5)
beta=TRUE
result.dml <- RDRwate(A=A, Y=Y, X=X, beta=beta, v1=v1, v2=v2,
                      method="DML", n.folds=5,
                      ps.library=c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
                      out.library=c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
                      seed=1)
result.dml
```

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat[at]gmail[dot]com, if you have any questions).

## Reference
Please cite the following paper:

Wang, Y., Liu, Y., & Yang, S. (2024). Rate doubly robust estimation for weighted average treatment effects. Under review. 
