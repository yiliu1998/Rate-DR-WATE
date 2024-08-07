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

Next, we can implement our package using code as follows. We consider evaluating the following estimands: ATE (overall), ATT (treated), ATC (controls), ATO (overlap weights), ATM (matching weights), ATEN (entropy), and two ATB's (beta weights, with $\nu_1=\nu_2\in\{3,5\}$). As an example, we consider applying our method on the dataset `df.hete`. We include `SL.glm`, `SL.glm.interaction` and `SL.glmnet` in the ensemble library of `SuperLearner` for getting our nuisance function estimates. 

As an important remark, we do not currently allow **missing data** (in all covariates, treatment and outcome) in the dataset. Therefore, if your data have missing values, we suggest either doing a complete-case analysis or imputing missing data first. Once the data is cleaned and prepared, it is ready to run our methods. The first step is to extract the binary treatment vector, outcome vector and covariates matrix (can also be a data frame). 

```{r}
A <- df.hete$A
Y <- df.hete$Y
X <- df.hete %>% select(-A, -Y)
```

Then, plug-in these arguments to our main function `RDRwate()`: 

```{r}
if(1==1)
v1 <- c(3,5)
v2 <- c(3,5)
beta=TRUE
result.eif <- RDRwate(A=A, Y=Y, X=X, beta=beta, v1=v1, v2=v2,
                      method="EIF", 
                      ps.library=c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
                      out.library=c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
                      seed=1)
print(result.eif)
```

```{r}
weights       Est    Std.Err
1   overall 10.234560 0.09159221
2   treated 10.717672 0.12000207
3   control  9.729105 0.13605199
4   overlap 10.160023 0.10447144
5  matching 10.058661 0.13527458
6   entropy 10.174391 0.09968131
7 beta(3,3) 10.061675 0.53400202
8 beta(5,5) 10.029801 0.87891157
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
print(result.dml)
```

```{r}
$result.dml.1
    weights       Est   Std.Err
1   overall 10.236479 0.0916399
2   treated 10.722442 0.1200593
3   control  9.723038 0.1359380
4   overlap 10.158679 0.1051220
5  matching 10.042332 0.1356507
6   entropy 10.175380 0.1002631
7 beta(3,3) 10.094221 0.5432445
8 beta(5,5) 10.123608 0.9244998

$result.dml.2
    weights       Est    Std.Err
1   overall 10.236479 0.09166915
2   treated 10.721589 0.12013919
3   control  9.728936 0.13602753
4   overlap 10.160243 0.10514657
5  matching 10.044138 0.13573247
6   entropy 10.176528 0.10029061
7 beta(3,3) 10.106283 0.53525898
8 beta(5,5) 10.146469 0.88103159
```

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions).

## Reference
Please cite the following paper:

Wang, Y., Liu, Y., & Yang, S. (2024). Rate doubly robust estimation for weighted average treatment effects. Under review. 
