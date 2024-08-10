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

```r
library(RDRwate)
```

The following packages are also required:

```r
library(dplyr)
library(ggplot2)
library(SuperLearner)
library(glmnet)
library(caret)
```

Running the following code, two sets of pre-generated observational data (`df.homo` and `df.hete`) will show up in the R environment. 

```r
load(system.file("data", "example.Rdata", package="RDRwate"))
```

A quick look to the generated data:

```r
head(df.hete)
```

```r
X1         X2         X3          X4 X5           X6           X7          X8
1 -0.564754438 -0.8032184  1.6256195 0.050354579  4 -2.259017750 2.535584e-03  1.81448472
2  1.299500029  0.4762622  0.4348662 0.301498110  2  2.599000059 9.090111e-02  1.23780539
3  0.001802605  3.4165443  0.1118126 0.444590718  2  0.003605209 1.976609e-01  0.01231736
4 -0.570083502  2.6136722 -1.0676598 0.545845427  2 -1.140167005 2.979472e-01 -2.98002280
5 -0.500492410 -1.0546837 -0.3603949 0.004596238  1 -0.500492410 2.112541e-05  0.52786117
6  0.512164377  2.1591574  0.3206966 0.266608928  1  0.512164377 7.108032e-02  1.10584351
         Y A
1 46.70110 0
2 45.36502 0
3 47.03340 0
4 39.80492 1
5 50.92973 1
6 45.43626 0
```

We also allow users to specify different methods for nuisance function estimation and prediction in our RDR estimators of WATE. The complete list of these methods can be found using the following code. 

```r
SuperLearner::listWrappers()
```

Next, we can implement our package using code as follows. We consider evaluating the following estimands: ATE (overall), ATT (treated), ATC (controls), ATO (overlap weights), ATM (matching weights), ATEN (entropy), and two ATB's (beta weights, with $\nu_1=\nu_2\in\{3,5\}$). As an example, we consider applying our method on the dataset `df.hete`. We include `SL.glm`, `SL.glm.interaction` and `SL.glmnet` in the ensemble library of `SuperLearner` for getting our nuisance function estimates. 

As an important remark, we do not currently allow **missing data** (in all covariates, treatment and outcome) in the dataset. Therefore, if your data have missing values, we suggest either doing a complete-case analysis or imputing missing data first. Once the data is cleaned and prepared, it is ready to run our methods. The first step is to extract the binary treatment vector, outcome vector and covariates matrix (can also be a data frame). 

```r
A <- df.hete$A
Y <- df.hete$Y
X <- df.hete %>% select(-A, -Y)
```

Then, plug-in these arguments to our main function `RDRwate()`: 

```r
v1 <- c(3,5)
v2 <- c(3,5)
beta=TRUE
result.eif <- RDRwate(A=A, Y=Y, X=X, beta=beta, v1=v1, v2=v2,
                      method="EIF", 
                      ps.library=c("SL.glm", "SL.glmnet"),
                      out.library=c("SL.glm", "SL.glmnet"),
                      seed=1)
print(result.eif)
```

```r
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

We apply the DML-based methods using 10 sample splits (i.e., repeating cross-fitting 10 times by different random sample splitting) and 5 folds for cross-fitting. 

```r
v1 <- c(3,5)
v2 <- c(3,5)
beta=TRUE
result.dml <- RDRwate(A=A, Y=Y, X=X, beta=beta, v1=v1, v2=v2,
                      method="DML", n.folds=5, n.split=10, 
                      ps.library=c("SL.glm", "SL.glmnet"),
                      out.library=c("SL.glm", "SL.glmnet"),
                      seed=1)
print(result.dml)
```

```r
$result.dml.1
    weights       Est   SE.mean  SE.median
1   overall 10.238393 0.0917362 0.09172764
2   treated 10.725430 0.1203668 0.12035407
3   control  9.727280 0.1360394 0.13605126
4   overlap 10.162204 0.1056018 0.10558813
5  matching 10.058743 0.1360952 0.13607341
6   entropy 10.176439 0.1007624 0.10072919
7 beta(3,3) 10.041305 0.5411684 0.53918889
8 beta(5,5)  9.989887 0.9185246 0.90626480

$result.dml.2
    weights       Est    SE.mean  SE.median
1   overall 10.235929 0.09182739 0.09195195
2   treated 10.689941 0.12135043 0.12037325
3   control  9.758034 0.13511888 0.13590928
4   overlap 10.150839 0.10492282 0.10503028
5  matching 10.038948 0.13566287 0.13621179
6   entropy 10.166876 0.10021070 0.10029119
7 beta(3,3)  9.789291 0.53855173 0.53947741
8 beta(5,5)  9.422549 0.90116910 0.89831503
```

The printed results contain both mean and median strategies for standard error estimation, based on <a href="https://www.aeaweb.org/articles?id=10.1257/aer.p20171038">Chernozhukov et al. (2017)</a>. In this example, DML-1 and DML-2 do not make too obvious difference, and both mean and median stratgies for the standard error are also similar. The two DML-based estimators also have similar results to the EIF-based estimator on all WATEs. 

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions).

## Reference
Please cite the following paper:

Wang, Y., Liu, Y., & Yang, S. (2024). Rate doubly robust estimation for weighted average treatment effects. Under review. 
