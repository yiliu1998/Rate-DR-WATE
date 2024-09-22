# WATE
`WATE` is an R package implementing rate doubly robust (RDR) estimators and their variance estimations for weighted average treatment effects (WATEs). The WATE is a class of general causal estimands targeting different population defined by propensity score weights (see [Li et al.(2018)](https://www.tandfonline.com/doi/full/10.1080/01621459.2016.1260466) as a representative paper for this concept). The package contains three RDR estimators and two classical (naive) estimators. We allow users implement the estimation and inference for a number of popular WATEs, including ATE, ATT, ATC, ATO (overlap weights), ATM (matching weights), ATEN (entropy weights), and ATB (beta family weights). The three proposed RDR estimators include: 

* the efficient influence function (EIF)-based estimator; and
* two double/debiased machine learning (DML)-based estimators (denoted as DML-1 and DML-2), which are technically by applying the sample-splitting and cross-fitting ([Chernozhukov et al. (2018)](https://academic.oup.com/ectj/article/21/1/C1/5056401)) to the EIF-based estimator.

The two classical (naive) estimators are the propensity score weighting estimator and the outcome imputation estimator. They are referred to as naive-1 and naive-2 estimators in our paper, respectively. 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("yiliu1998/WATE")
```

## Demonstration
You can download and run this Rmd file ([click here](https://github.com/yiliu1998/Rate-DR-WATE/tree/main/vignettes)) in your R Studio after downloading the package, which gives an illustrative example of our package. We also demonstrate the R code in it as follows. 

Load our package: 

```r
library(WATE)
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

Next, we can implement our package using code as follows. We consider evaluating the following estimands: ATE (overall), ATT (treated), ATC (controls), ATO (overlap weights), ATM (matching weights), ATEN (entropy), and two ATB's (beta weights, with $\nu_1=\nu_2\in$ {3,5}). As an example, we consider applying our method on the dataset `df.hete`. We include `SL.glm`, `SL.glm.interaction` and `SL.glmnet` in the ensemble library of `SuperLearner` for getting our nuisance function estimates. 

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

We apply the two DML-based methods using a 5-fold cross-fitting with 10 different sample-splittings (i.e., repeating cross-fitting 10 times by **different** random sample splits). 

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
    weights  Est.mean    SE.mean Est.median  SE.median
1   overall 10.237493 0.09167897  10.237459 0.09167677
2   treated 10.723602 0.12023261  10.723158 0.12022327
3   control  9.727349 0.13599869   9.727354 0.13600613
4   overlap 10.162038 0.10496583  10.161708 0.10496063
5  matching 10.058577 0.13540477  10.052414 0.13542817
6   entropy 10.176429 0.10022423  10.176269 0.10020905
7 beta(3,3) 10.048625 0.54010367  10.034852 0.53827577
8 beta(5,5)  9.990965 0.91675641   9.958111 0.90402204

$result.dml.2
    weights  Est.mean    SE.mean Est.median  SE.median
1   overall 10.234441 0.09174239  10.210011 0.09187881
2   treated 10.686658 0.12112849  10.655065 0.12017883
3   control  9.758435 0.13507440   9.790092 0.13587548
4   overlap 10.150758 0.10414168  10.118703 0.10435990
5  matching 10.034635 0.13498922  10.068590 0.13511639
6   entropy 10.166961 0.09954993  10.152039 0.09960444
7 beta(3,3)  9.796277 0.53732307   9.738101 0.53824587
8 beta(5,5)  9.417174 0.89896963   9.521300 0.89602009
```

The printed results of DML-based methods contain both mean and median-correction strategies for standard error estimation, which is based on the practical recommendation in [Chernozhukov et al. (2018)](https://academic.oup.com/ectj/article/21/1/C1/5056401). In this example, DML-1 and DML-2 do not make too outstanding difference, and both mean and median stratgies for the standard error are also similar. However, in some cases, the median strategy should be more robust because there could be extreme variances in a small number of sample-splittings. The two DML-based estimators also have similar results to the EIF-based estimator on all WATEs on this data set. 

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out yi.liu.biostat@gmail.com, if you have any questions).

## Reference
Please cite the following paper:

Wang, Y., Liu, Y., & Yang, S. (2024). Rate doubly robust estimation for weighted average treatment effects. Under review. 
