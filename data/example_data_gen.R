n <- 10000
set.seed(4399)
expit <- function(x) 1/(1+exp(-x))

#### Generate covariates
X1 <- rnorm(n)
X2 <- rnorm(n, sd=2)
X3 <- rt(n, df=3)
X4 <- rbeta(n, shape1=0.5, shape2=1.5)
X5 <- rbinom(n, 10, 0.3)
X6 <- X1*X5
X7 <- X4^2
X8 <- X2*X6

X <- cbind(X1,X2,X3,X4,X5,X6,X7,X8)

#### Generate the treatment
beta <- 0.015*(1:8)
ps <- expit(-0.2+X%*%beta)
# plot(ps)
# mean(ps)
A <- rbinom(n, 1, prob=ps)
# mean(A)

#### Generate the outcomes
alpha <- 0.05*(1:8)
# constant outcome
Y.c <- 45 + X%*%alpha + 10*A + rnorm(n)
# heterogeneous outcome
Y.h <- 45 + X%*%alpha + (10/3)*A*(X1-X2+X5) + rnorm(n)

df.homo <- data.frame(X, Y=Y.c, A=A)
df.hete <- data.frame(X, Y=Y.h, A=A)

save(file="example.Rdata", df.hete, df.homo)
