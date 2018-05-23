Parameter Estimation and Machine Learning
================
Yue Shi, Ph.D candidate, University of Washington
5/8/2018

-   [Maximum Likelihood Estimation (MLE)](#maximum-likelihood-estimation-mle)
    -   [Example: Halitosis example](#example-halitosis-example)
-   [Maximum a Posteriori (MAP) Estimation](#maximum-a-posteriori-map-estimation)
    -   [beta distribution](#beta-distribution)
    -   [Example:Thumbtack example (unfair coin)](#examplethumbtack-example-unfair-coin)
-   [Linear Regression](#linear-regression)
    -   [Simulation example](#simulation-example)
-   [Model selection and cross validation](#model-selection-and-cross-validation)
    -   [Ridge regression or L2 regularized linear regression](#ridge-regression-or-l2-regularized-linear-regression)
    -   [LASSO regression or L1 regularized linear regression](#lasso-regression-or-l1-regularized-linear-regression)
    -   [LOOCV (leave one out cross validation)](#loocv-leave-one-out-cross-validation)
    -   [Simulation study](#simulation-study)
    -   [Insulin example](#insulin-example)

Maximum Likelihood Estimation (MLE)
-----------------------------------

Likelihood is:
*L*(*θ* : *D**a**t**a*)=*P*(*D**a**t**a*|*θ*)
 Maximum likelihood Estimation is to find *θ*<sup>\*</sup> to maximize the likelihood. This *θ*<sup>\*</sup> is the maximum likelihood estimate of *θ*.

#### Example: Halitosis example

Halitosis, colloquially called bad breath, is a symptom in which a noticeably unpleasant odor is present on the exhaled breath. Halitosis is partly genetically determined. The genotype aa has a 40% chance of getting the disease, and the other two possible genotypes, AA and Aa, each has a 10% chance of getting the disease. We want to estimate the frequency of the A allele.

If the gene frequency of the A allele is p, and that of a is 1-p, then the frequency of the disease is expected to be (if the genotypes are in Hardy-Weinberg proportions as a result of random mating):

*F* = *p*<sup>2</sup> × (0.1)+2*p*(1 − *p*)×(0.1)+(1 − *p*)<sup>2</sup> × (0.4)

Now suppose we observe 1000 individuals and find that the 182 of them have the disease. Using a binomial distribution, the probability that 182 out of 1000 have the disease is the binomial class probability for 182 out of 1000 when the probability of the event is F (which is a function of p). This is

$$\\frac{1000!}{182!~818!} F^{182} (1-F)^{818}$$

``` r
##Given p, calcualte the probability of getting Halitosis
H=function(x){
  hal=0.1*x^2 + 0.1*2*x*(1-x) + 0.4*(1-x)^2 
  return(hal)
}

##Given H, calcualte the log likelihood of getting the data. 
LD=function(x){
  LL=182*log(x)+818*log(1-x)
  return(LL)
}

##Now let's find p which maximize the LD. 
fp=seq(0,1,0.001)
fLD=LD(H(fp))
mle=fp[which.max(fLD)]
mle
```

    ## [1] 0.477

``` r
plot(fLD~fp, xlab="Allele frequency of A", ylab="Data Likelihood")
grid(10,10)
abline(v=mle,col="red")
```

![](estimation_files/figure-markdown_github/unnamed-chunk-1-1.png)

Maximum a Posteriori (MAP) Estimation
-------------------------------------

Posterior probability is:
$$
P(\\theta | D) = \\frac{P(D|\\theta)P(\\theta)}{P(D)} = \\frac{P(D|\\theta)P(\\theta)}{\\int P(D|\\theta)P(\\theta)d\\theta}
$$
 Since the demoninator is not a function of *θ*, therefore we can ignore the denominator， then
$$
\\begin{eqnarray}
P(\\theta | D) &\\propto& P(D|\\theta)P(\\theta)\\\\
Posterior &\\propto& Likelihood \* Prior
\\end{eqnarray}
$$
 MAP estimation is to find *θ* that maximizes posterior *P*(*θ*|*D*).
The difference between MLE and MAP is that:
For MLE: find *θ* that maximizes *l**o**g**P*(*D*|*θ*); whereas
For MAP: find *θ* that maximizes *l**o**g**P*(*D*|*θ*)+*l**o**g**P*(*θ*).

The difference between MAP and Bayesian estimation is that:
MAP ignore the demoninator which integrate all of the possible *θ*, which makes Bayesian method much more computational demanding. MAP is must faster and popular now.

#### beta distribution

**beta distribution**is parameterized by two shape constraints *α* and *β*. It does the job nicely for expressing the prior belief of probability which is restricted in the range \[0,1\].

$$
P(p)=\\frac{1}{B(\\alpha,\\beta)}p^{\\alpha-1}(1-p)^{\\beta-1}
$$
 where $B(,) is the beta function.

When both *α* and *β* are greater than zero, it has the following properties:
$$
mean=\\frac{a}{a+b}
$$

$$
mode=\\frac{\\alpha-1}{\\alpha+\\beta-2}
$$
$$
variance=\\frac{\\alpha\\beta}{(\\alpha+\\beta)^2(\\alpha+\\beta+1)}
$$

#### Example:Thumbtack example (unfair coin)

``` r
nh <- 100 #the number of heads
nt <- 50 #the number of tails
alpha <- 1000 #alpha and beta are hyperparameters
beta <- 1000
logLikelihood <- function(p, nh, nt){
  return(nh*log(p)+nt*log(1-p))
}

logPosterior <- function(p, nh, nt, alpha, beta){
  return((nh+alpha-1)*log(p)+(nt+beta-1)*log(1-p))
}
p=seq(0,1,0.01)
LL=logLikelihood(p,nh,nt)
LP=logPosterior(p,nh,nt,alpha,beta)
par(mfrow=c(1,2))
plot(LL~p, ylab="Likelihood",main="Likelihood")
plot(LP~p, ylab="Likelihood",main="Posterior (alpha=1000, beta=1000)")
```

![](estimation_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
LL.max=p[which.max(LL)]
LL.max
```

    ## [1] 0.67

``` r
nh/(nh+nt)
```

    ## [1] 0.6666667

``` r
LP.max=p[which.max(LP)]
LP.max
```

    ## [1] 0.51

``` r
(nh+alpha-1)/(nh+nt+alpha+beta-1)
```

    ## [1] 0.5114007

Let's minimize the influence of priors

``` r
alpha=30
beta=30
LL=logLikelihood(p,nh,nt)
LP=logPosterior(p,nh,nt,alpha,beta)
par(mfrow=c(1,2))
plot(LL~p, ylab="Likelihood",main="Likelihood")
plot(LP~p, ylab="Likelihood",main="Posterior (alpha=30, beta=30)")
```

![](estimation_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
LL.max=p[which.max(LL)]
LL.max
```

    ## [1] 0.67

``` r
nh/(nh+nt)
```

    ## [1] 0.6666667

``` r
LP.max=p[which.max(LP)]
LP.max
```

    ## [1] 0.62

``` r
(nh+alpha-1)/(nh+nt+alpha+beta-1)
```

    ## [1] 0.6172249

Linear Regression
-----------------

Linear regression is used to predict the value of an outcome variable (continuous) Y based on one or more variables X. The aim is to establish a linear relationship (a mathmatical formula) between the predictor variables and the response variable.
*Y* = *β*<sub>0</sub> + *β*<sub>1</sub>*X*<sub>1</sub> + *β*<sub>2</sub>*X*<sub>2</sub> + *ϵ*
 *ϵ* is the error term, the part of Y the regression model is unble to explain, and it follows *N*(0, *σ*<sup>2</sup>). The goal of linear regression is to find *β* to minimize the sum of squared errors (SSE).
$$
SSE=\\sum\_{i}^{n}(y\_i-\\hat{y\_i})^2
$$

$\\hat{y\_{i}}$ is the fitted value for observation i.
F-statistic is a measure of goodness of fit with associated p value for the model.

#### Simulation example

``` r
X1=rnorm(100,mean=10, sd=5)
X2=rnorm(100,mean=-3, sd=4)
Y=1.1+0.02*X1+0.3*X2+rnorm(100,mean=0, sd=3)
```

Method 1: lsfit function in R

``` r
X=cbind(X1,X2)
rr=lsfit(X,Y)
ls.print(rr)
```

    ## Residual Standard Error=2.6631
    ## R-Square=0.1545
    ## F-statistic (df=2, 97)=8.8643
    ## p-value=3e-04
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   0.8276  0.5969  1.3866   0.1687
    ## X1          0.0133  0.0504  0.2645   0.7919
    ## X2          0.2851  0.0677  4.2084   0.0001

``` r
## quadratic regression
X3=X1^2
Y2=1.1+0.02*X1+0.3*X3+rnorm(100,mean=0,sd=3)
X=cbind(X1,X3)
rr=lsfit(X,Y2)
ls.print(rr)
```

    ## Residual Standard Error=2.7654
    ## R-Square=0.9934
    ## F-statistic (df=2, 97)=7298.914
    ## p-value=0
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   1.7380  0.7991  2.1749   0.0321
    ## X1         -0.2679  0.1610 -1.6638   0.0994
    ## X3          0.3156  0.0078 40.6355   0.0000

Method 2: lm function in R (you don't need to combine X1 and X2)

``` r
fit=lm(Y~X1+X2)
summary(fit)
```

    ## 
    ## Call:
    ## lm(formula = Y ~ X1 + X2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -6.4195 -1.7700 -0.0709  1.8400  5.0961 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.82760    0.59685   1.387    0.169    
    ## X1           0.01332    0.05036   0.265    0.792    
    ## X2           0.28508    0.06774   4.208 5.74e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.663 on 97 degrees of freedom
    ## Multiple R-squared:  0.1545, Adjusted R-squared:  0.1371 
    ## F-statistic: 8.864 on 2 and 97 DF,  p-value: 0.0002913

``` r
rr=lsfit(X,Y)
ls.print(rr)
```

    ## Residual Standard Error=2.861
    ## R-Square=0.0242
    ## F-statistic (df=2, 97)=1.2029
    ## p-value=0.3048
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   1.0564  0.8267  1.2778   0.2044
    ## X1         -0.2504  0.1666 -1.5030   0.1361
    ## X3          0.0124  0.0080  1.5461   0.1253

``` r
fit$coefficients
```

    ## (Intercept)          X1          X2 
    ##  0.82760142  0.01332246  0.28507984

``` r
fit$fitted.values
```

    ##            1            2            3            4            5 
    ## -0.020267535  0.330257196 -0.005981663 -0.840783905 -0.546879687 
    ##            6            7            8            9           10 
    ##  0.112866080  1.127121104  0.295464007  0.296731416  0.394665409 
    ##           11           12           13           14           15 
    ##  0.070348110 -1.817340579  1.278478165 -0.143623832 -0.308462740 
    ##           16           17           18           19           20 
    ## -0.834446982  0.164418935 -0.050498299  0.221401722 -2.629158974 
    ##           21           22           23           24           25 
    ## -0.909815271 -0.257053802 -1.622895100  0.318016936 -0.173081647 
    ##           26           27           28           29           30 
    ##  2.817524266  1.158443673  0.726554833  2.051738983 -1.208066882 
    ##           31           32           33           34           35 
    ##  0.232901922  0.482001650 -0.381504899  0.053550850 -0.553346110 
    ##           36           37           38           39           40 
    ## -0.070375401  1.923809863 -0.085501960  2.280230300  1.885068689 
    ##           41           42           43           44           45 
    ##  1.852793997 -2.145750148 -1.034006041 -0.837950269  0.473230794 
    ##           46           47           48           49           50 
    ##  1.188932781 -1.218475631 -0.461351409 -0.298574257 -1.352671833 
    ##           51           52           53           54           55 
    ##  1.974411330  0.248088061 -0.152710080  1.253359820  2.855897238 
    ##           56           57           58           59           60 
    ## -0.274824782  1.029717173  0.908236398  1.972237483  1.214391648 
    ##           61           62           63           64           65 
    ## -0.371743325  0.863179080  0.443284182 -0.634852739 -0.668130675 
    ##           66           67           68           69           70 
    ##  1.622458637  0.007999123  0.226402126  0.279484323  0.842497134 
    ##           71           72           73           74           75 
    ## -0.238508390 -1.613746979 -0.237816229  0.243464155 -0.522689792 
    ##           76           77           78           79           80 
    ## -1.812646971 -1.779259199  0.742580697 -0.738415732  0.197099428 
    ##           81           82           83           84           85 
    ##  0.629650457  0.888209045 -1.870932026 -0.471798295 -1.065823280 
    ##           86           87           88           89           90 
    ##  0.204973761 -0.021303406 -0.523238581 -0.007755391  0.047922420 
    ##           91           92           93           94           95 
    ## -0.629046015 -1.487939884  1.670161676  1.826493928  2.449756932 
    ##           96           97           98           99          100 
    ## -1.099893052  1.233285420  0.617454196  1.847935385  0.527921707

``` r
fit$residuals
```

    ##            1            2            3            4            5 
    ## -4.857908514 -1.882323415  1.547277625 -1.234414830 -1.125401670 
    ##            6            7            8            9           10 
    ##  1.068239692 -3.696430813  2.264661041 -0.185153924  2.003082706 
    ##           11           12           13           14           15 
    ## -1.483266451 -2.655215617  4.997542180 -1.879513622  1.439758051 
    ##           16           17           18           19           20 
    ## -0.055506030  3.766437293  0.485269722 -0.349873568 -2.681227733 
    ##           21           22           23           24           25 
    ## -0.476717372 -3.106622214 -2.502041237 -2.717781582  1.090428921 
    ##           26           27           28           29           30 
    ##  0.859449296 -1.204235950 -1.393542713 -6.419489575  3.965602372 
    ##           31           32           33           34           35 
    ##  2.142220795  1.785601165  0.006764883 -3.946864472  1.538517863 
    ##           36           37           38           39           40 
    ##  0.472326650  0.199020562  1.312067257  1.532978892 -0.086365982 
    ##           41           42           43           44           45 
    ## -1.183909331  0.677190085  1.101751802  5.008689566  2.887746436 
    ##           46           47           48           49           50 
    ## -3.808226127 -0.425788760 -6.301072551  4.965040121  2.856863305 
    ##           51           52           53           54           55 
    ## -0.873594625 -3.453810283  0.454563633 -2.427591544  2.619956698 
    ##           56           57           58           59           60 
    ##  4.649040205 -1.266955696 -1.038493652 -1.413434492  3.904666034 
    ##           61           62           63           64           65 
    ## -0.954321092 -2.542974163  0.458073654  2.037926369  2.343771631 
    ##           66           67           68           69           70 
    ##  3.671280758 -1.578102456 -5.094546150  5.096112130 -2.437656064 
    ##           71           72           73           74           75 
    ##  2.525928233  0.629030680 -2.691183291  0.794415384  0.201905530 
    ##           76           77           78           79           80 
    ## -2.114805832  5.031102694  2.690137328  2.267918785  1.624361981 
    ##           81           82           83           84           85 
    ## -0.549540350 -1.733459889 -3.411266629 -1.471610425  1.230201190 
    ##           86           87           88           89           90 
    ## -0.523149464 -1.447457494  2.743855114 -0.124977302  0.726204090 
    ##           91           92           93           94           95 
    ##  0.778386499 -1.693367524  3.766779794 -2.415384605  4.015162081 
    ##           96           97           98           99          100 
    ## -1.479480287  3.460840907 -3.667172987 -3.366211900 -2.266707438

``` r
anova(fit)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Y
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## X1         1   0.12   0.125  0.0176    0.8947    
    ## X2         1 125.61 125.610 17.7110 5.745e-05 ***
    ## Residuals 97 687.94   7.092                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(fit$residuals~fit$fitted.values)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-6-1.png)

Model selection and cross validation
------------------------------------

Linear regression has a few problems, especially when there are more variables than the number of samples. It will tend to overfit the data. Therefore, we need more constraints on the number of variables. Then the problem is to find *β* that minimize the cost function.

### Ridge regression or L2 regularized linear regression

The cost function is defined as:
$$
f(\\beta\_0, \\beta\_1, \\cdot\\cdot\\cdot, \\beta\_p)=\\sum\_{i=1}^{n}\[y\_i-(\\beta\_0+\\beta\_1x\_{1i}+\\cdot\\cdot\\cdot+\\beta\_px\_{pi})\]^2+\\lambda\\sum\_{j=1}^p\\beta\_j^2
$$
 Note that L2 regularization term (sum of squared *β*) encourages to choose *β* that have small manitude. *λ* is the tuning parameter, the higher the *λ*, the higher the penalty.

### LASSO regression or L1 regularized linear regression

The cost function is defined as:
$$
f(\\beta\_0, \\beta\_1, \\cdot\\cdot\\cdot, \\beta\_p)=\\sum\_{i=1}^{n}\[y\_i-(\\beta\_0+\\beta\_1x\_{1i}+\\cdot\\cdot\\cdot+\\beta\_px\_{pi})\]^2+\\lambda\\sum\_{j=1}^p|\\beta\_j|
$$
 Note that L1 regularization term (sum of absolute *β*) encourages many *β* values to be set to zero. Same as L2 regularization, *λ* is the tunning parameter. L1 regularization will give more extreme *β* values, since it will bring small *β* values to zero and keeps the large values.

### LOOCV (leave one out cross validation)

You want to find the model with the smallest test error. What is test error?
$$
f(\\beta\_0, \\beta\_1, \\cdot\\cdot\\cdot, \\beta\_p)=\\sum\_{i=1}^{n}\[y\_i-(\\beta\_0+\\beta\_1x\_{1i}+\\cdot\\cdot\\cdot+\\beta\_px\_{pi})\]^2
$$
 In order to get the true test error, you can use cross validation, which partition the data into K equally (or nearly equally) sized segments or folds. Then perform k iterations of training and validation such that each iteration a different fold of the data is held-out for validation while the remaining k-1 folds are used for learning. The true test erros if the average test error of the k-fold iterations. In data mining and machine learning, 10-fold cross-validation (k=10) is the most common. LOOCV is a special case of k-fold cross-validation where k equals the number of samples in the data. The accuracy estimate obtained using LOOCV is known to be almost unbiased. It is widely used when the available data are very rare.

### Simulation study

``` r
## Generate 3 variables (or features X1, X2, X3)
set.seed(1)
X1=rnorm(50,mean=10, sd=5)
X2=rnorm(50,mean=-3, sd=4)
X3=rnorm(50,mean=2, sd=2)
X=cbind(X1,X2,X3)

## Generate observation values
Y <- 1.1 + 0.02*X1 + 0.3*X2 + 0.5*X3 + rnorm(50,mean=0,sd=3) 
```

Let's try ridge regression (L2).

``` r
library(MASS)
ridge=lm.ridge(Y~X1+X2+X3, lambda=seq(0,20,1)) 
## When there are a lot of features, it will be easier to do lm.ridge(Y~X)
ridge$GCV ## get the generalized cross validation test errors (same as LOOCV). 
```

    ##         0         1         2         3         4         5         6 
    ## 0.1776816 0.1772726 0.1769474 0.1766973 0.1765142 0.1763910 0.1763213 
    ##         7         8         9        10        11        12        13 
    ## 0.1762997 0.1763210 0.1763806 0.1764747 0.1765995 0.1767518 0.1769287 
    ##        14        15        16        17        18        19        20 
    ## 0.1771274 0.1773458 0.1775815 0.1778327 0.1780976 0.1783747 0.1786625

We want to find *λ* to give the smallest GCV. If you couldn't get a "well" curve, increase the lambda range.

``` r
plot(seq(0,20,1),ridge$GCV)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
which.min(ridge$GCV)
```

    ##  7 
    ##  8

When *λ* = 7 , the model has the smallest GCV. What about the coefficients with *λ* = 7?

``` r
ridge$coef[,8]
```

    ##         X1         X2         X3 
    ## -0.6610979  0.9197593  1.3680452

``` r
##Compare with the real values when we generated Y. 
```

Get the L2 regularization term (the penalty term)

``` r
colSums((ridge$coef)^2)
```

    ##        0        1        2        3        4        5        6        7 
    ## 4.097200 3.938470 3.788791 3.647488 3.513947 3.387610 3.267968 3.154555 
    ##        8        9       10       11       12       13       14       15 
    ## 3.046947 2.944754 2.847619 2.755213 2.667234 2.583403 2.503464 2.427180 
    ##       16       17       18       19       20 
    ## 2.354331 2.284713 2.218139 2.154434 2.093435

``` r
plot(seq(0,20,1),colSums((ridge$coef)^2))
```

![](estimation_files/figure-markdown_github/unnamed-chunk-11-1.png)

L2 regularization term decreases as lambda value increases.

### Insulin example

``` r
a <- read.table(header = T, file="http://www.cs.washington.edu/homes/suinlee/genome560/data/mice.txt")  

X <- as.matrix(a[,2:7])    # X is sex, weight, length, Triglyceride, Total Cholesterol, FFA
                           # as.matrix is for 
Y <- a[,8]                 # Y is the insulin level

all=cbind(X,Y)
pairs(all) ## pairwise correlation
```

![](estimation_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
ridge_lm=lm.ridge(Y~X,lambda=seq(0,20,1))
plot(seq(0,20,1),ridge_lm$GCV)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
which.min(ridge_lm$GCV)
```

    ##  9 
    ## 10

``` r
ridge_lm$coef[,10]
```

    ##        Xsex   Xweight_g  Xlength_cm     XTrigly XTotal_Chol        XFFA 
    ## -0.11701872  0.23404427  0.01042570  0.03482851  0.03578169  0.07695456
