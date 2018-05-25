Parameter Estimation and Machine Learning
================
Yue Shi, Ph.D candidate, University of Washington
5/8/2018

-   [Maximum Likelihood Estimation (MLE)](#maximum-likelihood-estimation-mle)
    -   [Example: Halitosis example](#example-halitosis-example)
    -   [Example: model selection based on maximum likelihood](#example-model-selection-based-on-maximum-likelihood)
-   [Maximum a Posteriori (MAP) Estimation](#maximum-a-posteriori-map-estimation)
    -   [beta distribution](#beta-distribution)
    -   [Example:Thumbtack example (unfair coin)](#examplethumbtack-example-unfair-coin)
-   [Ordinary Linear Regression](#ordinary-linear-regression)
    -   [Simulation example](#simulation-example)
    -   [Example: quantitative trait loci analysis for cholesterol levels](#example-quantitative-trait-loci-analysis-for-cholesterol-levels)
-   [Model selection and cross validation](#model-selection-and-cross-validation)
    -   [Ridge regression or L2 regularized linear regression](#ridge-regression-or-l2-regularized-linear-regression)
    -   [LASSO regression or L1 regularized linear regression](#lasso-regression-or-l1-regularized-linear-regression)
    -   [LOOCV (leave one out cross validation)](#loocv-leave-one-out-cross-validation)
    -   [Simulation study](#simulation-study)
    -   [Example: Insulin](#example-insulin)
    -   [Example: Cholesterol](#example-cholesterol)

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

#### Example: model selection based on maximum likelihood

we will implement an algorithm for selecting among various structures of the regulatory network. Specifically, we will focus on two possible models of the galactose regulatory network in S. cerevisiae. We will select a model based on the expression data on these three genes measured across S. cerevisiae individuals.

Likelihood function of Model 1:
*L*(*D*|*θ*)=*P*(*G**a**l*80)*P*(*G**a**l*4|*G**a**l*80)*P*(*G**a**l*2|*G**a**l*4)
 Likelihood function of Model 2:
*L*(*D*|*θ*)=*P*(*G**a**l*4)*P*(*G**a**l*80)*P*(*G**a**l*2|*G**a**l*4, *G**a**l*80)
 Import the data

``` r
a <- read.table(header = T, file="https://sites.google.com/a/cs.washington.edu/genome560-spr18/disc-gal80-gal4-gal2.txt") 
head(a)
```

    ##     EXP X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19
    ## 1 Gal80  1  1  1  0  1  0  0  0  0   0   0   1   0   0   1   0   0   0   1
    ## 2  Gal4  1  0  0  1  1  1  1  1  1   1   1   0   1   0   0   0   1   1   0
    ## 3  Gal2  1  1  1  1  1  1  1  1  1   1   1   1   1   0   0   0   1   0   0
    ##   X20 X21 X22 X23 X24 X25 X26 X27 X28 X29 X30 X31 X32 X33 X34 X35 X36 X37
    ## 1   0   1   0   0   0   1   0   1   1   1   0   0   0   0   0   0   1   0
    ## 2   1   1   1   1   1   0   0   1   0   0   1   1   1   0   1   1   0   0
    ## 3   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   0   0
    ##   X38 X39 X40 X41 X42 X43 X44 X45 X46 X47 X48 X49 X50 X51 X52 X53 X54 X55
    ## 1   0   1   0   0   0   0   0   0   0   1   0   1   1   1   1   1   0   1
    ## 2   1   0   1   1   0   1   1   1   1   0   1   0   0   0   1   0   0   1
    ## 3   1   1   1   1   0   1   1   1   1   0   0   0   0   0   0   0   0   1
    ##   X56 X57 X58 X59 X60 X61 X62 X63 X64 X65 X66 X67 X68 X69 X70 X71 X72 X73
    ## 1   0   0   0   0   0   1   1   0   0   0   1   1   0   1   1   0   1   1
    ## 2   1   1   0   1   0   0   1   1   1   0   1   0   0   0   1   0   0   1
    ## 3   1   1   1   1   0   1   1   0   1   0   0   0   1   0   0   1   0   1
    ##   X74 X75 X76 X77 X78 X79 X80 X81 X82 X83 X84 X85 X86 X87 X88 X89 X90 X91
    ## 1   1   1   1   1   0   0   1   0   1   1   0   1   1   0   1   0   0   1
    ## 2   1   0   0   0   1   1   0   0   0   1   1   0   0   0   1   1   1   1
    ## 3   1   0   1   0   0   1   0   0   0   0   0   0   1   0   0   1   1   0
    ##   X92 X93 X94 X95 X96 X97 X98 X99 X100 X101 X102 X103 X104 X105 X106 X107
    ## 1   1   1   0   0   1   1   1   1    1    0    0    1    1    1    1    0
    ## 2   1   1   0   1   0   0   0   1    1    0    1    0    0    0    0    1
    ## 3   1   0   0   1   1   0   1   0    1    0    0    0    0    1    1    1
    ##   X108 X109 X110 X111 X112
    ## 1    1    1    0    0    1
    ## 2    0    0    0    0    0
    ## 3    0    1    0    1    1

``` r
dim(a)
```

    ## [1]   3 113

``` r
a=as.data.frame(a)
```

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

![](estimation_files/figure-markdown_github/unnamed-chunk-4-1.png)

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

Ordinary Linear Regression
--------------------------

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

    ## Residual Standard Error=2.7104
    ## R-Square=0.2968
    ## F-statistic (df=2, 97)=20.4701
    ## p-value=0
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   1.1104  0.6962  1.5948   0.1140
    ## X1          0.0649  0.0624  1.0405   0.3007
    ## X2          0.4240  0.0664  6.3833   0.0000

``` r
## quadratic regression
X3=X1^2
Y2=1.1+0.02*X1+0.3*X3+rnorm(100,mean=0,sd=3)
X=cbind(X1,X3)
rr=lsfit(X,Y2)
ls.print(rr)
```

    ## Residual Standard Error=2.5857
    ## R-Square=0.9924
    ## F-statistic (df=2, 97)=6354.503
    ## p-value=0
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   0.4274  1.0400  0.4110   0.6820
    ## X1          0.0956  0.1984  0.4821   0.6308
    ## X3          0.2971  0.0092 32.4603   0.0000

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
    ## -7.5524 -2.0715 -0.2669  1.7833  6.3085 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.11038    0.69623   1.595    0.114    
    ## X1           0.06493    0.06240   1.040    0.301    
    ## X2           0.42402    0.06643   6.383 5.96e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.71 on 97 degrees of freedom
    ## Multiple R-squared:  0.2968, Adjusted R-squared:  0.2823 
    ## F-statistic: 20.47 on 2 and 97 DF,  p-value: 3.832e-08

``` r
rr=lsfit(X,Y)
ls.print(rr)
```

    ## Residual Standard Error=3.2281
    ## R-Square=0.0026
    ## F-statistic (df=2, 97)=0.1243
    ## p-value=0.8833
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   0.7820  1.2984  0.6022   0.5484
    ## X1          0.0527  0.2477  0.2126   0.8321
    ## X3         -0.0038  0.0114 -0.3351   0.7383

``` r
fit$coefficients
```

    ## (Intercept)          X1          X2 
    ##  1.11038061  0.06492612  0.42402272

``` r
fit$fitted.values
```

    ##            1            2            3            4            5 
    ##  0.675939000 -1.218431892  0.630026561  2.637190290  1.665659875 
    ##            6            7            8            9           10 
    ## -0.660759198  2.219068374  0.865326162  1.372706891  1.557933969 
    ##           11           12           13           14           15 
    ##  0.472073831  5.266000396  1.035851982 -0.602412758  1.702279136 
    ##           16           17           18           19           20 
    ##  4.731268342  1.353433380  0.863510252  1.337789747 -0.006479785 
    ##           21           22           23           24           25 
    ## -0.722341226 -1.749298996  2.369887159 -0.372951316  1.296716791 
    ##           26           27           28           29           30 
    ## -0.779464667  0.282471433 -0.353433908  0.273445625 -0.756325353 
    ##           31           32           33           34           35 
    ##  0.897222364 -0.528437879 -1.346944795  2.030968300  2.503966491 
    ##           36           37           38           39           40 
    ##  1.385708358 -2.204867385  0.201494544 -3.650773117  1.495890144 
    ##           41           42           43           44           45 
    ## -2.428036810 -0.077329110 -1.227960891  0.247436037  0.176046530 
    ##           46           47           48           49           50 
    ##  3.786718224  0.205028798  4.032430209 -0.246740548  1.052999484 
    ##           51           52           53           54           55 
    ##  0.230861249  2.861204914  3.087686041  1.667373882  0.634836596 
    ##           56           57           58           59           60 
    ##  0.520180673 -0.013389458 -0.359357935 -2.886798372  1.904823535 
    ##           61           62           63           64           65 
    ## -1.692195712 -0.892727392 -2.441439932 -0.037325462  3.902807962 
    ##           66           67           68           69           70 
    ##  0.997728468 -0.766240689  0.338084925  0.791105909  2.763943738 
    ##           71           72           73           74           75 
    ##  0.979640683  3.159115829  3.016844899  1.111728195  2.463801315 
    ##           76           77           78           79           80 
    ##  1.819822709  0.996404569  4.446262836  3.373110252 -0.780151474 
    ##           81           82           83           84           85 
    ##  1.636995263  3.633597443  2.583664913 -2.122656658  0.884929541 
    ##           86           87           88           89           90 
    ##  0.168248751  0.129256607  2.954204474  0.726279813 -0.242500076 
    ##           91           92           93           94           95 
    ## -2.010474431  0.011159729  3.599128843  1.079299401  2.239993112 
    ##           96           97           98           99          100 
    ##  0.917166356  2.574870300  0.861092187 -0.168130373  1.131609614

``` r
fit$residuals
```

    ##          1          2          3          4          5          6 
    ## -0.1227163  2.2105757 -0.3668267 -2.0707059 -2.5771797  4.7762153 
    ##          7          8          9         10         11         12 
    ##  1.1290842 -0.8611271 -2.3999196  1.5330244  1.2765874 -1.1861641 
    ##         13         14         15         16         17         18 
    ##  1.8968649  5.4782977  5.9334310  2.9929938  1.9372422 -2.0740660 
    ##         19         20         21         22         23         24 
    ## -0.6654859  1.8814830  2.1009534  0.8830376  2.6886955  0.7472612 
    ##         25         26         27         28         29         30 
    ## -4.1143561 -2.0302047 -0.1669608 -3.1223606 -1.1038898 -2.3915177 
    ##         31         32         33         34         35         36 
    ## -7.5524002 -2.2607620  0.7590359 -3.0208028  6.3085012  1.7291649 
    ##         37         38         39         40         41         42 
    ## -0.8885444 -1.9755237  0.7399349 -3.0360516 -0.8128978  3.7847215 
    ##         43         44         45         46         47         48 
    ## -3.3627290  5.9891374  2.3219513 -3.0613489  0.8173157 -1.3899545 
    ##         49         50         51         52         53         54 
    ##  2.5002218  0.0739791 -2.8318932 -0.9127084 -2.1761638  3.9830984 
    ##         55         56         57         58         59         60 
    ## -1.1164758 -1.1274848  1.5536022 -1.7241155 -3.4627406 -1.1081483 
    ##         61         62         63         64         65         66 
    ## -0.6475314 -4.0991172  0.1625953  1.4425727 -1.3644095  2.1194403 
    ##         67         68         69         70         71         72 
    ##  0.8503586  1.3782278  1.6023139  0.9394864  1.9071155  4.8362930 
    ##         73         74         75         76         77         78 
    ## -3.8600115 -1.2167093 -3.0396970 -3.3928540  0.2617827 -0.5144449 
    ##         79         80         81         82         83         84 
    ## -3.6181306 -2.1783475  1.3472899  5.3773622 -1.1120397  1.9445107 
    ##         85         86         87         88         89         90 
    ## -2.2896053  1.4328731 -0.7674490  2.8694063 -2.4497877  0.8689489 
    ##         91         92         93         94         95         96 
    ##  0.4164216  2.3252374  2.6945473  1.7505347 -1.4643661  5.0741099 
    ##         97         98         99        100 
    ## -2.5324201 -5.2717330 -1.6874814 -1.0474785

``` r
anova(fit)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Y
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## X1         1   1.42   1.420  0.1934    0.6611    
    ## X2         1 299.34 299.344 40.7468 5.957e-09 ***
    ## Residuals 97 712.60   7.346                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(fit$residuals~fit$fitted.values)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-7-1.png)

#### Example: quantitative trait loci analysis for cholesterol levels

We are given the genotype and phenotype data from 334 mouse individuals. The genotype data measure binary genotype values of 1333 genetic markers for each mouse, and the phenotype data measure the normalized blood cholesterol levels. Given these data, we want to find the quantitative trait loci (QTLs) that contribute to elevated cholesterol level.

``` r
pheno=read.table("phenotype.txt", header=T)
geno=read.table("genotype.txt", header=T)
pheno=pheno[,-1]
Y=t(pheno)
colnames(Y)="chol"
dim(pheno) ## sample size: 334
```

    ## [1]   1 334

``` r
## we need reshape the geno data in the similar fashion as pheno.
X=geno[,-1]
names=geno[,1]
X=as.data.frame(t(X))
colnames(X)=names
dim(X) ## 334 samples,1333 loci
```

    ## [1]  334 1333

Let's do simple linear regression for each genetic marker and compute the SSE (sum of squared error) for each marker. Find out which marker has the minimum SSE. Remember, residual means the difference between the fitted value and the true value.

``` r
sse=vector()
r.sq=vector()
for(i in 1:ncol(X)){
  lm=lm(Y~X[,i])
  sse[i]=sum((lm$residuals)^2)
  r.sq[i]=summary(lm)$r.squared
}
length(sse)
```

    ## [1] 1333

``` r
length(r.sq)
```

    ## [1] 1333

``` r
which.min(sse)
```

    ## [1] 847

``` r
which.max(sse)
```

    ## [1] 598

``` r
which.min(r.sq)
```

    ## [1] 598

``` r
which.max(r.sq)
```

    ## [1] 847

``` r
max(r.sq)
```

    ## [1] 0.1239932

``` r
plot(r.sq)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
lm.598=lm(Y~X[,598])
lm.847=lm(Y~X[,847])

par(mfrow=c(1,2))
plot(Y~X[,598], main="Marker 598")
abline(lm.598)
plot(Y~X[,847], main="Marker 847")
abline(lm.847)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-9-2.png)

It turns out Marker 847 has the minimum SSE, Marker 598 has the maximum SSE. It is confirmed by *R*<sup>2</sup>. Marker 847 has the maximum *R*<sup>2</sup>, whereas Marker 598 has the minimum *R*<sup>2</sup>.

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

![](estimation_files/figure-markdown_github/unnamed-chunk-12-1.png)

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

![](estimation_files/figure-markdown_github/unnamed-chunk-14-1.png)

L2 regularization term decreases as lambda value increases.

### Example: Insulin

``` r
a <- read.table(header = T, file="http://www.cs.washington.edu/homes/suinlee/genome560/data/mice.txt")  

X <- as.matrix(a[,2:7])    # X is sex, weight, length, Triglyceride, Total Cholesterol, FFA
                           # as.matrix is for 
Y <- a[,8]                 # Y is the insulin level

all=cbind(X,Y)
pairs(all) ## pairwise correlation
```

![](estimation_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
ridge_lm=lm.ridge(Y~X,lambda=seq(0,20,1))
plot(seq(0,20,1),ridge_lm$GCV)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-16-1.png)

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

### Example: Cholesterol

Use the same data as above: quantitative trait loci analysis for cholesterol levels, use *R*<sup>2</sup> (coefficient of determination) to compare the ridge regression and the ordinary linear regression.

$$
R^2=1-\\frac{SSE}{SST}
$$
 SSE is the sum of squared error, whereas SST is the total sum of squares
$$
SSE=\\sum\_{i=1}^{n}(y\_i-\\hat{y\_i})^2
$$
$$
SST=\\sum\_{i=1}^{n}(y\_i-\\bar{y})^2
$$
 Note that $\\bar{y}$ is the mean of the observed data.

``` r
pheno=read.table("phenotype.txt", header=T)
geno=read.table("genotype.txt", header=T)
pheno=pheno[,-1]
Y=t(pheno)
colnames(Y)="chol"

X=geno[,-1]
names=geno[,1]
X=as.data.frame(t(X))
colnames(X)=names
X=as.matrix(X)
library(MASS)
lambdas=10^seq(-6,6,0.1)
length(lambdas)
```

    ## [1] 121

``` r
ridge_chol=lm.ridge(Y~X, lambda=lambdas) ##lm.ridge require X to be a matrix. 
select(ridge_chol)
```

    ## modified HKB estimator is -9.0968e-27 
    ## modified L-W estimator is -259.5528 
    ## smallest value of GCV  at 3162.278

``` r
lambda.best=which.min(ridge_chol$GCV)
lambda.best
```

    ## 3.162278e+03 
    ##           96

``` r
chol.best=lm.ridge(Y~X,lambda=lambda.best)

## get fitted values for lm.rige
## Note: there is indeed no predict method for ridgelm object. You will have to do it by hands. 
y.pred=as.matrix(cbind(const=1,X)) %*% coef(chol.best)
SSE=sum((Y-y.pred)^2)
SST=sum((Y-chol.best$ym)^2)
r.ridge=1-(SSE/SST)
r.ridge
```

    ## [1] 0.7972105

When we used only one marker and did single linear regression, Marker 847 gave us the maximum *R*<sup>2</sup> of 0.124. With ridge regression, we got *R*<sup>2</sup> of 0.797. Much better!
