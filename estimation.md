Parameter Estimation and Machine Learning
================
Yue Shi, Ph.D candidate, University of Washington
5/8/2018

-   [Maximum Likelihood Estimation (MLE)](#maximum-likelihood-estimation-mle)
    -   [Example: Halitosis example](#example-halitosis-example)
    -   [Example: Model selection based on maximum likelihood](#example-model-selection-based-on-maximum-likelihood)
-   [Maximum a Posteriori Estimation (MAP)](#maximum-a-posteriori-estimation-map)
    -   [beta distribution](#beta-distribution)
    -   [Example:Thumbtack example (unfair coin)](#examplethumbtack-example-unfair-coin)
-   [Ordinary Linear Regression](#ordinary-linear-regression)
    -   [Example: Simulation data](#example-simulation-data)
    -   [Example: Quantitative trait loci analysis for cholesterol levels](#example-quantitative-trait-loci-analysis-for-cholesterol-levels)
-   [Model selection and cross validation](#model-selection-and-cross-validation)
    -   [Ridge regression or L2 regularized linear regression](#ridge-regression-or-l2-regularized-linear-regression)
    -   [LASSO regression or L1 regularized linear regression](#lasso-regression-or-l1-regularized-linear-regression)
    -   [LOOCV (leave one out cross validation)](#loocv-leave-one-out-cross-validation)
    -   [Example: Simulation data](#example-simulation-data-1)
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

#### Example: Model selection based on maximum likelihood

we will implement an algorithm for selecting among various structures of the regulatory network. Specifically, we will focus on two possible models of the galactose regulatory network in S. cerevisiae. We will select a model based on the expression data on these three genes measured across S. cerevisiae individuals.

Likelihood function of Model 1:
*L*(*D*|*θ*)=*P*(*G**a**l*80)*P*(*G**a**l*4|*G**a**l*80)*P*(*G**a**l*2|*G**a**l*4)
 Likelihood function of Model 2:
*L*(*D*|*θ*)=*P*(*G**a**l*4)*P*(*G**a**l*80)*P*(*G**a**l*2|*G**a**l*4, *G**a**l*80)
 Import the data

``` r
a <- read.table(header = T, file="https://sites.google.com/a/cs.washington.edu/genome560-spr18/disc-gal80-gal4-gal2.txt") 
names=a[,1]
a=a[,-1]
rownames(a)=names
dim(a) ## 3 genes, 112 samples
```

    ## [1]   3 112

``` r
a=t(a)
```

Remember: MLE solution is the same as the sample mean. We can calculate each component of the likelihood function from the data, and compare the log likelihood for Model 1 and Model 2.

Calculate *P*(*G**a**l*80) (the probablity of Gal 80 has high expression level) and *P*(*G**a**l*4) (the probablity of Gal 4 has high expression level).

``` r
p.80=sum(a[,1])/nrow(a)
p.4=sum(a[,2])/nrow(a)
```

Calculate *P*(*G**a**l*4|*G**a**l*80), *P*(*G**a**l*2|*G**a**l*4) and *P*(*G**a**l*2|*G**a**l*4, *G**a**l*80)

``` r
p.4H80H=sum(a[,1]==1 & a[,2]==1)/sum(a[,1])
p.4L80H=sum(a[,1]==1 & a[,2]==0)/sum(a[,1])
p.4H80L=sum(a[,1]==0 & a[,2]==1)/(nrow(a)-sum(a[,1]))
p.4L80L=sum(a[,1]==0 & a[,2]==0)/(nrow(a)-sum(a[,1]))


p.2H4H=sum(a[,2]==1 & a[,3]==1)/sum(a[,2])
p.2L4H=sum(a[,2]==1 & a[,3]==0)/sum(a[,2])
p.2H4L=sum(a[,2]==0 & a[,3]==1)/(nrow(a)-sum(a[,2]))
p.2L4L=sum(a[,2]==0 & a[,3]==0)/(nrow(a)-sum(a[,2]))

p.2H80H4H=sum(a[,1]==1 & a[,2]==1 & a[,3]==1)/sum(a[,1]==1 & a[,2]==1)
p.2L80H4H=sum(a[,1]==1 & a[,2]==1 & a[,3]==0)/sum(a[,1]==1 & a[,2]==1)
p.2H80H4L=sum(a[,1]==1 & a[,2]==0 & a[,3]==1)/sum(a[,1]==1 & a[,2]==0)
p.2L80H4L=sum(a[,1]==1 & a[,2]==0 & a[,3]==0)/sum(a[,1]==1 & a[,2]==0)
p.2H80L4H=sum(a[,1]==0 & a[,2]==1 & a[,3]==1)/sum(a[,1]==0 & a[,2]==1)
p.2L80L4H=sum(a[,1]==0 & a[,2]==1 & a[,3]==0)/sum(a[,1]==0 & a[,2]==1)
p.2H80L4L=sum(a[,1]==0 & a[,2]==0 & a[,3]==1)/sum(a[,1]==0 & a[,2]==0)
p.2L80L4L=sum(a[,1]==0 & a[,2]==0 & a[,3]==0)/sum(a[,1]==0 & a[,2]==0)
```

``` r
LL.m1=sum(a[,1])*log(p.80)+
(nrow(a)-sum(a[,1]))*(1-log(p.80))+
sum(a[,1]==1 & a[,2]==1)*log(p.4H80H)+
sum(a[,1]==1 & a[,2]==0)*log(p.4L80H)+
sum(a[,1]==0 & a[,2]==1)*log(p.4H80L)+
sum(a[,1]==0 & a[,2]==0)*log(p.4L80L)+
sum(a[,2]==1 & a[,3]==1)*log(p.2H4H)+
sum(a[,2]==1 & a[,3]==0)*log(p.2L4H)+
sum(a[,2]==0 & a[,3]==1)*log(p.2H4L)+
sum(a[,2]==0 & a[,3]==0)*log(p.2L4L)

LL.m2=sum(a[,1])*log(p.80)+
(nrow(a)-sum(a[,1]))*(1-log(p.80))+
sum(a[,2])*log(p.4)+
(nrow(a)-sum(a[,2]))*(1-log(p.4))+
sum(a[,1]==1 & a[,2]==1 & a[,3]==1)*log(p.2H80H4H)+
sum(a[,1]==1 & a[,2]==1 & a[,3]==0)*log(p.2L80H4H)+
sum(a[,1]==1 & a[,2]==0 & a[,3]==1)*log(p.2H80H4L)+
sum(a[,1]==1 & a[,2]==0 & a[,3]==0)*log(p.2L80H4L)+
sum(a[,1]==0 & a[,2]==1 & a[,3]==1)*log(p.2H80L4H)+
sum(a[,1]==0 & a[,2]==1 & a[,3]==0)*log(p.2L80L4H)+
sum(a[,1]==0 & a[,2]==0 & a[,3]==1)*log(p.2H80L4L)+
sum(a[,1]==0 & a[,2]==0 & a[,3]==0)*log(p.2L80L4L)

LL.m1
```

    ## [1] -80.05638

``` r
LL.m2
```

    ## [1] 44.27889

Since *L*(*M**o**d**e**l*2 : *D**a**t**a*) is greater than *L*(*M**o**d**e**l*1 : *D**a**t**a*), therefore Model 2 is favored over Model 1.

Maximum a Posteriori Estimation (MAP)
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

![](estimation_files/figure-markdown_github/unnamed-chunk-6-1.png)

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

![](estimation_files/figure-markdown_github/unnamed-chunk-7-1.png)

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

#### Example: Simulation data

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

    ## Residual Standard Error=2.6967
    ## R-Square=0.3078
    ## F-statistic (df=2, 97)=21.5651
    ## p-value=0
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   0.0864  0.6895  0.1253   0.9006
    ## X1          0.1027  0.0588  1.7479   0.0836
    ## X2          0.3767  0.0594  6.3454   0.0000

``` r
## quadratic regression
X3=X1^2
Y2=1.1+0.02*X1+0.3*X3+rnorm(100,mean=0,sd=3)
X=cbind(X1,X3)
rr=lsfit(X,Y2)
ls.print(rr)
```

    ## Residual Standard Error=3.0624
    ## R-Square=0.9908
    ## F-statistic (df=2, 97)=5202.706
    ## p-value=0
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   2.9980  1.0986  2.7290   0.0075
    ## X1         -0.3681  0.1967 -1.8713   0.0643
    ## X3          0.3181  0.0087 36.3679   0.0000

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
    ## -6.7393 -1.5400 -0.1868  1.5868  6.6911 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.08638    0.68950   0.125   0.9006    
    ## X1           0.10274    0.05878   1.748   0.0836 .  
    ## X2           0.37673    0.05937   6.345 7.09e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.697 on 97 degrees of freedom
    ## Multiple R-squared:  0.3078, Adjusted R-squared:  0.2935 
    ## F-statistic: 21.57 on 2 and 97 DF,  p-value: 1.785e-08

``` r
rr=lsfit(X,Y)
ls.print(rr)
```

    ## Residual Standard Error=3.2013
    ## R-Square=0.0245
    ## F-statistic (df=2, 97)=1.2172
    ## p-value=0.3005
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept  -1.5142  1.1485 -1.3185   0.1904
    ## X1          0.2219  0.2056  1.0790   0.2832
    ## X3         -0.0058  0.0091 -0.6326   0.5285

``` r
fit$coefficients
```

    ## (Intercept)          X1          X2 
    ##  0.08638082  0.10273825  0.37672948

``` r
fit$fitted.values
```

    ##            1            2            3            4            5 
    ## -0.754040710  0.852239797  0.059455013  1.286553989 -0.375907220 
    ##            6            7            8            9           10 
    ##  0.356424336  0.063511603  1.129013061 -1.140687853 -0.511118728 
    ##           11           12           13           14           15 
    ## -0.402449377  3.968900166  0.965253152 -0.354707118  0.154780860 
    ##           16           17           18           19           20 
    ## -1.999824207  0.376737648  0.510347412 -0.065311979  1.712147897 
    ##           21           22           23           24           25 
    ## -1.772492576 -2.115091558  1.242370215  0.762659736 -0.571674964 
    ##           26           27           28           29           30 
    ##  1.529818157 -1.948452569  1.992508524 -3.382857996 -2.035764814 
    ##           31           32           33           34           35 
    ## -1.804414178 -1.890056106  1.832969490  3.269942951 -0.425813310 
    ##           36           37           38           39           40 
    ## -2.034588282  1.255799400  1.527681216  1.595798601 -0.598150057 
    ##           41           42           43           44           45 
    ## -2.154553118  1.104945645 -2.904832576 -0.769896845 -2.873032100 
    ##           46           47           48           49           50 
    ## -0.190840640  1.557686150 -0.447395837 -1.161298991 -1.406109506 
    ##           51           52           53           54           55 
    ##  3.041911990 -0.239968233  2.108762547  3.956558917  1.752702150 
    ##           56           57           58           59           60 
    ##  0.540612024  2.680690024 -2.321981453 -2.870361718 -2.312004199 
    ##           61           62           63           64           65 
    ##  1.477298102 -1.005901304  3.091157960 -1.555526721 -3.177614794 
    ##           66           67           68           69           70 
    ## -0.334208531 -4.162235670  1.717371130 -0.002983777 -0.495004091 
    ##           71           72           73           74           75 
    ##  2.208840657 -2.278063811  0.534164311 -0.132887695 -1.270020715 
    ##           76           77           78           79           80 
    ##  0.772859702  2.842390033  0.919497035 -0.204067964  1.574898720 
    ##           81           82           83           84           85 
    ## -1.328422764 -0.942037909  0.558905683  0.927676137  0.637693575 
    ##           86           87           88           89           90 
    ##  0.668792582  1.820020782  0.342605380  1.503010395  3.129932642 
    ##           91           92           93           94           95 
    ##  0.232452135 -1.762860682  2.973852927  1.063499601  0.440495772 
    ##           96           97           98           99          100 
    ##  0.460890906 -1.639315743 -3.073406645 -2.974481076  1.731744015

``` r
fit$residuals
```

    ##          1          2          3          4          5          6 
    ##  0.9068832  6.6911177  0.6171562 -3.4957205 -0.1730067  2.5308735 
    ##          7          8          9         10         11         12 
    ## -0.6609824  1.1694370 -3.6879147 -0.9994137  0.4250052 -2.1877276 
    ##         13         14         15         16         17         18 
    ##  2.4787566 -2.6239326  6.0186112  2.1491064  3.5396122 -1.0151518 
    ##         19         20         21         22         23         24 
    ##  0.6107878  1.5608568 -1.7883827 -2.8992153  1.9864484  1.6646519 
    ##         25         26         27         28         29         30 
    ##  4.6019273  1.9450075 -3.3652972 -2.0789109 -0.3025080  2.1038441 
    ##         31         32         33         34         35         36 
    ## -4.4919726  3.7837454 -0.2754491  1.9132797 -1.4234406 -1.5178726 
    ##         37         38         39         40         41         42 
    ## -1.0894095 -1.0510206 -0.8385669  2.9300713  0.2071389  2.0864957 
    ##         43         44         45         46         47         48 
    ##  0.3268249 -0.4304703  2.4027358 -1.2395949 -0.9113377  1.2386378 
    ##         49         50         51         52         53         54 
    ## -2.3856852 -6.0825407  1.0031774  3.6341138  4.3301870  1.1018317 
    ##         55         56         57         58         59         60 
    ## -3.0463013 -4.0775992 -4.2332713  1.2074113 -2.4050525  3.2267260 
    ##         61         62         63         64         65         66 
    ## -4.1793871 -0.9605932  0.2880749  1.3143361 -0.9016654 -1.6063500 
    ##         67         68         69         70         71         72 
    ## -5.8683010 -6.7393014 -1.2611233 -1.6332521  3.8204370  0.1726141 
    ##         73         74         75         76         77         78 
    ## -0.8353533  4.3736969 -1.3348155  0.1127017  1.0425812  1.9175052 
    ##         79         80         81         82         83         84 
    ## -2.8760678 -2.8420326  1.3121663  4.8776927 -3.4135833 -0.2074437 
    ##         85         86         87         88         89         90 
    ##  1.2189353  1.5127400 -2.8876042 -0.2006822 -0.3640998 -0.7723966 
    ##         91         92         93         94         95         96 
    ##  1.8696075 -0.4030820 -0.6627004 -3.3817694  1.5040190  1.2382656 
    ##         97         98         99        100 
    ## -1.1195507  0.9039698  6.4384634  0.9186376

``` r
anova(fit)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Y
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## X1         1  20.85  20.846  2.8666   0.09364 .  
    ## X2         1 292.80 292.800 40.2635 7.092e-09 ***
    ## Residuals 97 705.39   7.272                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(fit$residuals~fit$fitted.values)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-10-1.png)

#### Example: Quantitative trait loci analysis for cholesterol levels

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

![](estimation_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
lm.598=lm(Y~X[,598])
lm.847=lm(Y~X[,847])

par(mfrow=c(1,2))
plot(Y~X[,598], main="Marker 598")
abline(lm.598)
plot(Y~X[,847], main="Marker 847")
abline(lm.847)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-12-2.png)

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

### Example: Simulation data

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

![](estimation_files/figure-markdown_github/unnamed-chunk-15-1.png)

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

![](estimation_files/figure-markdown_github/unnamed-chunk-17-1.png)

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

![](estimation_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
ridge_lm=lm.ridge(Y~X,lambda=seq(0,20,1))
plot(seq(0,20,1),ridge_lm$GCV)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-19-1.png)

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
