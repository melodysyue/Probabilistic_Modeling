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
-   [Logistic regression](#logistic-regression)
    -   [Example: Probability of admission](#example-probability-of-admission)
-   [Model selection and cross validation](#model-selection-and-cross-validation)
    -   [Ridge regression or L2 regularized linear regression](#ridge-regression-or-l2-regularized-linear-regression)
    -   [LASSO regression or L1 regularized linear regression](#lasso-regression-or-l1-regularized-linear-regression)
    -   [LOOCV (leave one out cross validation)](#loocv-leave-one-out-cross-validation)
    -   [Example: Simulation data](#example-simulation-data-1)
    -   [Example: Insulin](#example-insulin)
    -   [Example: Quantitative trait loci analysis for cholesterol levels](#example-quantitative-trait-loci-analysis-for-cholesterol-levels-1)
    -   [Example: Prostate cancer](#example-prostate-cancer)

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

    ## Residual Standard Error=3.1703
    ## R-Square=0.0695
    ## F-statistic (df=2, 97)=3.6224
    ## p-value=0.0304
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   1.4152  0.7757  1.8244   0.0712
    ## X1         -0.0134  0.0675 -0.1985   0.8431
    ## X2          0.2082  0.0778  2.6740   0.0088

``` r
## quadratic regression
X3=X1^2
Y2=1.1+0.02*X1+0.3*X3+rnorm(100,mean=0,sd=3)
X=cbind(X1,X3)
rr=lsfit(X,Y2)
ls.print(rr)
```

    ## Residual Standard Error=2.9134
    ## R-Square=0.9904
    ## F-statistic (df=2, 97)=5024.191
    ## p-value=0
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   2.1079  1.1191  1.8835   0.0626
    ## X1         -0.0816  0.2339 -0.3491   0.7278
    ## X3          0.3015  0.0112 26.9033   0.0000

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
    ## -8.2611 -1.7342 -0.2275  1.6679  7.8275 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  1.41521    0.77573   1.824   0.0712 . 
    ## X1          -0.01340    0.06751  -0.198   0.8431   
    ## X2           0.20815    0.07784   2.674   0.0088 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.17 on 97 degrees of freedom
    ## Multiple R-squared:  0.0695, Adjusted R-squared:  0.05031 
    ## F-statistic: 3.622 on 2 and 97 DF,  p-value: 0.03039

``` r
rr=lsfit(X,Y)
ls.print(rr)
```

    ## Residual Standard Error=3.2841
    ## R-Square=0.0015
    ## F-statistic (df=2, 97)=0.0731
    ## p-value=0.9296
    ## 
    ##           Estimate Std.Err t-value Pr(>|t|)
    ## Intercept   0.5787  1.2615  0.4587   0.6475
    ## X1          0.0404  0.2636  0.1531   0.8786
    ## X3         -0.0030  0.0126 -0.2405   0.8105

``` r
fit$coefficients
```

    ## (Intercept)          X1          X2 
    ##  1.41521420 -0.01339771  0.20815063

``` r
fit$fitted.values
```

    ##           1           2           3           4           5           6 
    ##  1.30553564  1.16023687 -0.08947948  1.51654120 -0.07873965  1.91418038 
    ##           7           8           9          10          11          12 
    ##  1.72268051  1.01886303  1.07450762  0.44089727  0.02581148 -0.29050799 
    ##          13          14          15          16          17          18 
    ##  0.81088580  2.12683769  0.08418170  1.50010527 -0.17472938  0.23322615 
    ##          19          20          21          22          23          24 
    ## -0.91813183 -0.54973632 -0.52799892  0.05631037  0.39388312  1.76486788 
    ##          25          26          27          28          29          30 
    ## -0.32222663  0.40320167 -0.19819361  1.11675084  1.20266625  0.49730140 
    ##          31          32          33          34          35          36 
    ##  0.49820652  0.60950864  1.10351390  0.18199793 -0.13696928  0.63841066 
    ##          37          38          39          40          41          42 
    ##  0.96832391 -0.57125394  0.30111258  1.12275276  2.19577156 -0.61441515 
    ##          43          44          45          46          47          48 
    ##  0.76392432  0.36737786 -0.56139503  0.34965453 -0.42658522 -0.74948767 
    ##          49          50          51          52          53          54 
    ##  0.10232253  0.22893775 -0.26150315  1.73000856  1.54560698  1.89866501 
    ##          55          56          57          58          59          60 
    ##  0.90762356  1.59848388  1.46049350  2.24612033  0.01259383  0.56512344 
    ##          61          62          63          64          65          66 
    ##  0.33336038  0.92100275  0.22247694  0.89046708 -0.65608552  0.08644721 
    ##          67          68          69          70          71          72 
    ##  0.30520816  0.32789886 -1.31212173  0.51842354  1.84186309  1.19390762 
    ##          73          74          75          76          77          78 
    ##  0.56075206 -0.22535107  1.55784283 -0.20646345  0.03455207 -0.44807093 
    ##          79          80          81          82          83          84 
    ##  1.44858170  1.72538498 -0.88767864  0.40393734  1.85558924 -0.46161631 
    ##          85          86          87          88          89          90 
    ##  1.16859047 -0.06911211  0.56736532  0.56334905  1.37565883  0.42599257 
    ##          91          92          93          94          95          96 
    ##  2.12620599  1.26992300  1.80832779  2.92690852  1.47049025  0.94363533 
    ##          97          98          99         100 
    ## -0.40177256  0.71872416  0.27459577  0.72167487

``` r
fit$residuals
```

    ##           1           2           3           4           5           6 
    ## -3.94060126  3.40928089 -2.67782214  4.73936481 -1.97280259  1.80222820 
    ##           7           8           9          10          11          12 
    ## -4.50389364 -3.02051373 -0.01339654  6.33773180  1.04094596  2.10946237 
    ##          13          14          15          16          17          18 
    ##  0.97162750 -0.91948540 -0.69449391  0.83283158 -1.50022949 -0.75183156 
    ##          19          20          21          22          23          24 
    ##  0.82010405  2.32491076 -3.62332563 -2.22216198 -3.85985143 -1.36759427 
    ##          25          26          27          28          29          30 
    ##  0.63795276 -2.38348977  3.46941620  3.58270999 -4.43234659  0.02754703 
    ##          31          32          33          34          35          36 
    ## -1.44628663  7.82754332 -0.69299500 -2.42730814 -0.35999297  3.46840897 
    ##          37          38          39          40          41          42 
    ##  3.97985804  0.18490653 -1.92708675 -6.50706550  0.82179923 -0.38528476 
    ##          43          44          45          46          47          48 
    ## -0.78851512  3.65870478  4.40304405  1.62309783  1.24296297 -1.69311713 
    ##          49          50          51          52          53          54 
    ##  1.00352143 -8.26113022 -3.92363058 -2.95659454  0.71315110 -0.17758038 
    ##          55          56          57          58          59          60 
    ## -5.58683935  1.33683518  5.10884535 -1.85749664 -2.01942314 -0.13972097 
    ##          61          62          63          64          65          66 
    ##  3.46050208  6.01757412  4.63859914 -5.53939240 -1.21145577  3.07520939 
    ##          67          68          69          70          71          72 
    ## -6.88793968  5.32236578 -3.07290097  0.81632638 -1.41101906 -5.32652050 
    ##          73          74          75          76          77          78 
    ## -1.11999500  3.83896187  1.56899778 -1.07117304  0.56440938  3.75667270 
    ##          79          80          81          82          83          84 
    ## -1.60770821  5.03726943  2.69688755  0.68722960 -0.27750598 -1.33293513 
    ##          85          86          87          88          89          90 
    ##  1.53205630  1.06364537  1.37198397 -4.54767246  0.12611780 -1.11787075 
    ##          91          92          93          94          95          96 
    ##  5.10102492 -1.12807089 -0.53920822 -1.07550355  0.92142545 -1.39774062 
    ##          97          98          99         100 
    ## -1.18049002  4.19012604 -0.41081134 -3.97635634

``` r
anova(fit)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Y
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## X1         1   0.95   0.952  0.0948 0.758880   
    ## X2         1  71.86  71.864  7.1501 0.008797 **
    ## Residuals 97 974.93  10.051                    
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

Logistic regression
-------------------

Logistic regression, also called a logit model, is used to model a dichotomous categorical outcome variables (binary). In the logit model, the log odds of the outcome is modeled as a linear combination of the predictor variables.

#### Example: Probability of admission

Example: A researcher is interested in how variables, such as GRE (Graduate Record Exam scores), GPA (grade point average) and prestige of the undergraduate institution, effect admission into graduate school. The response variable, admit/don’t admit, is a binary variable. A rank of 1 has the highest prestige, and rank of 4 has the lowest.

``` r
mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
## view the first few rows of the data
head(mydata)
```

    ##   admit gre  gpa rank
    ## 1     0 380 3.61    3
    ## 2     1 660 3.67    3
    ## 3     1 800 4.00    1
    ## 4     1 640 3.19    4
    ## 5     0 520 2.93    4
    ## 6     1 760 3.00    2

``` r
summary(mydata)
```

    ##      admit             gre             gpa             rank      
    ##  Min.   :0.0000   Min.   :220.0   Min.   :2.260   Min.   :1.000  
    ##  1st Qu.:0.0000   1st Qu.:520.0   1st Qu.:3.130   1st Qu.:2.000  
    ##  Median :0.0000   Median :580.0   Median :3.395   Median :2.000  
    ##  Mean   :0.3175   Mean   :587.7   Mean   :3.390   Mean   :2.485  
    ##  3rd Qu.:1.0000   3rd Qu.:660.0   3rd Qu.:3.670   3rd Qu.:3.000  
    ##  Max.   :1.0000   Max.   :800.0   Max.   :4.000   Max.   :4.000

``` r
xtabs(~admit+rank, data=mydata)
```

    ##      rank
    ## admit  1  2  3  4
    ##     0 28 97 93 55
    ##     1 33 54 28 12

``` r
mylogit=glm(admit~gre+gpa+rank, data=mydata,family="binomial")
summary(mylogit)
```

    ## 
    ## Call:
    ## glm(formula = admit ~ gre + gpa + rank, family = "binomial", 
    ##     data = mydata)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5802  -0.8848  -0.6382   1.1575   2.1732  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -3.449548   1.132846  -3.045  0.00233 ** 
    ## gre          0.002294   0.001092   2.101  0.03564 *  
    ## gpa          0.777014   0.327484   2.373  0.01766 *  
    ## rank        -0.560031   0.127137  -4.405 1.06e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 499.98  on 399  degrees of freedom
    ## Residual deviance: 459.44  on 396  degrees of freedom
    ## AIC: 467.44
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
confint(mylogit) ##get CIs using profiled log-likelihood
```

    ## Waiting for profiling to be done...

    ##                     2.5 %       97.5 %
    ## (Intercept) -5.7109591680 -1.260314066
    ## gre          0.0001715446  0.004461385
    ## gpa          0.1415710585  1.428341503
    ## rank        -0.8149612229 -0.315479733

``` r
confint.default(mylogit) ##get CIs using standard errors
```

    ##                     2.5 %       97.5 %
    ## (Intercept) -5.6698857745 -1.229211021
    ## gre          0.0001539942  0.004433925
    ## gpa          0.1351569663  1.418870181
    ## rank        -0.8092153067 -0.310847467

``` r
mydata$rank=factor(mydata$rank)
mylogit=glm(admit~gre+gpa+rank, data=mydata,family="binomial")
summary(mylogit)
```

    ## 
    ## Call:
    ## glm(formula = admit ~ gre + gpa + rank, family = "binomial", 
    ##     data = mydata)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.6268  -0.8662  -0.6388   1.1490   2.0790  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -3.989979   1.139951  -3.500 0.000465 ***
    ## gre          0.002264   0.001094   2.070 0.038465 *  
    ## gpa          0.804038   0.331819   2.423 0.015388 *  
    ## rank2       -0.675443   0.316490  -2.134 0.032829 *  
    ## rank3       -1.340204   0.345306  -3.881 0.000104 ***
    ## rank4       -1.551464   0.417832  -3.713 0.000205 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 499.98  on 399  degrees of freedom
    ## Residual deviance: 458.52  on 394  degrees of freedom
    ## AIC: 470.52
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
newdata1 <- with(mydata, data.frame(gre = mean(gre), gpa = mean(gpa), rank = factor(1:4)))
newdata1$rankP <- predict(mylogit, newdata = newdata1, type = "response")
newdata1
```

    ##     gre    gpa rank     rankP
    ## 1 587.7 3.3899    1 0.5166016
    ## 2 587.7 3.3899    2 0.3522846
    ## 3 587.7 3.3899    3 0.2186120
    ## 4 587.7 3.3899    4 0.1846684

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

![](estimation_files/figure-markdown_github/unnamed-chunk-16-1.png)

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

![](estimation_files/figure-markdown_github/unnamed-chunk-18-1.png)

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

![](estimation_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
ridge_lm=lm.ridge(Y~X,lambda=seq(0,20,1))
plot(seq(0,20,1),ridge_lm$GCV)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-20-1.png)

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

### Example: Quantitative trait loci analysis for cholesterol levels

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

### Example: Prostate cancer

``` r
library(ElemStatLearn) ## contains the dta
data(prostate)
head(prostate)
```

    ##       lcavol  lweight age      lbph svi       lcp gleason pgg45       lpsa
    ## 1 -0.5798185 2.769459  50 -1.386294   0 -1.386294       6     0 -0.4307829
    ## 2 -0.9942523 3.319626  58 -1.386294   0 -1.386294       6     0 -0.1625189
    ## 3 -0.5108256 2.691243  74 -1.386294   0 -1.386294       7    20 -0.1625189
    ## 4 -1.2039728 3.282789  58 -1.386294   0 -1.386294       6     0 -0.1625189
    ## 5  0.7514161 3.432373  62 -1.386294   0 -1.386294       6     0  0.3715636
    ## 6 -1.0498221 3.228826  50 -1.386294   0 -1.386294       6     0  0.7654678
    ##   train
    ## 1  TRUE
    ## 2  TRUE
    ## 3  TRUE
    ## 4  TRUE
    ## 5  TRUE
    ## 6  TRUE

``` r
str(prostate)
```

    ## 'data.frame':    97 obs. of  10 variables:
    ##  $ lcavol : num  -0.58 -0.994 -0.511 -1.204 0.751 ...
    ##  $ lweight: num  2.77 3.32 2.69 3.28 3.43 ...
    ##  $ age    : int  50 58 74 58 62 50 64 58 47 63 ...
    ##  $ lbph   : num  -1.39 -1.39 -1.39 -1.39 -1.39 ...
    ##  $ svi    : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ lcp    : num  -1.39 -1.39 -1.39 -1.39 -1.39 ...
    ##  $ gleason: int  6 6 7 6 6 6 6 6 6 6 ...
    ##  $ pgg45  : int  0 0 20 0 0 0 0 0 0 0 ...
    ##  $ lpsa   : num  -0.431 -0.163 -0.163 -0.163 0.372 ...
    ##  $ train  : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...

``` r
pairs(prostate)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
table(prostate$gleason) ## very few data points for 8 and 9
```

    ## 
    ##  6  7  8  9 
    ## 35 56  1  5

``` r
prostate$gleason=ifelse(prostate$gleason==6, 0 ,1) ## group 7,8,9 together
table(prostate$gleason)
```

    ## 
    ##  0  1 
    ## 35 62

``` r
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
prostate.cor=cor(prostate)
corrplot(prostate.cor)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-22-2.png)

``` r
str(prostate)
```

    ## 'data.frame':    97 obs. of  10 variables:
    ##  $ lcavol : num  -0.58 -0.994 -0.511 -1.204 0.751 ...
    ##  $ lweight: num  2.77 3.32 2.69 3.28 3.43 ...
    ##  $ age    : int  50 58 74 58 62 50 64 58 47 63 ...
    ##  $ lbph   : num  -1.39 -1.39 -1.39 -1.39 -1.39 ...
    ##  $ svi    : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ lcp    : num  -1.39 -1.39 -1.39 -1.39 -1.39 ...
    ##  $ gleason: num  0 0 1 0 0 0 0 0 0 0 ...
    ##  $ pgg45  : int  0 0 20 0 0 0 0 0 0 0 ...
    ##  $ lpsa   : num  -0.431 -0.163 -0.163 -0.163 0.372 ...
    ##  $ train  : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...

``` r
train=subset(prostate,train==TRUE)[,1:9]
str(train)
```

    ## 'data.frame':    67 obs. of  9 variables:
    ##  $ lcavol : num  -0.58 -0.994 -0.511 -1.204 0.751 ...
    ##  $ lweight: num  2.77 3.32 2.69 3.28 3.43 ...
    ##  $ age    : int  50 58 74 58 62 50 58 65 63 63 ...
    ##  $ lbph   : num  -1.39 -1.39 -1.39 -1.39 -1.39 ...
    ##  $ svi    : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ lcp    : num  -1.39 -1.39 -1.39 -1.39 -1.39 ...
    ##  $ gleason: num  0 0 1 0 0 0 0 0 0 1 ...
    ##  $ pgg45  : int  0 0 20 0 0 0 0 0 0 30 ...
    ##  $ lpsa   : num  -0.431 -0.163 -0.163 -0.163 0.372 ...

``` r
test=subset(prostate,train==FALSE)[,1:9]
str(test)
```

    ## 'data.frame':    30 obs. of  9 variables:
    ##  $ lcavol : num  0.737 -0.777 0.223 1.206 2.059 ...
    ##  $ lweight: num  3.47 3.54 3.24 3.44 3.5 ...
    ##  $ age    : int  64 47 63 57 60 69 68 67 65 54 ...
    ##  $ lbph   : num  0.615 -1.386 -1.386 -1.386 1.475 ...
    ##  $ svi    : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ lcp    : num  -1.386 -1.386 -1.386 -0.431 1.348 ...
    ##  $ gleason: num  0 0 0 1 1 0 0 1 0 0 ...
    ##  $ pgg45  : int  0 0 0 5 20 0 0 20 0 0 ...
    ##  $ lpsa   : num  0.765 1.047 1.047 1.399 1.658 ...

``` r
table(prostate$train)
```

    ## 
    ## FALSE  TRUE 
    ##    30    67

Model with OLS

``` r
library(leaps)
subfit=regsubsets(lpsa~., data=train)
sub.sum=summary(subfit)
which.min(sub.sum$bic)
```

    ## [1] 3

``` r
plot(subfit, scale="bic") 
```

![](estimation_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
## It shows lcavol, lweight and gleason are selected in the best 3-variable model, which gives the minimum BIC. 

ols=lm(lpsa~lcavol+lweight+gleason,data=train)
plot(ols$fitted.values, train$lpsa, xlab="prediction from train dataset",ylab="observation")
```

![](estimation_files/figure-markdown_github/unnamed-chunk-23-2.png)

``` r
pred=predict(ols,newdata=test)
plot(pred, test$lpsa, xlab="Prediction from test dataset", ylab="observation")
```

![](estimation_files/figure-markdown_github/unnamed-chunk-23-3.png)

``` r
MSE.ols=mean((test$lpsa-pred)^2)
MSE.ols
```

    ## [1] 0.5084126

Model with Ridge regression

``` r
x.test=as.matrix(test[,1:8])
y.test=test[,9]
x=as.matrix(train[,1:8])
y=train[,9]
library(glmnet)
```

    ## Warning: package 'glmnet' was built under R version 3.4.4

    ## Loading required package: Matrix

    ## Loading required package: foreach

    ## Loaded glmnet 2.0-16

``` r
ridge=cv.glmnet(x,y,alpha=0)
ridge.bs.lambda=ridge$lambda.min
ridge.best=glmnet(x,y,alpha=0, lambda=ridge.bs.lambda)
pred.ridge=predict(ridge.best, newx=x.test)
MSE.ridge=mean((y.test-pred.ridge)^2)
MSE.ridge
```

    ## [1] 0.4785047

Model with Lasso regression

``` r
lasso=cv.glmnet(x,y,alpha=1)
lasso.bs.lambda=lasso$lambda.min
lasso.best=glmnet(x,y,alpha=1, lambda=ridge.bs.lambda)
pred.lasso=predict(lasso.best, newx=x.test)
MSE.lasso=mean((y.test-pred.lasso)^2)
MSE.lasso
```

    ## [1] 0.4446758

Model with regression tree

``` r
set.seed(1)
library(rpart)
tree=rpart(lpsa~., data=train)
library(partykit)
```

    ## Warning: package 'partykit' was built under R version 3.4.4

    ## Loading required package: grid

    ## Loading required package: libcoin

    ## Loading required package: mvtnorm

``` r
plot(as.party(tree))
```

![](estimation_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
pred.tree=predict(tree, newx=x.test)
MSE.tree=mean((y.test-pred.tree)^2)
```

    ## Warning in y.test - pred.tree: longer object length is not a multiple of
    ## shorter object length

``` r
MSE.tree
```

    ## [1] 1.854035

Note: regression tree is very unstable and tends to overfit, since it loses the information from the continuous variable. Such instability can be solved with bagging method, which is basically a collection of decision tress, and take the average of the output from these trees.

Model with random forest regression

``` r
library(randomForest)
```

    ## randomForest 4.6-12

    ## Type rfNews() to see new features/changes/bug fixes.

``` r
set.seed(1)
rf=randomForest(lpsa~.,data=train)
plot(rf)
```

![](estimation_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
which.min(rf$mse)
```

    ## [1] 82

``` r
rf.best=randomForest(lpsa~.,data=train, ntree=which.min(rf$mse))
importance(rf.best)
```

    ##         IncNodePurity
    ## lcavol      28.239422
    ## lweight     16.504314
    ## age          6.708906
    ## lbph         6.312513
    ## svi          7.873221
    ## lcp          8.666698
    ## gleason      4.965704
    ## pgg45        7.991093

``` r
pred.rf=predict(rf.best,newdata=x.test)
MSE.rf=mean((y.test-pred.rf)^2)
MSE.rf
```

    ## [1] 0.5219353

Note: A primary disadvantage of random forests is that the results are not easily interpretable: that is, if you would like to draw conclusions about the meaning of the classification model, random forests may not be the best choice.

Difference between decision tree modeling and linear regression: Decision Tree and Linear Regression are both supervised learning algorithms. While decision tree is easy to interpret, linear regression is only good when relationships between variables are linear and also when you need to also find the marginal effect. When the data is not linear, you should choose decision tress and random forests.
