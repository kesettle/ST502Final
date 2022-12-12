ST 502 Final
================
Katelyn Settlemyre, Julia Farrell
2022-12-04

## Introduction

This report will explore the Chi-Square test for homogeneity in detail.
We will derive the likelihood ratio test (LRT) statistic used to conduct
this test and explain the Pearson Chi-Square statistic that can be used
as an approximation. Once the theory is well established, we then
conduct a simulation. The goal is to determine how well the asymptotic
rejection region performs at controlling the alpha level of the Pearson
Chi-Square test and to determine the power of the asymptotic test when
comparing certain alternative hypotheses.

The Chi-Square test for homogeneity is used in a specific case:
comparing $J$ multinomials with $I$ classes (with $I,J \in \mathbb{N}$),
where the researcher is interested in determining if the probabilities
of each cell are the same across every multinomial.

## Data Example

*Use hospital data given to conduct a* $\chi^2$ *test for homogeneity.*

``` r
#create and print matrix of hospital data
rows <- rbind(a=c(41, 27, 51), b=c(36,3,40), c=c(169, 106, 109))
(hospDat <- matrix(data=rows, nrow = 3, ncol = 3, 
                   dimnames = list(c("A", "B", "C"), 
                                   c("Surgical Site Infections", 
                                     "Pneumonia Infections", 
                                     "Bloodstream Infections"))))
```

      Surgical Site Infections Pneumonia Infections Bloodstream Infections
    A                       41                   27                     51
    B                       36                    3                     40
    C                      169                  106                    109

``` r
summary(hospDat)
```

     Surgical Site Infections Pneumonia Infections Bloodstream Infections
     Min.   : 36.0            Min.   :  3.00       Min.   : 40.00        
     1st Qu.: 38.5            1st Qu.: 15.00       1st Qu.: 45.50        
     Median : 41.0            Median : 27.00       Median : 51.00        
     Mean   : 82.0            Mean   : 45.33       Mean   : 66.67        
     3rd Qu.:105.0            3rd Qu.: 66.50       3rd Qu.: 80.00        
     Max.   :169.0            Max.   :106.00       Max.   :109.00        

Now we will conduct a Chi-Square test for homogeneity using this sample
data.  
$H_0$: The distribution of infections is the same for each hospital  
$H_1$: The distribution of infections is the same for each hospital.

``` r
x = chisq.test(hospDat)
x
```


        Pearson's Chi-squared test

    data:  hospDat
    X-squared = 30.696, df = 4, p-value = 3.531e-06

``` r
x$expected
```

      Surgical Site Infections Pneumonia Infections Bloodstream Infections
    A                 50.29897             27.80756               40.89347
    B                 33.39175             18.46048               27.14777
    C                162.30928             89.73196              131.95876

The p-value of our chi-square test statistic is extremely small,
indicating that we reject the null hypothesis. The data in this example
provides support for the alternative hypothesis that the multinomials
from these hospitals are not homogeneous.

## Deriving the Likelihood Ratio Test

The goal is to derive the likelihood ratio test for a generalized case
comparing $J$ independent multinomial distributions, each with $I$
categories.

Let there be $J$ independent multinomial distributions, each with $I$
categories, where $I,J\in \mathbb{N}$. We want to test the hypothesis
$H_0 = \pi_{11}=\pi_{12}=...=\pi_{1J}, \pi_{21}=\pi_{22}=...=\pi_{1J}, \pi_{I1}=\pi_{I2}=...=\pi_{IJ}$
vs. $H_a:$ at least one probability differs

To derive the likelihood ratio test, we will initially look at the
likelihood function, which is just the product of the $J$ multinomials.

Expected Counts Under $H_0$:
$L(\pi_{ij}'s)=\prod_{j=1}^{J}{\binom{0}{0}\cdot \pi_{ij}^{n_{ij}}}\cdot\pi_{2j}^{n_{2j}}\cdot ... \cdot\pi_{IJ}^{n_{IJ}}$
$\propto \prod_{j=1}^{J}\prod_{i=1}^{I}{\pi_{ij}^{n_{ij}}}$ subject to
constraints \*\*\*

Under the null hypothesis,$\pi_{11}=\pi_{12}=...=\pi_{1J}$, so we will
replace this with the common value $\pi_1$. Similarly, we will continue
forward considering the common values $\pi_1,..., \pi_I$. There is one
restriction on these probabilities to make them valid, which is
$\pi_1+\pi_2+...+ \pi_I=1$

Degrees of freedom for reference distribution:  
$dim(\Omega) = J\cdot (I-1)$  
$dim(\omega) = (I-1)$  
$df = J(I-1) - (I-1) = (I-1)(J-1)$

.  
. More derivation steps  
.  
.  
The likelihood ratio test for the homogeneity hypothesis being tested is
therefore given by
$LRT=2\sum_{j=1}^J\sum_{i=1}^I{Obs_{ij} \cdot ln(\frac{Obs_{ij}}{Exp_{ij}})}$

## Simulation

Goal:  
- Determine how well the asymptotic rejection region performs at
controlling $\alpha$.  
- Determine the power of the asymptotic test when comparing certain
alternative situations.

Process:  
For the simulation process, we will do the case of 2 multinomials each
with 3 classes. We will generate many, many tables of these multinomial
pairs, each of varying sample sizes, and from 3 cases of probabilities:
the equal probability case, and two mixed cases. After generating our
tables, we will then conduct a theoretical $\chi^2$ test and obtain a
proportion of times we reject $H_0$ over all tests for each combination
of sample sizes under each probability case. This is our approximated
$\alpha$. Similarly, we will look at the power of the test by comparing
multinomials of varying sample sizes with different probability cases:
equal case vs the first mixed case, equal case vs the second mixed case,
first mixed case vs the second mixed case.  
For both the $\alpha$ control and the power inspection, we will graph
the comparisons.

``` r
#set seed for reproduction
set.seed(17)

#sample sizes
n1 <- c(20,30,50,100)
n2 <- n1

#probabilities - three cases being tested
p1 <- c(1/3, 1/3, 1/3) #equal
p2 <- c(0.1, 0.3, 0.6) #mixed 1
p3 <- c(0.1, 0.1, 0.8) #mixed 2


#Use N=50000 random tables (start lower, end up here)
N <- 50000

#number of classes & multinomials being compared
#classes I = 3   
#multinomials J = 2   

# function to create two multinomials,using sample size and probabilities input
### sampleSize: (int)
### probs: (list of doubles) length of p split determines number of classes. All = 3 for this simulation.
### returns the approximate alpha values from simulated multinomials

#create function to generate tables like in hospital example
multGenerate <- function(size1, size2, prob) {
  #generate multinomials
  mult1 <- rmultinom(1, size1, prob)
  mult2 <- rmultinom(1, size2, prob)
  #combine multinomials into table
  multTab <- cbind(mult1, mult2)
  #transpose table to have multinomials as the rows
  multTrans = t(multTab)
  #return table
  return(multTrans)
}

#testing in wrapper function
wrapAlpha <- function(size1, size2, prob){
  #generate table
  a <- multGenerate(size1, size2, prob)
  #compute chisq and compare to theoretical cutoff
  #here theoretical alpha = 0.05, but can easily update function to pick alpha
  #df = (I-1)(J-1) = (3-1)(2-1) = 2
  b <- chisq.test(a)
  c <- isTRUE(b$statistic > qchisq(0.95, 2))
  #return T/F value
  return(c)
}

##test wrapper; replicate many many times and store as vector
#alph <- replicate(5000, wrapAlpha(n1[2], n2[2],p1))
##take the mean of the vector to approximate alpha level
#mean(alph)

#set up vectors to store approximated alpha levels
p1Alphas <- c()
p2Alphas <- c()
p3Alphas <- c()

#for each combination of sample sizes, approximate alpha with each set of probabilities
for (i in 1:4) {
  for(j in 1:4){
    p1Alphas <- append(p1Alphas, mean(replicate(N,wrapAlpha(n1[i], n2[j], p1))))
    p2Alphas <- append(p2Alphas, mean(replicate(N,wrapAlpha(n1[i], n2[j], p2))))
    p3Alphas <- append(p3Alphas, mean(replicate(N,wrapAlpha(n1[i], n2[j], p3))))
  }
}


#combine all into data frame
n_1 <- c(20,20,20,20,30,30,30,30,50,50,50,50,100,100,100,100)
n_2 <- c(20,30,50,100,20,30,50,100,20,30,50,100,20,30,50,100)
alP1 <- replicate(16, "equal")
alP2 <- replicate(16, "mixed1")
alP3 <- replicate(16, "mixed2")
d1<-data.frame(alpha=p1Alphas, n1=n_1, n2=n_2, p=alP1)
d2<-data.frame(alpha=p2Alphas, n1=n_1, n2=n_2, p=alP2)
d3<-data.frame(alpha=p3Alphas, n1=n_1, n2=n_2, p=alP3)
alphaDat <- rbind(d1,d2,d3)

#plot
ggplot(alphaDat, aes(x=n2,y=alpha, group=p))+geom_line(aes(color=p))+facet_grid(.~n1)
```

![](data:image/png;base64,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)<!-- -->

``` r
#Power inspection
#compare 1-2, 1-3, 2-3

#Similar to multGenerate, but allows for 2 different probabilities
pwrCompare <- function(size1, size2, prob1, prob2) {
  #generate multinomials
  mult1 <- rmultinom(1, size1, prob1)
  mult2 <- rmultinom(1, size2, prob2)
  #combine multinomials into matrix
  multTab <- cbind(mult1, mult2)
  #transpose matrix
  multTrans = t(multTab)
  return(multTrans)
}

wrapPWR <- function(size1, size2, prob1, prob2){
  #generate table
  a <- pwrCompare(size1, size2, prob1, prob2)
  #compute chisq and compare to theoretical cutoff
  #here theoretical alpha = 0.05, but can easily update function to pick alpha
  #df = (I-1)(J-1) = (3-1)(2-1) = 2
  b <- chisq.test(a)
  c <- isTRUE(b$statistic > qchisq(0.95, 2))
  #return T/F value
  return(c)
}

##test wrapper; replicate many many times and store as vector
#pwr <- replicate(5000, wrapPWR(n1[3], n2[3],p1,p2))
##take the mean of the vector to approximate alpha level
#mean(pwr)

#set up vectors to store approximated power
pwrp1p2 <- c()
pwrp1p3 <- c()
pwrp2p3 <- c()

#for each combination of sample sizes, approximate power with each pair of probabilities
for (i in 1:4) {
  for(j in 1:4){
    pwrp1p2 <- append(pwrp1p2, mean(replicate(N,wrapPWR(n1[i], n2[j], p1,p2))))
    pwrp1p3 <- append(pwrp1p3, mean(replicate(N,wrapPWR(n1[i], n2[j], p1,p3))))
    pwrp2p3 <- append(pwrp2p3, mean(replicate(N,wrapPWR(n1[i], n2[j], p2,p3))))
  }
}


#combine all into data frame
be12 <- replicate(16, "1v2")
be13 <- replicate(16, "1v3")
be23 <- replicate(16, "2v3")
b1<-data.frame(power=pwrp1p2, n1=n_1, n2=n_2, p=be12)
b2<-data.frame(power=pwrp1p3, n1=n_1, n2=n_2, p=be13)
b3<-data.frame(power=pwrp2p3, n1=n_1, n2=n_2, p=be23)
powerDat <- rbind(b1,b2,b3)

#create graphs
ggplot(powerDat, aes(x=n2,y=power, group=p))+geom_line(aes(color=p))+facet_grid(.~n1)
```

![](data:image/png;base64,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)<!-- -->
