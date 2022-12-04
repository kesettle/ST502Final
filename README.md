ST 502 Final
================
Katelyn Settlemyre, Julia Farrell
2022-12-04

## Introduction

*What is the general idea of testing for homogeneity? Perhaps also
explain the purpose and outline of this project.*

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

    ##   Surgical Site Infections Pneumonia Infections Bloodstream Infections
    ## A                       41                   27                     51
    ## B                       36                    3                     40
    ## C                      169                  106                    109

## Deriving the Likelihood Ratio Test

*Derive the LRT test for homogeneity; some work already done on page 151
of our notes. Pages 129, 132, 137, 139 may also be useful?*

## Simulation

*Simulate test with data. Potentially make use of code already given in
notes. Two multinomial case only, with 3 categories in each multinomial.
Plot summaries of simulations. Set seed for reproduction purposes.*

Goal:

Process:

``` r
#set seed for reproduction
set.seed(17)
```
