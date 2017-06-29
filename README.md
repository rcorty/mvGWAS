
<!-- README.md is generated from README.Rmd. Please edit that file -->
mvGWAS
======

The goal of the mvGWAS package is to provide an easy-to-use interface for geneticists to conduct a mean-variance genome wide association study on a large number of genetic loci. This type of GWAS test each locus for association with phenotype mean and variance.

At each locus, there are three questions (equivalently, three null hypotheses to test). Does an individual's genotype at this locus influence...

1.  the expected value of the phenotype?
2.  the residual variance of the phenotype?
3.  at least one of the expected value or the residual variance of the phenotype?

This package uses the double generalized linear model (DGLM) to address each of these questions. Github doesn't support math mode, so, if you want to see the model, check out [a related manuscript on BioArXiV](http://biorxiv.org/content/biorxiv/early/2017/06/20/149377.full.pdf) for now.

Installation
------------

You can install mvGWAS from github with:

``` r
install.packages("devtools")
devtools::install_github("rcorty/mvGWAS")
```

and load it with:

``` r
library(mvGWAS)
```
