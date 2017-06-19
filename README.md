
<!-- README.md is generated from README.Rmd. Please edit that file -->
GWASviaDGLM
===========

The goal of GWASviaDGLM is to provide an easy-to-use interface for geneticists to conduct an mvGWAS -- a mean-variance genome wide association study. This study tests all the loci in the genome for association with the mean and variance of a phenotype.

At each locus, there are three questions (equivalently, three null hypotheses to test).

1.  Does an individual's genotype at this locus influence the expected value of his phenotype?
2.  Does an individual's genotype at this locus influence the residual variance of his phenotype?
3.  Does an individual's genotype at this locus influence influence at least one of the expected value or the residual variance of his phenotype?

This package uses the double generalized linear model (DGLM) to address each of these questions. This model can be written:

$$
y\_i \\sim \\text{N}(\\mu\_i, \\text{exp}(2\\nu\_i))\\\\
\\mu\_i = x\_i^T\\beta + q\_i^T\\alpha\\\\
\\nu\_i = z\_i^T\\gamma + q\_i^T\\theta
$$

where *y*<sub>*i*</sub> is the phenotype of individual *i*, *μ*<sub>*i*</sub> is its estimated expected value, 2exp*ν*<sub>*i*</sub> is its estimated residual variance,

Installation
------------

You can install GWASviaDGLM from github with:

``` r
# install.packages("devtools")
devtools::install_github("rcorty/GWASviaDGLM")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```
