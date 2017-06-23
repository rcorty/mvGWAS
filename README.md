
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

Create an mvGWAS object
-----------------------

Once you've installed the package and its dependencies, the first step is to create an mvGWAS object. You'll need two things on your filesystem:

1.  a directory with some .VCF or .VCF.GZ files in it. For now, set it up so there are just a few (&lt;10) files each of which has only a few variants (&lt;1000) on a few participants (&lt;100).
2.  a file with some phenotype data, either a .RDS file or a .CSV file. The first column of the data.frame (that's stored as a .RDS or that results from `readr::read_csv(phenotype_file)`) must be "ID".

``` r
my_gwas <- mvGWAS$new(phenotype_file = 'my_phenotype_file.csv',
                      genotype_directory = 'genotype_dir')
```

Note that `mvGWAS` is a generator function for a reference class. That means that in the example above, `my_gwas` is where all the action is going to be going forward. (It is an object of class "mvGWAS".)

Conduct a scan using your local computer
----------------------------------------

Since there's not much data, we can do the whole analysis on the local machine. To conduct a scan we use the `conduct_scan` function. It takes two arguments: `mean_formula` and `var_formula`. They use the [model formula system](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/formula.html) that is typical in R.

`mean_formula` is a two-sided formula. Its right side is just the phenotype to be analyzed. Its left side is covariates and keywords, separated by `+`. There are three keywords available:

1.  GT -- the most likely genotype, a factor with up to 3 levels, which becomes a two-column design matrix
2.  DS -- the expected alternate allele dosage, a numeric in \[0, 2\]
3.  GP -- the probability of each genotype, three numeric quantities in \[0, 1\], which become a two-column design matrix

`var_formula` is a one-sided formula. It has no left side. Its right side follows the same rules as `mean_formula`.

``` r
my_gwas$conduct_scan(mean_formula = my_trait ~ PC1 + PC2 + sex + age + DS,
                     var_formula = ~ PC1 + PC2 + sex + GT)
```

By default, `conduct_scan()` runs on the local machine, using as many cores as are available (based on `parallel::detectCores()`). If you have a multi-core computer but don't want to use them all in `conduct_scan()`, you can pass in the maximum number of cores to use as `num_cores`. For example:

``` r
my_gwas$conduct_scan(mean_formula = my_trait ~ PC1 + PC2 + sex + age + DS,
                     var_formula = ~ PC1 + PC2 + sex + GT,
                     num_cores = 2)
```

If `conduct_scan()` returns `TRUE`, it worked without error and the results are accessible as `my_gwas$results`.

View the results
----------------

``` r
my_gwas$qq_plot()
my_gwas$manhattan_plot()
```
