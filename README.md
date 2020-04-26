# tempoSig
Mutation Signature Extraction using Maximum Likelihood

## Overview
**tempoSig** implements maximum likelihood-based extraction of mutational signature proportions of a set of mutation count data under a known set of input signature lists. It uses the same algorithm as in [mutational-signatures](https://github.com/mskcc/mutation-signatures), re-implemented in R/C++ to achieve a speed-up of ~400x. The speed-up enables fast estimation of p-values via permutation-based sampling. The basic object (S4 class) can store input data, reference signature list, output exposure of samples, and p-values. Utilities for plotting and file ouput are also included. 

## Installation
Compilation requires GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/). In Ubuntu Linux,

    $ sudo apt-get install libgsl-dev
    
In OS-X,

    $ brew install gsl

Install [Rcpp](https://cran.r-project.org/package=Rcpp) if not installed already, as well as [gtools](https://cran.r-project.org/package=gtools).

**tempoSig** can then be installed in R by

    > devtools::install_github("mskcc/tempoSig")
