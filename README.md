# tempoSig
Mutation Signature Extraction using Maximum Likelihood
## Overview
**tempoSig** implements maximum likelihood-based extraction of mutational signature proportions of a set of mutation count data under a known set of input signature lists. It uses the same algorithm as in [mutational-signatures](https://github.com/mskcc/mutation-signatures), re-implemented in R/C++ to achieve a speed-up of ~400x. The speed-up enables a fast estimation of permutation-based p-values. The basic object (S4 class) can store input data, reference signature list, output exposure of samples, and p-values. Utilities for plotting and file ouput are also included. 

