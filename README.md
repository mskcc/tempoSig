# tempoSig
Mutation Signature Extraction using Maximum Likelihood

## Overview
**tempoSig** implements maximum likelihood-based extraction of mutational signature proportions of a set of mutation count data under a known set of input signature lists. It uses the same algorithm as in [mutational-signatures](https://github.com/mskcc/mutation-signatures), re-implemented in R/C++ to achieve a speed-up of ~400x. The speed-up enables fast estimation of p-values via permutation-based sampling. The basic object (S4 class) can store input data, reference signature list, output exposure of samples, and p-values. Utilities for plotting and file ouput are also included. 

## Algorithm
Input data is of the form of catalog matrix:

Mutation context | Tumor Sample Barcode 1 | Tumor Sample Barcode 2
---------------- | ---------------------- | ----------------------
A[C>A]A          |                      0 |                      1
A[C>A]C          |                      3 |                      0
A[C>A]G          |                      2 |                      5
A[C>A]T          |                      0 |                      2
C[C>A]A          |                      0 |                      0
C[C>A]C          |                      1 |                      1
C[C>A]G          |                      0 |                      0
C[C>A]T          |                      1 |                      1

Mutation context is the set of categories to which mutation data from sequencing experiments have been classified. Trinucleotide contexts shown above are the pyrimidine bases before and after mutation flanked by upstream and downstream nucleotides (96 in total). Each column corresponding to **Tumor_Sample_Barcode** contains non-negative counts of single nucleotide variants (SNVs). This trinucleotide matrix can be generated from [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) files using the R package [maftools](https://bioconductor.org/packages/maftools/). MAF files can be generated from VCF files using [vcf2maf](https://github.com/mskcc/vcf2maf).

The other input is the set of reference signature proportions:

Mutation context | SBS2 | SBS3 | SBS4
---------------- | ---- | ---- | ----
A[C>A]A          | 6e-7 | 0.02 | 0.04
A[C>A]C          | 2e-3 | 0.02 | 0.03

Both [version 2](https://github.com/mskcc/tempoSig/edit/master/inst/extdata/COSMIC_SNV_signatures_v2.txt) and [version 3](https://github.com/mskcc/tempoSig/edit/master/inst/extdata/COSMIC_SBS_signatures-v3.txt) tables of [COSMIC signature lists](https://cancer.sanger.ac.uk/cosmic/signatures) are included.

The "refitting" (as opposed to de novo discovery) of signature propotion solves the non-negative matrix factorization problem

    X = W * H
    
where **X** is the catalog matrix, **W** is the signature matrix (assumed to be known and fixed), and **H** is the exposure matrix of dimension (no. of reference signatures x no. of samples). The maximum likelihood estimation (MLE) algorithm formulates this problem in terms of a multinomial statistical model with observed counts **X** of categorical groups (trinucleotide contexts) mixed by given fixed proportions **W**. **tempoSig** uses the quasi-Newton [Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm](https://www.gnu.org/software/gsl/doc/html/multimin.html) for multi-dimensional optimization.
The output matrix is the transpose of **H**:

Tumor Sample Barcode | SBS2 | SBS3 | SBS4
-------------------- | ---- | ---- | ----
SAMPLE_1             | 0.00 | 0.01 | 0.25
SAMPLE_2             | 0.01 | 0.12 | 0.20

Each row is a vector of proportions that add up to 1.

To estimate p-values of significance of each proportion, the input data for each sample (columns of **X**) are randomly shuffled by permutation to sample the null distribution. The exposure vectors inferred from these null samples are compared with the observed vector from the original data, with the p-value defined as the proportion of null samples whose proprotions are higher than the observed values. 

## Installation
Compilation requires GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/). In Ubuntu Linux,

    $ sudo apt-get install libgsl-dev
    
In OS-X,

    $ brew install gsl

Install [Rcpp](https://cran.r-project.org/package=Rcpp) if not installed already, as well as [gtools](https://cran.r-project.org/package=gtools).

**tempoSig** can then be installed in R by

    > devtools::install_github("mskcc/tempoSig")
