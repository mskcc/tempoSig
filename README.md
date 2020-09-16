# tempoSig
Mutational Signature Extraction using Maximum Likelihood

## Overview
**tempoSig** implements maximum likelihood-based extraction of mutational signature proportions of a set of mutation count data under a known set of input signature lists. The basic algorithm is the same as in [mutation-signatures](https://github.com/mskcc/mutation-signatures), but re-implemention in R/C++ here enables a substantial speed-up of the order of ~100x. This speed-up allows for the fast estimation of p-values via permutation-based sampling. The basic object (S4 class) can store input data, reference signature list, output exposure of samples, and p-values. Utilities for plotting and file ouput are also included. 

## Algorithm
Input data are of the form of catalog matrix:

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

Mutation context | Signature.1 | Signature.2 | Signature.3 | Signature.4
---------------- | ---- | ---- | ---- | ----
A[C>A]A          | 9e-4 | 6e-7 | 0.02 | 0.04
A[C>A]C          | 2e-3 | 1e-4 | 0.02 | 0.03

Both [version 2](https://github.com/mskcc/tempoSig/blob/master/inst/extdata/cosmic_snv_signatures_v2.txt) and [version 3](https://github.com/mskcc/tempoSig/blob/master/inst/extdata/cosmic_sigProfiler_SBS_signatures.txt) tables of [COSMIC signature lists](https://cancer.sanger.ac.uk/cosmic/signatures) are included.

The "refitting" (as opposed to de novo discovery) of signature propotion solves the non-negative matrix factorization problem

    X = W * H
    
where **X** is the catalog matrix, **W** is the signature matrix (assumed to be known and fixed), and **H** is the exposure matrix of dimension (no. of reference signatures x no. of samples). The maximum likelihood estimation (MLE) algorithm formulates this problem in terms of a multinomial statistical model with observed counts **X** of categorical groups (trinucleotide contexts) mixed by given fixed proportions **W**. **tempoSig** uses the quasi-Newton [Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm](https://www.gnu.org/software/gsl/doc/html/multimin.html) for multi-dimensional optimization.
The output matrix is the transpose of **H**:

Tumor Sample Barcode | Signature.1 | Signature.2 | Signature.3 | Signature.4
-------------------- | ---- | ---- | ---- | ----
SAMPLE_1             | 0.37 | 0.01 | 0.00 | 0.00
SAMPLE_2             | 0.31 | 0.18 | 0.00 | 0.00

Each row is a vector of proportions that add up to 1.

To estimate p-values of significance of each proportion, the signature profile of each reference signature (columns of **W**) are randomly shuffled by permutation to sample the null distribution. The exposure vectors inferred from these null samples are compared with the observed vector from the original data, with the p-value defined as the fraction of null samples whose proprotions are higher than the observed values.

## Installation
Compilation requires GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/). In Ubuntu Linux,

    $ sudo apt-get install libgsl-dev
    
In OS-X,

    $ brew install gsl

Dependencies include [Rcpp](https://cran.r-project.org/package=Rcpp), [gtools](https://cran.r-project.org/package=gtools), [argparse](https://cran.r-project.org/package=argparse), and [coneproj](https://cran.r-project.org/package=coneproj). Installing **tempoSig** via

    > devtools::install_github("mskcc/tempoSig")

will also install dependencies.

## Quick start with command-line interface

### Main inference
If you are not interested in interactive usages with more flexibility and functionality, or want to use **tempoSig** as a part of a pipeline, use the command-line script [tempoSig.R](https://github.com/mskcc/tempoSig/blob/master/exec/tempoSig.R). If you installed **tempoSig** as an R package using `install_github`, find the path via

    > system.file('exec', 'tempoSig.R', package = 'tempoSig')
   
If you cloned the repository, the file is located at the `./exec` subdirectory of the github main directory. We denote this package directory path as `PKG_PATH`. The command syntax is

    $ $PKG_PATH/exec/tempoSig.R -h
     usage: ./tempoSig.R [-h]
                    [--cosmic_v2 | --cosmic_v3 | --cosmic_v3_SA | --cosmic_v3_exome]
                    [--sigfile SIGFILE]
                    [--pvalue]
                    [--nperm NPERM]
                    [--seed SEED]
                    [--pv.out PV.OUT] 
                    [--cbio]
                    CATALOG OUTPUT

     Fit mutational catalog to signatures

     positional arguments:
       CATALOG            input catalog data file
       OUTPUT             output file name

     optional arguments:
       -h, --help         show this help message and exit
       --cosmic_v2        use COSMIC v2 reference signatures (default)
       --cosmic_v3        use COSMIC v3 reference signatures
       --cosmic_v3_SA     use COSMIC v3 SigAnalyzer reference signatures
       --cosmic_v3_exome  use COSMIC v3 exome reference signatures
       --sigfile SIGFILE  custom input reference signature file; overrides
                          --cosmic_v2/3
       --pvalue           estimate p-values (default FALSE)
       --nperm NPERM      number of permutations for p-value estimation; default
                          1000
       --seed SEED        random number seed
       --pv.out PV.OUT    p-value output file
       --cbio             output in cBioPortal format (default FALSE)
     
Only two arguments are mandatory: `CATALOG` and `OUTPUT`, each specifying the paths of input catalog data and output file to be written. Both are tab-delimited text files with headers. See [tcga-brca_catalog.txt](https://github.com/mskcc/tempoSig/blob/master/inst/extdata/tcga-brca_catalog.txt) for a catalog file example. For instance,

    $ $PKG_PATH/exec/tempoSig.R $PKG_PATH/extdata/tcga-brca_catalog.txt output.txt
    
fits catalog data for 10 samples in `tcga-brca_catalog.txt` to [COSMIC v2 signatures](https://github.com/mskcc/tempoSig/edit/master/inst/extdata/cosmic_snv_signatures_v2.txt) (default). The output file `output.txt` has the following format:

Sample Name    | Number of Mutations | Signature.1 | Signature.2 | Signature.3 | Signature.4
-------------- | ------------------- | ----------- | ----------- | ----------- | -----------
TCGA.BH.A0EI   | 18                  | 0.61        | 0.01        | 0.00        | 0.00
TCGA.E9.A22B   | 50                  | 0.51        | 0.22        | 0.00        | 0.00
TCGA.OL.A5RV   | 10                  | 0.41        | 0.00        | 0.23        | 0.00 

The following will use the [COSMIC v3 signatures](https://github.com/mskcc/tempoSig/edit/master/inst/extdata/cosmic_sigProfiler_SBS_signatures.txt):

    $ $PKG_PATH/exec/tempoSig.R  --cosmic_v3 $PKG_PATH/extdata/tcga-brca_catalog.txt output_v3.txt

The output is similar, with the columns corresponding to 67 signatures:

Sample Name    | Number of Mutations | SBS.1.      | SBS.2       | SBS.3
-------------- | ------------------- | ----------- | ----------- | -----------
TCGA.BH.A0EI   | 18                  | 0.373       | 8.5e-3      | 0
TCGA.E9.A22B   | 50                  | 0.310       | 0.180       | 0
TCGA.OL.A5RV   | 10                  | 0.337       | 0           | 0 

One can use a custom reference signature list (in the same format as the default version 3 file) via the optional argument `--sigfile SIGFILE`.

### Catalog matrix generation

If a MAF file contains the column `Ref_Tri` [trinucleotide contexts surrounding the mutation site; use `make_trinuc_maf.py` script in [mutation-signatures](https://github.com/mskcc/mutation-signatures)], the catalog matrix can also be generated using [maf2cat()](https://github.com/mskcc/tempoSig/blob/master/man/maf2cat.Rd) or its command-line wrapper:

    $ ./maf2cat2.R -h
    usage: ./maf2cat2.R [-h] MAF CATALOG

    Construct mutational catalog from MAF file with Ref_Tri column

    positional arguments:
      MAF         input MAF file
      CATALOG     output catalog file

    optional arguments:
      -h, --help  show this help message and exit

### P-value estimation
Optionally, statistical significance of the set of proportions (rows in the exposure output) can be estimated by permutation sampling. For each signature, the exposure inference is repeated multiple times after permutation of the reference signature profile. P-values are the fractions of permuted replicates whose proportions (**H0**) are not lower than those of the original (**H1**). The p-value estimation is turned on by the argument `--pvalue`. The number of permutations is 1,000 by default and can be set with `--nperm NPERM`. The default output has the format:

Sample Name          | Number of Mutations | Signature.1.observed | Signature.1.pvalue   | Signature.2.observed | Signature.2.pvalue
-------------------- | ------------------- | -------------------- | -------------------- | -------------------- | ------------------
TCGA.BH.A0EI         | 18                  | 0.61                 | 0                    | 0.066                | 0.05              
TCGA.E9.A22B         | 50                  | 0.51                 | 0                    | 0.20                 | 0              
TCGA.OL.A5RV         | 10                  | 0.41                 | 0.6                  | 1.2e-11              | 0.55        

Note that p-value of 0 indicates that out of `NPERM` samples, none exceeded **H1**, and therefore must be interpreted as *P* < 1/`NPERM`. Alternatively, one can have two output files, one for exposure and the other for p-vaues, by specifiying the argument `--pv.out PV.OUT`. The exposure output `OUTPUT` is in the same format as that without p-value computation. The p-value output `PV.OUT` has the analogous format with columns for each signature p-values.

## Documentation

For more detailed interactive usage and extra options, see the [vignettes](http://htmlpreview.github.io/?https://github.com/mskcc/tempoSig/blob/master/old/tempoSig.html).

## Benchmark

In **Fig. 1**, the function [`simulateSpectra()`](man/simulateSpectra.Rd) was used with a breast cancer-like signature exposure profile (version 2) to generate simulated data of sample size 100 and mean TMB of 500 per genome. The **tempoSig** and [mutation-signatures](https://github.com/mskcc/mutation-signatures) were used to estimate exposures, which are identical. These estimates are broadly consistent with true exposure values. The accuracy improves with increasing mutation loads.

<figure>
<img src="https://github.com/mskcc/tempoSig/blob/master/old/msig_vs_cli.png" align="center" height="432" width="432"/>
    <figcaption> Fig. 1: Comparison of exposure fits from tempoSig and mutation-signatures. Data were simulated with breast cancer-like signatures. True exposure values are also shown along the y-axis. </figcaption>
</figure>

<br>
</br>

In **Fig. 2**, **tempoSig** and **mutation-signatures** running times were measured for varing sample sizes (fixed TMB of 500). The command-line version of **tempoSig** has overhead costs in speed compared to interactive R runs using the function [`extractSig()`](man/extractSig.Rd), which decreases in relative magnitude with increasing sample size. Overall, the efficiency gain of **tempoSig** is up to ~400x.

<br>

<figure>
<img src="https://github.com/mskcc/tempoSig/blob/master/old/time_cli.png" align="center" height="432" width="432"/>
    <figcaption> Fig. 2: Running time of tempoSig and mutation-signatures on data of varying sizes. Other conditions same as in Fig. 1. </figcaption>
</figure>

<br>

In **Fig. 3**, the overall accuracy of exposure proportions inferred from data sets simulated with breast cancer-like signature proportions were compared to true values using cosine similarity (higher the better; ranges from 0 to 1). Five other existing algorithms for refitting ([deconstructSigs](https://cran.r-project.org/package=deconstructSigs), [YAPSA](https://www.bioconductor.org/packages/YAPSA/), [MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns/), [MutationalCone](https://doi.org/10.1371/journal.pone.0221235), and [decompTumor2Sig](https://www.bioconductor.org/packages/decompTumor2Sig/)) were applied to the same data sets and their results compared to `tempoSig`. Although those more recent than `deconstructSigs`, one of the earliest refitting algorithms, exhibited similar performance, the maximum likelihood-based inference (`tempoSig` and `mutation-signatures`) consistently outperformed all others.

<br>
<figure>
<img src="https://github.com/mskcc/tempoSig/blob/master/old/cosim6.png" align="center" height="480" width="600"/>
    <figcaption> Fig. 3: Accuracy comparison of exposures predicted from simulated data of varying mutation loads with six refitting algorithms. </figcaption>
</figure>
