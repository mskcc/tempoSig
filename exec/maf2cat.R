#!/usr/bin/env Rscript
library('argparse')
suppressMessages(library('BSgenome.Hsapiens.UCSC.hg19'))
parser <- ArgumentParser(description = 'Construct mutational catalog from MAF file')
parser$add_argument('maf', metavar = 'MAF',
                    help = 'input MAF file')
parser$add_argument('catalog', metavar = 'CATALOG', 
                    help = 'output catalog file')
args <- parser$parse_args()

maf <- maftools::read.maf(args$maf)
catalog <- t(maftools::trinucleotideMatrix(maf = maf, prefi = 'chr', add = TRUE, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')$nmf_matrix)
write.table(catalog, file = args$catalog, quote = F, sep = '\t')
