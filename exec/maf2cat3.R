#!/usr/bin/env Rscript
# Convert maf file without Ref_Tri column into catalog (assumes hg19)
library('argparse')
library('tempoSig')
library('BSgenome.Hsapiens.UCSC.hg19')
#maf2cat3 <- function(maf, ref.genome, fix.chr = TRUE, progress.bar = TRUE){}

parser <- ArgumentParser(description = 'Construct mutational catalog from MAF file with Ref_Tri column')
parser$add_argument('maf', metavar = 'MAF', help = 'input MAF file')
parser$add_argument('catalog', metavar = 'CATALOG', help = 'output catalog file')
args <- parser$parse_args()

catalog <- maf2cat3(maf = args$maf, ref.genome=BSgenome.Hsapiens.UCSC.hg19)
write.table(catalog, file = args$catalog, quote = F, sep = '\t')
