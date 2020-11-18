#!/usr/bin/env Rscript
library(maftools)

args = commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

var.annovar.maf <- annovarToMaf(annovar = input_file, Center = 'BCM', refBuild = 'GRCh38.p13', table = 'refGene')
write.table(var.annovar.maf,file=output_file,sep="\t",quote=FALSE,row.names=FALSE)
