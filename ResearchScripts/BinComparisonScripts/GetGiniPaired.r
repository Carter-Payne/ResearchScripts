#SCOPE pipeline, from the docs: https://bioconductor.org/packages/release/bioc/vignettes/SCOPE/inst/doc/SCOPE_vignette.html
args = commandArgs(trailingOnly=TRUE)#args[1]=bin size
library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)

#bamfolder <- system.file("extdata", package = "WGSmapp")
bamfolder <- /path/to/bam/folder
bamFile <- list.files(bamfolder, pattern = '*.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = "hg38", resolution=as.integer(args[1]))
ref_raw <- bambedObj$ref

mapp <- get_mapp(ref_raw, hgref = "hg38")
#head(mapp)

gc <- get_gc(ref_raw, hgref = "hg38")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
#ref_raw

coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = 'paired-end', hgref = "hg38")
Y_raw <- coverageObj$Y

QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = Y_raw, sampname_raw = sampname_raw, ref_raw = ref_raw, QCmetric_raw = QCmetric_raw)

Y <- qcObj$Y

Gini <- get_gini(Y)
file <- paste(args[2],"gini.csv", sep="")
write.table( x = data.frame(Gini), file = file, sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

