#!/usr/bin/env Rscript
#Takes arg[1] as the name of the CN and Y files, and arg[2] as the gini coefficient produced by GetGiniSingle.r.
args = commandArgs(trailingOnly=TRUE)
library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
bamfolder <- /folder/to/bam/files
bamFile <- list.files(/folder/to/bam/files, pattern = '*.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = "hg38",resolution=as.integer(args[1]))
ref_raw <- bambedObj$ref
mapp <- get_mapp(ref_raw, hgref = "hg38")
gc <- get_gc(ref_raw, hgref = "hg38")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = 'single-end', hgref = "hg38")
Y_raw <- coverageObj$Y
QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = Y_raw, sampname_raw = sampname_raw, ref_raw = ref_raw, QCmetric_raw = QCmetric_raw)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric
Gini <- get_gini(Y)
normObj <- normalize_codex2_ns_noK(Y_qc = Y,gc_qc = ref$gc,norm_index = which(Gini<=as.numeric(args[2])))
ploidy <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref)
normObj.scope <- normalize_scope_foreach(Y_qc = Y, gc_qc = ref$gc,K = 1, ploidyInt = ploidy,norm_index = which(Gini<=as.numeric(args[2])), T = 1:5,beta0 = normObj$beta.hat)
Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
fGC.hat <- normObj.scope$fGC.hat[[which.max(normObj.scope$BIC)]]
Ploidies <- initialize_ploidy(Y, Yhat, ref, maxPloidy = 6, minPloidy = 1.5, minBinWidth = 5, SoS.plot = FALSE)
chrs <- unique(as.character(seqnames(ref)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs
for (chri in chrs) {
	message('\n', chri, '\n')
	segment_cs[[chri]] <- segment_CBScs(Y = Y,Yhat = Yhat,sampname = colnames(Y),ref = ref,chr = chri,mode = "integer", max.ns = 1)
}
iCN <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
EndCN <- paste(args[1],'_CN.txt',sep="")
EndY <- paste(args[1],'_Y.txt',sep="")
write.csv(iCN, EndCN)
write.csv(Y, EndY)

