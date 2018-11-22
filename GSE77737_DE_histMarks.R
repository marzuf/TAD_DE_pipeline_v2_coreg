SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

vview <- function(dt) {
  if(nrow(dt) >= 5 & ncol(dt) >= 5)
    dt[1:5,1:5]
}

startTime <- Sys.time()

cat("START> Rscript GSE77737_DE_histMarks.R\n")
# Rscript GSE77737_DE_histMarks.R <histMark> <signalType>

# data are stored in:
# /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss/H3K27ac/H3K27ac_MSI_countMatrix.Rdata

printAndLog <- function(mytxt, myfile="") {
  cat(mytxt)
  cat(mytxt, file = myfile, append = T)
}

## HARD-CODED ######################################
curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"
chromoLevels <- paste0("chr", c(1:22, "X"))
pvalThresh <- 0.05
DESeq_fitType <- "local"
# I used local because I got the message
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# and "local" is implemented in both DESeq and DESeq2
####################################################

outFold <- file.path(
  setDir,
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/GSE77737_DE_HISTMARKS",
  curr_dataset
)
system(paste0("mkdir -p ", outFold))

histMark <- "H3K27ac"
signalType <- "mean"

# Rscript GSE77737_DE_histMarks.R H3K27ac mean
# Rscript GSE77737_DE_histMarks.R <histMark> <signalType>
args <- commandArgs(trailingOnly = TRUE)
histMark <- args[1]
signalType <- args[2]

logFile <- file.path(outFold, paste0(histMark, "_", signalType, "_", cond1, "_", cond2, "_logFile.txt"))
system(paste0("rm -f ", logFile))

# HARD CODED MAIN FOLDER
mainFold <- file.path(setDir, 
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss",
  histMark)

cond1DT <- eval(parse(text = load(
      file.path(mainFold, paste0(histMark, "_", cond1, "_countMatrix.Rdata")))))
vview(cond1DT)
nrow(cond1DT)

cond2DT <- eval(parse(text = load(
  file.path(mainFold, paste0(histMark, "_", cond2, "_countMatrix.Rdata")))))
vview(cond2DT)
nrow(cond2DT)

tmp1dt <- cond1DT[, c("peak_name", "peak_chromo", "peak_start", "peak_end")]
tmp2dt <- cond2DT[, c("peak_name", "peak_chromo", "peak_start", "peak_end")]
stopifnot(all.equal(tmp1dt, tmp2dt))
peakCoord_DT <- tmp1dt
stopifnot(!any(duplicated(peakCoord_DT$peak_name)))
rownames(peakCoord_DT) <- peakCoord_DT$peak_name

signal_cond1DT <- cond1DT[, grep(paste0("_", signalType, "$"), colnames(cond1DT))]
vview(signal_cond1DT)
signal_cond2DT <- cond2DT[, grep(paste0("_", signalType, "$"), colnames(cond2DT))]
vview(signal_cond2DT)
colnames(signal_cond1DT) <- gsub(paste0("_", signalType), "", colnames(signal_cond1DT))
colnames(signal_cond2DT) <- gsub(paste0("_", signalType), "", colnames(signal_cond2DT))
rownames(signal_cond1DT) <- cond1DT$peak_name
rownames(signal_cond2DT) <- cond2DT$peak_name

peakSignalCount_DT <- merge(peakCoord_DT, signal_cond1DT, by="row.names")
rownames(peakSignalCount_DT) <- peakSignalCount_DT$Row.names
stopifnot(rownames(peakSignalCount_DT) == peakSignalCount_DT$peak_name )
peakSignalCount_DT$Row.names <- NULL
peakSignalCount_DT <- merge(peakSignalCount_DT, signal_cond2DT, by="row.names")
vview(peakSignalCount_DT)
rownames(peakSignalCount_DT) <- peakSignalCount_DT$Row.names
stopifnot(rownames(peakSignalCount_DT) == peakSignalCount_DT$peak_name )
peakSignalCount_DT$Row.names <- NULL
head(peakSignalCount_DT)

stopifnot(ncol(peakSignalCount_DT) == (ncol(peakCoord_DT) + ncol(signal_cond1DT) + ncol(signal_cond2DT)) )
stopifnot(nrow(peakSignalCount_DT) == nrow(peakCoord_DT))
stopifnot(nrow(peakSignalCount_DT) == nrow(signal_cond1DT))
stopifnot(nrow(peakSignalCount_DT) == nrow(signal_cond2DT)) 

peakSignalCount_DT$peak_name <- as.character(peakSignalCount_DT$peak_name)
peakSignalCount_DT$peak_chromo <- as.character(peakSignalCount_DT$peak_chromo)
peakSignalCount_DT$peak_start <- as.numeric(as.character(peakSignalCount_DT$peak_start))
peakSignalCount_DT$peak_end <- as.numeric(as.character(peakSignalCount_DT$peak_end))
stopifnot(!any(is.na(peakSignalCount_DT)))
head(peakSignalCount_DT,2)
peakSignalCount_DT$peak_chromo <- factor(peakSignalCount_DT$peak_chromo, levels = chromoLevels[chromoLevels %in% peakSignalCount_DT$peak_chromo])
peakSignalCount_DT <- peakSignalCount_DT[order(as.numeric(peakSignalCount_DT$peak_chromo), peakSignalCount_DT$peak_start, peakSignalCount_DT$peak_end),]

peakCountMatrix <- peakSignalCount_DT
rownames(peakCountMatrix) <- peakCountMatrix$peak_name
peakCountMatrix$peak_name <- NULL
peakCountMatrix$peak_start <- NULL
peakCountMatrix$peak_end <- NULL
peakCountMatrix$peak_chromo <- NULL

cond1_samples <- colnames(peakCountMatrix)[grepl(paste0("_", cond1), colnames(peakCountMatrix))]
stopifnot(length(cond1_samples) > 0)
cond2_samples <- colnames(peakCountMatrix)[grepl(paste0("_", cond2), colnames(peakCountMatrix))]
stopifnot(length(cond2_samples) > 0)
stopifnot( (length(cond1_samples) + length(cond2_samples)) == ncol(peakCountMatrix)) 

peakCountMatrix <- peakCountMatrix[,c(cond1_samples, cond2_samples)]
myconditions <- factor(c(rep(cond1, length(cond1_samples)), rep(cond2, length(cond2_samples))))
sampleConditions <- setNames(as.character(myconditions), c(cond1_samples, cond2_samples))
  
stopifnot( (length(cond1_samples) + length(cond2_samples)) == ncol(peakCountMatrix)) 

vview(peakCountMatrix)
int_peakCountMatrix <- round(peakCountMatrix)
vview(int_peakCountMatrix)

############## TRY DESeq
library(DESeq)
# needs integer !!!
chipseq_cds_DESeq <- DESeq::newCountDataSet(int_peakCountMatrix, myconditions)

# 1) normalization
chipseq_cds_DESeq <- DESeq::estimateSizeFactors(chipseq_cds_DESeq)  # estimate size factors from the count data
sizeFactors(chipseq_cds_DESeq)

# if we divide each column by size factor, counts are brought to common scale, this is done in:
head(DESeq::counts(chipseq_cds_DESeq, normalized=T))

# 2) variance estimation
chipseq_cds_DESeq <- DESeq::estimateDispersions(chipseq_cds_DESeq, fitType = DESeq_fitType)  # estimate relationship between variance and mean (dispersion=square of coeff. of biological variation)
# variance of counts as sum of 1) sample-to-sample variation (dominates for highly expressed genes) 2) uncertainty counting reads (shot noise or Poisson noise)
# the function performs 3 steps 1) estimates dispersion value for each gene, 2) fits curve through estimate, 3) assign to each dipsersion value (choice between per-gene estimated or fitted value)
plotDispEsts(chipseq_cds_DESeq)
# looks similar to the plot of the vignette

# 3) calling differential expression
chipseq_DESeq_DE_results <- DESeq::nbinomTest(cds = chipseq_cds_DESeq, condA = cond1, condB = cond2)
head(chipseq_DESeq_DE_results)
# id  baseMean baseMeanA baseMeanB foldChange log2FoldChange      pval padj
# 1 peak1 5.9507654  7.727674 5.4661538  0.7073479   -0.499508222 0.3241936    1
# 2 peak2 1.4491898  1.463575 1.4452667  0.9874910   -0.018160455 1.0000000    1
# 3 peak3 0.9047068  1.104052 0.8503398  0.7701989   -0.376697115 0.9732335    1

# id=peak name
# baseMean=mean normalised counts all samples 
# baseMeanA  = mean sample condition A
# baseMeanB  = mean sample condition B
# foldChange, log2FoldChange = FC from condition A to condition B (cond. B over cond. A -> cond2 over cond1)
# pval, padj = pval and BH adj. pval

txt <- paste0("... DESeq - range pval:\t", paste0(round(range(chipseq_DESeq_DE_results$pval), 4), collapse=" - "), "\n")
printAndLog(mytxt = txt, myfile = logFile)
txt <- paste0("... DESeq - range adj. pval:\t", paste0(round(range(chipseq_DESeq_DE_results$padj), 4), collapse=" - "), "\n")
printAndLog(mytxt = txt, myfile = logFile)
txt <- paste0("... DESeq - number of DE peaks identified (adj. pval thresh =", pvalThresh, "):\t", sum(chipseq_DESeq_DE_results$padj <= pvalThresh), "/", nrow(chipseq_DESeq_DE_results), "\n" )
printAndLog(mytxt = txt, myfile = logFile)

# ... range pval:	 0 - 1 
# ... range adj. pval:	 0.2289 - 1 
# ... number of DE peaks identified with DESeq (adj. pval thresh = 0.05 ):	 0 / 169760 

txt <- paste0("... DESeq - Top 5 peaks - pval\n")
printAndLog(mytxt = txt, myfile = logFile)
write.table(chipseq_DESeq_DE_results[order(chipseq_DESeq_DE_results$pval),][1:5,], row.names=F, col.names=T, file ="", quote=F, sep="\t")
write.table(chipseq_DESeq_DE_results[order(chipseq_DESeq_DE_results$pval),][1:5,], row.names=F, col.names=T, file =logFile, quote=F, sep="\t", append=T)
txt <- paste0("\n")
printAndLog(mytxt = txt, myfile = logFile)
txt <- paste0("... DESeq - Top 5 peaks - adj. pval\n")
printAndLog(mytxt = txt, myfile = logFile)
write.table(chipseq_DESeq_DE_results[order(chipseq_DESeq_DE_results$padj),][1:5,], row.names=F, col.names=T, file ="", quote=F, sep="\t")
write.table(chipseq_DESeq_DE_results[order(chipseq_DESeq_DE_results$padj),][1:5,], row.names=F, col.names=T, file =logFile, quote=F, sep="\t", append =  T)

topPeaks_DESeq_pval <- chipseq_DESeq_DE_results[order(chipseq_DESeq_DE_results$pval),][1:5,"id"]
topPeaks_DESeq_adjPval <- chipseq_DESeq_DE_results[order(chipseq_DESeq_DE_results$padj),][1:5,"id"]

rm(chipseq_DESeq_DE_results)

############## TRY DESeq2
detach("package:DESeq", unload=TRUE)
library(DESeq2)

DESeq2_colData <- data.frame(sample = names(sampleConditions), condition = sampleConditions, stringsAsFactors = FALSE)
rownames(DESeq2_colData) <- DESeq2_colData$sample
DESeq2_colData <- DESeq2_colData[,c("condition"), drop=F]
# comparison will be last level over reference revel (-> cond2 over cond1)
DESeq2_colData$condition <- factor(DESeq2_colData$condition, levels = c(cond1, cond2))
# needs integer !!!
chipseq_cds_DESeq2 <- DESeq2::DESeqDataSetFromMatrix(countData = int_peakCountMatrix,
                                                     colData = DESeq2_colData,
                                                     design = ~condition)
# differential expression analysis
# DESeq performs a default analysis through the steps:
# 1) estimation of size factors: estimateSizeFactors, 2) estimation of dispersion: estimateDispersions, 3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
chipseq_DESeq2_DE <- DESeq2::DESeq(chipseq_cds_DESeq2, fitType = DESeq_fitType)
chipseq_DESeq2_DE_results <- results(chipseq_DESeq2_DE)
head(chipseq_DESeq2_DE_results)
# his object contains the results columns: baseMean, log2FoldChange, lfcSE, stat, pvalue and padj, and also includes metadata columns of variable information. 

# lfcSE -> standard error of the log2FoldChange
#               - For the Wald test, stat is the Wald statistic: the log2FoldChange divided by lfcSE, which is compared to a standard Normal distribution to generate a two-tailed pvalue. 
#               - For the likelihood ratio test (LRT), stat is the difference in deviance between the reduced model and the full model, 
#                 which is compared to a chi-squared distribution to generate a pvalue.

# log2 fold change (MLE): condition MSS vs MSI 
# Wald test p-value: condition MSS vs MSI 
chipseq_DESeq2_DE_results <- as.data.frame(chipseq_DESeq2_DE_results)
# there are some NA
chipseq_DESeq2_DE_results <- chipseq_DESeq2_DE_results[!is.na(chipseq_DESeq2_DE_results$pvalue),]
txt <- paste0("... DESeq2 - range pval:\t", paste0(round(range(chipseq_DESeq2_DE_results$pvalue), 4), collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... DESeq2 - range adj. pval:\t", paste0(round(range(chipseq_DESeq2_DE_results$padj), 4), collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... DESeq2 - number of DE peaks identified (adj. pval thresh =", pvalThresh, "):\t", sum(chipseq_DESeq2_DE_results$padj <= pvalThresh), "/", nrow(chipseq_DESeq2_DE_results), "\n" )
printAndLog(mytxt=txt, myfile=logFile)
# ... DESeq2 - range pval:	 0 - 1 
# ... DESeq2 - range adj. pval:	 0.9372 - 1 
# ... DESeq2 - number of DE peaks identified (adj. pval thresh = 0.05 ):	 0 / 169757 

txt <- paste0("... DESeq2 - Top 5 peaks - pval\n")
printAndLog(mytxt=txt, myfile=logFile)
write.table(chipseq_DESeq2_DE_results[order(chipseq_DESeq2_DE_results$pval),][1:5,], row.names=T, col.names=T, file ="", quote=F, sep="\t")
write.table(chipseq_DESeq2_DE_results[order(chipseq_DESeq2_DE_results$pval),][1:5,], row.names=T, col.names=T, file =logFile, append=T, quote=F, sep="\t")
txt <- paste0("\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... DESeq2 - Top 5 peaks - adj. pval\n")
printAndLog(mytxt=txt, myfile=logFile)
write.table(chipseq_DESeq2_DE_results[order(chipseq_DESeq2_DE_results$padj),][1:5,], row.names=T, col.names=T, file ="", quote=F, sep="\t")
write.table(chipseq_DESeq2_DE_results[order(chipseq_DESeq2_DE_results$padj),][1:5,], row.names=T, col.names=T, file =logFile, append=T, quote=F, sep="\t")
txt <- paste0("\n")
printAndLog(mytxt=txt, myfile=logFile)

topPeaks_DESeq2_pval <- rownames(chipseq_DESeq2_DE_results[order(chipseq_DESeq2_DE_results$pval),][1:5,])
topPeaks_DESeq2_adjPval <- rownames(chipseq_DESeq2_DE_results[order(chipseq_DESeq2_DE_results$padj),][1:5,])

rm(chipseq_DESeq2_DE_results)

############## TRY edgeR
detach("package:DESeq2", unload=TRUE)
library(edgeR)

edgeR_colData <- DESeq2_colData
colnames(edgeR_colData) <- "group"

chipseq_cds_edgeR <- DGEList(counts = peakCountMatrix,
                             samples = edgeR_colData) 

stopifnot(chipseq_cds_edgeR$samples$group.1 == chipseq_cds_edgeR$samples$group)
stopifnot(levels(chipseq_cds_edgeR$samples$group.1) == levels(chipseq_cds_edgeR$samples$group))
chipseq_cds_edgeR$samples$group.1 <- NULL

# store normalization factors
chipseq_cds_edgeR <- calcNormFactors(chipseq_cds_edgeR)
plotMDS(chipseq_cds_edgeR, col = as.numeric(chipseq_cds_edgeR$samples$group))

# skip the filter step for the moment
      # chipseq_cds_edgeR_cpm <- cpm(chipseq_cds_edgeR)
      # k <- rowSums(chipseq_cds_edgeR_cpm > 1 ) >= 2
      # chipseq_cds_edgeR_filter <- chipseq_cds_edgeR[k,]
      # chipseq_cds_edgeR_cpm <- cpm(chipseq_cds_edgeR_filter)
mydesign <- model.matrix(~0+group, data = chipseq_cds_edgeR$samples)
colnames(mydesign) <- levels(chipseq_cds_edgeR$samples$group)

# calculate dispersion estimates
chipseq_cds_edgeR <- estimateGLMCommonDisp(chipseq_cds_edgeR, mydesign)
chipseq_cds_edgeR <- estimateGLMTrendedDisp(chipseq_cds_edgeR, mydesign)
chipseq_cds_edgeR <- estimateGLMTagwiseDisp(chipseq_cds_edgeR, mydesign)

# fit the model
chipseq_cds_edgeR_fit <- glmFit(chipseq_cds_edgeR, mydesign)

plotBCV(chipseq_cds_edgeR)
# !!! very weird plot !!!

# conduct the LKHD ratio test
chipseq_cds_edgeR_DE <- glmLRT(chipseq_cds_edgeR_fit, coef=2)
chipseq_cds_edgeR_DE_results <- as.data.frame(topTags(chipseq_cds_edgeR_DE, n=Inf))

txt <- paste0("... edgeR - range pval:\t", paste0(round(range(chipseq_cds_edgeR_DE_results$PValue), 4), collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... edgeR - range adj. pval:\t", paste0(round(range(chipseq_cds_edgeR_DE_results$FDR), 4), collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... edgeR - number of DE peaks identified (adj. pval thresh =", pvalThresh, "):\t", sum(chipseq_cds_edgeR_DE_results$FDR <= pvalThresh), "/", nrow(chipseq_cds_edgeR_DE_results), "\n" )
printAndLog(mytxt=txt, myfile=logFile)

# ... edgeR - range pval:	 0 - 1 
# ... edgeR - range adj. pval:	 0 - 1 
# ... edgeR - number of DE peaks identified (adj. pval thresh = 0.05 ):	 166096 / 169757 

txt <- paste0("... edgeR - Top 5 peaks - pval\n")
printAndLog(mytxt=txt, myfile=logFile)
write.table(chipseq_cds_edgeR_DE_results[order(chipseq_cds_edgeR_DE_results$PValue),][1:5,], row.names=T, col.names=T, file ="", quote=F, sep="\t")
write.table(chipseq_cds_edgeR_DE_results[order(chipseq_cds_edgeR_DE_results$PValue),][1:5,], row.names=T, col.names=T, file =logFile, append=T, quote=F, sep="\t")
txt <- paste0("\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... edgeR - Top 5 peaks - adj. pval\n")
printAndLog(mytxt=txt, myfile=logFile)
write.table(chipseq_cds_edgeR_DE_results[order(chipseq_cds_edgeR_DE_results$FDR),][1:5,], row.names=T, col.names=T, file ="", quote=F, sep="\t")
write.table(chipseq_cds_edgeR_DE_results[order(chipseq_cds_edgeR_DE_results$FDR),][1:5,], row.names=T, col.names=T, file =logFile, append=T, quote=F, sep="\t")

topPeaks_edgeR_pval <- rownames(chipseq_cds_edgeR_DE_results[order(chipseq_cds_edgeR_DE_results$PValue),][1:5,])
topPeaks_edgeR_adjPval <- rownames(chipseq_cds_edgeR_DE_results[order(chipseq_cds_edgeR_DE_results$FDR),][1:5,])

rm(chipseq_cds_edgeR_DE_results)

txt <- paste0("> top 5 peaks pval - DESeq:\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("..." , paste0(topPeaks_DESeq_pval, collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... top 5 peaks adj. pval - DESeq:\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("..." , paste0(topPeaks_DESeq_adjPval, collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... top 5 peaks pval - DESeq2:\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("..." , paste0(topPeaks_DESeq2_pval, collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... top 5 peaks adj. pval - DESeq2:\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("..." , paste0(topPeaks_DESeq2_adjPval, collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... top 5 peaks pval - edgeR:\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("..." , paste0(topPeaks_edgeR_pval, collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("... top 5 peaks adj. pval - edgeR:\n")
printAndLog(mytxt=txt, myfile=logFile)
txt <- paste0("..." , paste0(topPeaks_edgeR_adjPval, collapse=" - "), "\n")
printAndLog(mytxt=txt, myfile=logFile)


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
