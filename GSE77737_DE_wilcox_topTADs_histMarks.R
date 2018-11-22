SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

vview <- function(dt) {
  if(nrow(dt) >= 5 & ncol(dt) >= 5)
    dt[1:5,1:5]
}

startTime <- Sys.time()

cat("START> Rscript GSE77737_DE_histMarks.R\n")
# Rscript GSE77737_DE_histMarks.R <histMark> <signalType>
# Rsscript GSE77737_DE_histMarks.R <histMark> <signalType>
# Rscript GSE77737_DE_histMarks.R H3K27ac mean
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
# pvalThresh <- 0.05
# DESeq_fitType <- "local"
# # I used local because I got the message
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# and "local" is implemented in both DESeq and DESeq2
####################################################

outFold <- file.path(
  setDir,
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/GSE77737_DE_wilcox_topTADs_HISTMARKS",
  curr_dataset
)
system(paste0("mkdir -p ", outFold))

histMark <- "H3K27ac"
signalType <- "mean"

# Rscript GSE77737_DE_wilcox_topTADs_histMarks.R H3K27ac mean
# Rscript GSE77737_DE_wilcox_topTADs_histMarks.R <histMark> <signalType>
args <- commandArgs(trailingOnly = TRUE)
histMark <- args[1]
signalType <- args[2]

logFile <- file.path(outFold, paste0(histMark, "_", signalType, "_", cond1, "_", cond2, "_wilcoxDE_logFile.txt"))
system(paste0("rm -f ", logFile))

# txt <- paste0("> SOME SETTINGS:\n")
# txt <- pastse0("... pvalThresh = ", pvalThresh, "\n")
# txt <- pastse0("... local = ", DESeq_fitType, "\n")

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

# iterate over the topTADs and perform wilcox tests TAD by TAD



# go TAD by TAD, extract the peaks that are within the TAD
# for each peak -> wilcoxon Test
pvalFile <- file.path(setDir,
                      "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/TCGAcrc_msi_mss/11_runEmpPvalCombined",
                       "emp_pval_combined.Rdata")
stopifnot(file.exists(pvalFile))
pvalComb <- eval(parse(text = load(pvalFile)))
adj_pvalComb <- p.adjust(pvalComb, method="BH")
adj_pvalComb_sorted <- sort(adj_pvalComb)
adj_pvalComb_rank <- rank(adj_pvalComb, ties="min")

topTADs <- names(adj_pvalComb_sorted)[1:10]

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TADpos_DT <- read.delim(TADpos_file, header=F,col.names=c("chromo", "region", "start", "end"))


for(topTAD in topTADs) {
  txt <- paste0("> ", histMark, " - ", signalType, "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... topTAD:\t", topTAD, "\n")
  printAndLog(txt, logFile)
  
  tad_start <- TADpos_DT$start[TADpos_DT$region == topTAD]
  tad_end <- TADpos_DT$end[TADpos_DT$region == topTAD]
  
  tad_peaks <- peakCoord_DT$peak_name[peakCoord_DT$peak_start >= tad_start &
                                        peakCoord_DT$peak_end <= tad_end
                                        ]
  txt <- paste0("...... # peaks:\t", length(tad_peaks), "\n")
  printAndLog(txt, logFile)
  
  if(length(tad_peaks) == 0) {
    wt <- wilcox.test(unlist(tad_countsMat[,cond1cols]), unlist(tad_countsMat[,cond2cols]))
    txt <- paste0("...... Wilcox's p-val:\t", "NA", "\n")
    printAndLog(txt, logFile)
  } else{
    tad_countsMat <- peakCountMatrix[tad_peaks,]
    cond1cols <- colnames(peakCountMatrix)[grep(cond1, colnames(peakCountMatrix))]
    cond2cols <- colnames(peakCountMatrix)[grep(cond2, colnames(peakCountMatrix))]
    wt <- wilcox.test(unlist(tad_countsMat[,cond1cols]), unlist(tad_countsMat[,cond2cols]))
    
    txt <- paste0("...... Mean ", cond1, "\t=\t", sprintf("%.2f", mean(unlist(tad_countsMat[,cond1cols]), na.rm=T)), "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... Mean ", cond2, "\t=\t", sprintf("%.2f", mean(unlist(tad_countsMat[,cond2cols]), na.rm=T)), "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... Wilcox's p-val:\t", sprintf("%2.2e",wt$p.value), "\n")
    printAndLog(txt, logFile)
  }
  

  
}

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
