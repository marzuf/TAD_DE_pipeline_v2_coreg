library(foreach)
library(doMC)
library(dplyr)
library(reshape2)

source("coreg_utils.R")

# Rscript histpeaks_ratios_topTADs.R

# in this version assign mRNA to TAD based on entrezID

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 20))

if(SSHFS) setwd(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"))

curr_dataset <- "TCGAcrc_msi_mss"
histMark <- "H3K27ac"

cond1 <- "MSI"
cond2 <- "MSS"

outFold <- file.path(setDir, 
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg",
                     "HISTPEAKS_STATS",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 500, 10)

all_hist_marks <- c("H3K27ac", "H3K4me1")

all_peak_size <- foreach(histMark = all_hist_marks) %do% {
  
  cat(paste0("> START: ", histMark, "\n"))
  # HARD CODED BINARY MATRIX
  cat("... load binary matrix\n")
  binaryMatrixF <- file.path(setDir, 
                             "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_BINARYMATRIX",
                             curr_dataset, 
                             histMark, 
                             paste0(histMark, "_", cond1, "_", cond2, "_binaryMatrix.Rdata"))
  stopifnot(file.exists(binaryMatrixF))
  binaryMatrix <- eval(parse(text = load(binaryMatrixF))) 
  head(binaryMatrix)
  
  binaryMatrix$peak_name <- as.character(binaryMatrix$peak_name)
  binaryMatrix$peak_chromo <- as.character(binaryMatrix$peak_chromo)
  binaryMatrix$peak_start <- as.numeric(as.character(binaryMatrix$peak_start))
  binaryMatrix$peak_end <- as.numeric(as.character(binaryMatrix$peak_end))
  
  binaryMatrix <- binaryMatrix[order(binaryMatrix$peak_chromo, binaryMatrix$peak_start, binaryMatrix$peak_end),]
  binaryMatrix$peak_size <- binaryMatrix$peak_end - binaryMatrix$peak_start + 1

  log10(binaryMatrix$peak_size)
  
} # end-iterate over histone marks

outFile <- file.path(outFold, paste0("density_histPeaks_size.",plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(all_peak_size, my_xlab = "peak size (log10 bp)", plotTit=paste0("Histone peak size"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

