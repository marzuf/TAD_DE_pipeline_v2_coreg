SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

library(foreach)
library(doMC)

registerDoMC(ifelse(SSHFS, 2, 20))

vview <- function(dt) {
  if(nrow(dt) >= 5 & ncol(dt) >= 5)
    dt[1:5,1:5]
}

startTime <- Sys.time()
 
cat("START> Rscript GSE77737_binaryMatrix_histMarks.R\n")
# Rscript GSE77737_binaryMatrix_histMarks.R <histMark> 

# use the input of build_allMergedPeaks_binaryMatrix.R
#(script located in /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017)
# -> input data are stored in:
# /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss/H3K27ac/H3K27ac_MSI_MSS_binaryMatrix.Rdata

## HARD-CODED ######################################
curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"
chromoLevels <- paste0("chr", c(1:22, "X"))

# HARD CODED NUMBER OF FILES TO CHECK
H3K27ac_nCheckCond1 <- 6
H3K4me1_nCheckCond1 <- 4
H3K27ac_nCheckCond2 <- 22
H3K4me1_nCheckCond2 <- 23

####################################################

outFold <- file.path(
  setDir,
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/GSE77737_BINARYMATRIX_HISTMARKS",
  curr_dataset
)
system(paste0("mkdir -p ", outFold))

histMark <- "H3K27ac"
# Rscript GSE77737_binaryMatrix_histMarks.R H3K27ac 
# Rscript GSE77737_binaryMatrix_histMarks.R <histMark> 
args <- commandArgs(trailingOnly = TRUE)
histMark <- args[1]

# HARD CODED MAIN FOLDER
mainFold <- file.path(setDir, 
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_BINARYMATRIX",
  curr_dataset, 
  histMark)
# //mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_BINARYMATRIX/TCGAcrc_msi_mss/H3K27ac/H3K27ac_MSI_MSS_binaryMatrix.Rdata

binaryMat <-  eval(parse(text = load(
  file.path(mainFold, paste0(histMark, "_", cond1, "_", cond2, "_binaryMatrix.Rdata")))))

assignedPeaksFile <- file.path(setDir,
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017",
  histMark, 
  paste0(histMark, "_mergedPeaks_assignedRegions.txt"))
stopifnot(file.exists(assignedPeaksFile))

peak2tad_DT <- read.delim(assignedPeaksFile, col.names=c("chromo", "start", "end", "peak", "region"), stringsAsFactors = FALSE)
peak2tad_DT <- peak2tad_DT[!is.na(peak2tad_DT$region),]
stopifnot(!any(is.na(peak2tad_DT)))

topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region

### FIRST TEST: TEST FOR EACH PEAK
peakTests_DT <- foreach(i_peak = seq_len(nrow(binaryMat)), .combine='rbind') %dopar% {
  count1_cond1 <- sum( binaryMat[i_peak, grepl(paste0("_", cond1), colnames(binaryMat))] == 1)
  count0_cond1 <- sum( binaryMat[i_peak, grepl(paste0("_", cond1), colnames(binaryMat))] == 0)
  count1_cond2 <- sum( binaryMat[i_peak, grepl(paste0("_", cond2), colnames(binaryMat))] == 1)
  count0_cond2 <- sum( binaryMat[i_peak, grepl(paste0("_", cond2), colnames(binaryMat))] == 0)
  stopifnot( (count1_cond1 + count0_cond1) == eval(parse(text = paste0(histMark, "_", "nCheckCond1"))))
  stopifnot( (count1_cond2 + count0_cond2) == eval(parse(text = paste0(histMark, "_", "nCheckCond2"))))
  testMatrix <- matrix(c(count1_cond1, count1_cond2, count0_cond1, count0_cond2), byrow=T, nrow=2)  
  FT <- fisher.test(testMatrix, alternative="greater")
  # same as
  # testMatrix <- matrix(c(count1_cond1, count0_cond1, count1_cond2, count0_cond2), byrow=T, nrow=2)  
  # fisher.test(testMatrix, alternative="greater")
  data.frame(
    peak_chromo = binaryMat[i_peak,"peak_chromo"],
    peak_start = binaryMat[i_peak,"peak_start"],
    peak_end = binaryMat[i_peak,"peak_end"],
    peak_name = binaryMat[i_peak,"peak_name"],
    FT_pval = round(as.numeric(FT$p.value), 4), 
    FT_or = round(as.numeric(FT$estimate),4),
    nPeaks_cond1 = count1_cond1,
    nPeaks_cond2 = count1_cond2,
    ratioPeaks_cond1 = round(count1_cond1/(count1_cond1+count0_cond1), 4),
    ratioPeaks_cond2 = round(count1_cond2/(count1_cond2+count0_cond2), 4),
    stringsAsFactors = F
  )
}
rownames(peakTests_DT) <- NULL
colnames(peakTests_DT) <- gsub("cond1", cond1, colnames(peakTests_DT))
colnames(peakTests_DT) <- gsub("cond2", cond2, colnames(peakTests_DT))
# stopifnot(nrow(peakTests_DT) == nrow(binaryMat))

outFile <- file.path(outFold, paste0(histMark, "_", cond1, "_", cond2, "_", "peakTests_DT.Rdata"))
save(peakTests_DT, file = outFile)
if(file.exists(outFile)) {
  cat(paste0("... written: ", outFile, "\n"))
} else {
  stop("error\n")
}

### SECOND TEST: ITERATE OVER topTADs
topTADs_histPeaks_FT <- foreach(tad = topTADs, .combine='rbind') %dopar% {
  # retrieve peaks mapping to the TAD
  tad_peaks <- peak2tad_DT$peak[peak2tad_DT$region == tad]
  
  if(length(tad_peaks) == 0){
    cat("... no peak in topTAD:", tad, "\n")
    return(data.frame(
      TAD = tad,
      nbrPeaks = 0,
      nTotPeaks_cond1 = NA,
      nTotPeaks_cond2 = NA,
      ratioPeaks_cond1 = NA,
      ratioPeaks_cond2 = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  tad_binaryMat <- binaryMat[binaryMat$peak_name %in% tad_peaks,]
  nTotPeaks_cond1 <- sum(tad_binaryMat[, grepl(paste0("_", cond1), colnames(tad_binaryMat))])
  ratioPeaks_cond1 <- nTotPeaks_cond1/prod(dim(tad_binaryMat[, grepl(paste0("_", cond1), colnames(tad_binaryMat))]))
  nTotPeaks_cond2 <- sum(tad_binaryMat[, grepl(paste0("_", cond2), colnames(tad_binaryMat))])
  ratioPeaks_cond2 <- nTotPeaks_cond2/prod(dim(tad_binaryMat[, grepl(paste0("_", cond2), colnames(tad_binaryMat))]))
  
  data.frame(
    TAD = tad,
    nbrPeaks = nrow(tad_binaryMat),
    nTotPeaks_cond1 = round(nTotPeaks_cond1,4),
    nTotPeaks_cond2 = round(nTotPeaks_cond2,4),
    ratioPeaks_cond1 = round(ratioPeaks_cond1, 4),
    ratioPeaks_cond2 = round(ratioPeaks_cond2, 4),
    stringsAsFactors = FALSE
  )
}
rownames(topTADs_histPeaks_FT) <- NULL
colnames(topTADs_histPeaks_FT) <- gsub("cond1", cond1, colnames(topTADs_histPeaks_FT))
colnames(topTADs_histPeaks_FT) <- gsub("cond2", cond2, colnames(topTADs_histPeaks_FT))

outFile <- file.path(outFold, paste0(histMark, "_", cond1, "_", cond2, "_", "topTADs_histPeaks_FT.Rdata"))
save(topTADs_histPeaks_FT, file = outFile)
if(file.exists(outFile)) {
  cat(paste0("... written: ", outFile, "\n"))
} else {
  stop("error\n")
}


# High frequency events were
# extrapolated by simple summation across 
# rows from these bit tables. Permutation analysis was used to determine which sets of VELs, by 
# frequency, were concordant across samples greater than would be expected by chance. 



########################################################################################################
########################################################################################################

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
