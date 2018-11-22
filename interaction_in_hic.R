# Rscript interaction_in_hic.R

library(reshape2)
library(methods)
library(Matrix)
# compare to v0: correct for distance !!!

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")

source("coreg_utils.R")

cat(paste0("> START ", "interaction_in_hic.R****",  "\n"))

inFold <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom/input_TopDom")
inFold <- file.path(setDir, "/mnt/nas_marie/TAD_DA_pipeline/TAD_call_pipeline_TopDom/input_TopDom")

tissue <- "SB"
binSize <- 40000
binSizeKb <- binSize/1000

nTop <- 5

outFold <- "INTERACTION_IN_HIC"
if(!SSHFS) system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0(tissue, "_HiC_TAD_interactions.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile=""

# load info about the peaks
####################################################

# HARD CODED MAIN FOLDER
mainFold <- file.path(setDir, 
                      "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017")
h3k27ac_refPeakFile <- file.path(mainFold, "H3K27ac", paste0("H3K27ac", "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(h3k27ac_refPeakFile))
h3k27ac_refPeak_DT <- read.delim(h3k27ac_refPeakFile, header=F, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"), stringsAsFactors = FALSE)
stopifnot(is.numeric(h3k27ac_refPeak_DT$peak_start))
stopifnot(is.numeric(h3k27ac_refPeak_DT$peak_end))

h3k4me1_refPeakFile <- file.path(mainFold, "H3K4me1", paste0("H3K4me1", "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(h3k4me1_refPeakFile))
h3k4me1_refPeak_DT <- read.delim(h3k4me1_refPeakFile, header=F, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"), stringsAsFactors = FALSE)
stopifnot(is.numeric(h3k4me1_refPeak_DT$peak_start))
stopifnot(is.numeric(h3k4me1_refPeak_DT$peak_end))

txt <- paste0("!!! HARD-CODED !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... binSize\t=\t", binSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... tissue\t=\t", tissue,"\n")
printAndLog(txt, logFile)
txt <- paste0("... nTop interactions\t=\t", nTop,"\n")
printAndLog(txt, logFile)

all_TADs <- list(
  c(TAD_name = "chr1_TAD150", chromo = "chr1", start = 89440001, end = 89920000),
  c(TAD_name = "chr12_TAD81", chromo = "chr12", start = 54160001, end = 54600000),
  c(TAD_name = "chr6_TAD58", chromo = "chr6", start = 32520001, end = 32840000)
)  

for(i in 1:length(all_TADs)) {
  
  curr_tad <- all_TADs[[i]][["TAD_name"]]
  curr_chromo <- all_TADs[[i]][["chromo"]]
  curr_start <- as.numeric(all_TADs[[i]][["start"]])
  curr_end <- as.numeric(all_TADs[[i]][["end"]])
  
  txt <- paste0("\n*** CONSENSUS TAD:\t", curr_tad, "\t", curr_start, " - ", curr_end, " ***\n")
  printAndLog(txt, logFile)
  
  tissueFile <- file.path(inFold, paste0(tissue, "_", curr_chromo, "_", binSizeKb, "k_matrix_pos_zero.txt"))
#  cat(paste0("tissueFile = ", tissueFile, "\n"))
  stopifnot(file.exists(tissueFile))  
  
  cat("... load Hi-C data\n")
  hic_DT <- read.delim(tissueFile, header=F, stringsAsFactors = FALSE)
  hic_DT_s <- hic_DT
  hic_DT[1:3,1:10]
  stopifnot( ncol(hic_DT) == nrow(hic_DT) + 3)
  
  hic_DT <- hic_DT[,-c(1:3)]
  stopifnot( ncol(hic_DT) == nrow(hic_DT) )
  
  # index of the start (1 -> 1, 40001 -> 2)
  startBin <- (curr_start-1)/binSize + 1
  # index of the start (40000 -> 1, 80000 -> 2)
  endBin <- curr_end/binSize

  stopifnot(endBin >= startBin)  
  
  sub_hic_DT <- hic_DT[startBin:endBin, startBin:endBin]
  cat("... subset Hi-C matrix\n")
  dim(sub_hic_DT)  
  
  # replace the diago with 0 
  diag(sub_hic_DT) <- 0
  
  stopifnot(dim(sub_hic_DT)[1] == dim(sub_hic_DT)[2])
  stopifnot(dim(sub_hic_DT)[1] > 0)
  # change here compared to v0: divide by diago mean to correct for distance
  notCorrected_sub_hic_DT <- sub_hic_DT
  colnames(notCorrected_sub_hic_DT) <- rownames(notCorrected_sub_hic_DT) <- 1:dim(sub_hic_DT)[1]
  notCorrected_sub_hic_DT <- melt(as.matrix(notCorrected_sub_hic_DT))
  colnames(notCorrected_sub_hic_DT) <- c("bin1", "bin2","value")
  notCorrected_sub_hic_DT$dist <- abs(notCorrected_sub_hic_DT$bin1-notCorrected_sub_hic_DT$bin2)
  notCorrected_sub_hic_DT <- notCorrected_sub_hic_DT[notCorrected_sub_hic_DT$dist > 0,]
  subDT <- notCorrected_sub_hic_DT[notCorrected_sub_hic_DT$dist == 1,]
  corrected_sub_hic_DT <- do.call(rbind, by(notCorrected_sub_hic_DT, notCorrected_sub_hic_DT$dist, function(subDT) {
    data.frame(bin1=subDT$bin1,
               bin2=subDT$bin2,
               value = subDT$value/mean(subDT$value),
               stringsAsFactors = FALSE)
  } ))
  correctedMat <- new("dgTMatrix", i = as.integer(corrected_sub_hic_DT$bin1-1),
                    j = as.integer(corrected_sub_hic_DT$bin2-1),
                    x = corrected_sub_hic_DT$value,
                  Dim = as.integer(dim(sub_hic_DT)))
  
  stopifnot(isSymmetric(correctedMat))
  diag(correctedMat) <- 0
  
  sub_hic_DT <- as.data.frame(as.matrix(correctedMat))
  
  j <- 0
  
  stopifnot(is.numeric(h3k27ac_refPeak_DT$peak_start))
  stopifnot(is.numeric(h3k27ac_refPeak_DT$peak_end))
  stopifnot(is.numeric(h3k4me1_refPeak_DT$peak_start))
  stopifnot(is.numeric(h3k4me1_refPeak_DT$peak_end))
  
  while(j < nTop) {
    
    ### MAX ONLY:
    max_interaction_row <- which.max(as.numeric(unlist(sub_hic_DT)))%/%ncol(sub_hic_DT) + 1
    max_interaction_col <- which.max(as.numeric(unlist(sub_hic_DT)))%%ncol(sub_hic_DT)
    if(max_interaction_col == 0) {
      max_interaction_col <- ncol(sub_hic_DT)
      max_interaction_row <- max_interaction_row - 1
    }
    stopifnot( sub_hic_DT[max_interaction_row, max_interaction_col] == max(sub_hic_DT) )
    stopifnot( sub_hic_DT[max_interaction_col, max_interaction_row] == max(sub_hic_DT) )
    
    binA_start <- curr_start + (max_interaction_row-1)*binSize
    binA_end <- binA_start + binSize - 1 
    binA_bin <- (binA_start-1)/binSize + 1
    
    binB_start <- curr_start + (max_interaction_col-1)*binSize
    binB_end <- binB_start + binSize - 1
    binB_bin <- (binB_start-1)/binSize + 1
    
    
    txt <- paste0("... ", curr_tad, " - interaction top ", j+1, "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... ", "interaction value = \t", max(sub_hic_DT), "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... ", "interaction binA = \t", binA_start , " - ", binA_end , "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... ", "interaction binB = \t", binB_start , " - ", binB_end , "\n")
    printAndLog(txt, logFile)
    
    # -> will not be true in this version because is hic_DT is not normalized for distance
    # stopifnot(hic_DT[binA_bin,binB_bin] ==  max(sub_hic_DT) )
    
    sub_hic_DT[max_interaction_row, max_interaction_col] <- 0
    sub_hic_DT[max_interaction_col, max_interaction_row] <- 0
    
    
    ### CHECK OVERLAP WITH PEAKS
    # H3K27ac peaks that fall within the bins
    # binA
    txt <- paste0("H3K27ac peaks in binA\n")
    printAndLog(txt, logFile)
    write.table(h3k27ac_refPeak_DT[
      h3k27ac_refPeak_DT$peak_chromo == curr_chromo &
      h3k27ac_refPeak_DT$peak_start >= binA_start &
      h3k27ac_refPeak_DT$peak_end <= binA_end,], file = logFile, col.names=T, row.names=F, sep="\t", quote=F, append=T)
    # binB
    txt <- paste0("H3K27ac peaks in binB\n")
    printAndLog(txt, logFile)
    write.table(h3k27ac_refPeak_DT[
      h3k27ac_refPeak_DT$peak_chromo == curr_chromo &
        h3k27ac_refPeak_DT$peak_start >= binB_start &
        h3k27ac_refPeak_DT$peak_end <= binB_end,], file = logFile, col.names=T, row.names=F, sep="\t", quote=F, append=T)
    
    # H3K4me1 peaks that fall within the bins
    # binA
    txt <- paste0("H3K4me1 peaks in binA\n")
    printAndLog(txt, logFile)
    write.table(h3k4me1_refPeak_DT[
      h3k4me1_refPeak_DT$peak_chromo == curr_chromo &
        h3k4me1_refPeak_DT$peak_start >= binA_start &
        h3k4me1_refPeak_DT$peak_end <= binA_end,], file = logFile, col.names=T, row.names=F, sep="\t", quote=F, append=T)
    # binB
    txt <- paste0("H3K4me1 peaks in binA\n")
    printAndLog(txt, logFile)
    write.table(h3k4me1_refPeak_DT[
      h3k4me1_refPeak_DT$peak_chromo == curr_chromo &
        h3k4me1_refPeak_DT$peak_start >= binB_start &
        h3k4me1_refPeak_DT$peak_end <= binB_end,], file = logFile, col.names=T, row.names=F, sep="\t", quote=F, append=T)
    
    
    j <- j+1
    
  }
  
  
}
