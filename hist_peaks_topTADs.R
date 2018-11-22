library(foreach)
library(doMC)
library(dplyr)
library(reshape2)


# Rscript hist_peaks_topTADs.R

# in this version assign mRNA to TAD based on entrezID

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 20))

if(SSHFS) setwd(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"))

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"
histMark <- "H3K27ac"
buildTable <- TRUE

cat(paste0("... buildTable = ", as.character(buildTable), "\n"))

outFold <- file.path(setDir, 
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg",
                     "HIST_PEAKS_TOPTADS",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "hist_peaks_logFile.txt")
system(paste0("rm -f ", logFile))

# TOP TADs
cat("... load topTADs data\n")
topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region


all_hist_marks <- c("H3K27ac", "H3K4me1")

if(buildTable) {
  allHist_topTADs_histPeaks_DT <- foreach(histMark = all_hist_marks, .combine='rbind') %do% {
    
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
    
    cond1_samples <- gsub(paste0("_", cond1), "", colnames(binaryMatrix)[grepl(paste0("_", cond1), colnames(binaryMatrix))])
    cond2_samples <- gsub(paste0("_", cond2), "", colnames(binaryMatrix)[grepl(paste0("_", cond2), colnames(binaryMatrix))])
    
    binaryMatrix$peak_name <- as.character(binaryMatrix$peak_name)
    binaryMatrix$peak_chromo <- as.character(binaryMatrix$peak_chromo)
    binaryMatrix$peak_start <- as.character(binaryMatrix$peak_start)
    binaryMatrix$peak_end <- as.character(binaryMatrix$peak_end)
    
    # HARD CODED PEAK COUNT MATRIX
    cat("... load count matrix - cond1\n")
    cond1_countMatrixF <-  file.path(setDir, 
                                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX",
                                     curr_dataset, 
                                     histMark, 
                                     paste0(histMark, "_", cond1, "_countMatrix.Rdata"))
    stopifnot(file.exists(cond1_countMatrixF))
    cond1_countMatrix <- eval(parse(text = load(cond1_countMatrixF)))
    head(cond1_countMatrix)
    
    cat("... load count matrix - cond2\n")
    cond2_countMatrixF <-  file.path(setDir, 
                                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX",
                                     curr_dataset, 
                                     histMark, 
                                     paste0(histMark, "_", cond2, "_countMatrix.Rdata"))
    stopifnot(file.exists(cond2_countMatrixF))
    cond2_countMatrix <- eval(parse(text = load(cond2_countMatrixF)))
    head(cond2_countMatrix)
    
    cat("... prepare count matrix\n")
    countMatrix <- merge(cond1_countMatrix, cond2_countMatrix, by=c("peak_name", "peak_chromo", "peak_start", "peak_end"))
    countMatrix <- countMatrix[,-grep("_mean0$", colnames(countMatrix))]
    countMatrix <- countMatrix[,-grep("_sum$", colnames(countMatrix))]
    head(countMatrix)
    
    countMatrix$peak_name <- as.character(countMatrix$peak_name)
    countMatrix$peak_chromo <- as.character(countMatrix$peak_chromo)
    countMatrix$peak_start <- as.character(countMatrix$peak_start)
    countMatrix$peak_end <- as.character(countMatrix$peak_end)
    
    # check that binary correspond to counts ?
    tmp_count <- cond1_countMatrix[,c("peak_name", "peak_chromo", "peak_start", "peak_end", paste0(cond1_samples[1], "_MSI_mean"))]
    tmp_bin <- binaryMatrix[,c("peak_chromo", "peak_start", "peak_end", "peak_name", paste0(cond1_samples[1], "_MSI"))]
    tmp_all <- merge(tmp_count, tmp_bin, by=c("peak_name", "peak_chromo", "peak_start", "peak_end"))
    boxplot(as.formula(paste0(cond1_samples[1], "_MSI_mean ~ ", cond1_samples[1], "_MSI")), data=tmp_all)
    
    tmpbin <- binaryMatrix[,-c(1:4)]
    tmpcount <- countMatrix[,-c(1:4)]
    colnames(tmpcount) <- gsub("_mean", "", colnames(tmpcount))
    stopifnot(colnames(tmpbin) == colnames(tmpcount))
    
    tmpbincount <- data.frame(bin_vect=as.numeric(as.matrix(tmpbin)), 
                              count_vect = as.numeric(as.matrix(tmpcount)))
    nrow(tmpbincount)
    boxplot(count_vect ~ bin_vect, data = tmpbincount)
    
    binaryMatrix <- binaryMatrix[order(binaryMatrix$peak_chromo, binaryMatrix$peak_start, binaryMatrix$peak_end),]
    countMatrix <- countMatrix[order(countMatrix$peak_chromo, countMatrix$peak_start, countMatrix$peak_end),]
    
    topTADs_histPeaks_DT <- foreach(i_tad = 1:nrow(topTADs_DT), .combine='rbind') %do% {
      
      
      
      
      curr_tad <- topTADs_DT$region[i_tad]
      
      cat("*** ", curr_tad, "***\n", file = logFile, append = TRUE)
      
      curr_chromo <- topTADs_DT$chromo[i_tad]
      curr_start <- topTADs_DT$start[i_tad]
      curr_end <- topTADs_DT$end[i_tad]
      # select peaks that are in curr_tad
      tad_countMat <- countMatrix[countMatrix$peak_chromo == curr_chromo &
                                    countMatrix$peak_start >= curr_start & 
                                    countMatrix$peak_end <= curr_end,]
      tad_binMat <- binaryMatrix[binaryMatrix$peak_chromo == curr_chromo &
                                   binaryMatrix$peak_start >= curr_start & 
                                   binaryMatrix$peak_end <= curr_end,]
      stopifnot(nrow(tad_countMat) == nrow(tad_binMat))
      stopifnot(tad_countMat$peak_name == tad_binMat$peak_name)
      stopifnot(ncol(tad_countMat) == (4 + length(cond1_samples) + length(cond2_samples)))
      stopifnot(ncol(tad_binMat) == (4 + length(cond1_samples) + length(cond2_samples)))
      
      tmpBin_DT <- tad_binMat
      colnames(tmpBin_DT) <- gsub(".+_(.+)_.+", "\\1", colnames(tmpBin_DT))
      tmpBin_DT <- tmpBin_DT[order(tmpBin_DT$peak_start, tmpBin_DT$peak_end),]
      tmpBin_DT$peak_name <- NULL
      
      write.table(tmpBin_DT, append=T, file = logFile, col.names=T, row.names=F, sep="\t", quote=F)
      cat("\n", file = logFile, append=T)
      
      tmpBin_DT <- tad_binMat
      tmpBin_DT <- melt(tmpBin_DT, 
                         id = c("peak_chromo", "peak_start", "peak_end", "peak_name"))
      tmpBin_DT$sample <- gsub(".+_(.+)_.+", "\\1", tmpBin_DT$variable)
      tmpBin_DT$status <- gsub(".+_.+_(.+)", "\\1", tmpBin_DT$variable)
      tmpBin_DT$variable <- NULL
      
      tmpBin_DT_agg <- aggregate(value ~ peak_chromo+peak_start+peak_end+peak_name+status, data = tmpBin_DT, FUN=sum)
      tmpBin_DT_agg <- tmpBin_DT_agg[order(tmpBin_DT_agg$peak_start, tmpBin_DT_agg$peak_end, tmpBin_DT_agg$status),]
      colnames(tmpBin_DT_agg)[colnames(tmpBin_DT_agg) == "value"] <- "sum_peaks" 
      write.table(tmpBin_DT_agg, append=T, file = logFile, col.names=T, row.names=F, sep="\t", quote=F)
      
      tmpBin_DT_aggLen <- aggregate(value ~ peak_chromo+peak_start+peak_end+peak_name+status, data = tmpBin_DT, FUN=length)
      tmpBin_DT_aggLen <- tmpBin_DT_aggLen[order(tmpBin_DT_aggLen$peak_start, tmpBin_DT_aggLen$peak_end, tmpBin_DT_aggLen$status),]
      colnames(tmpBin_DT_aggLen)[colnames(tmpBin_DT_aggLen) == "value"] <- "nSamples"
      
      tmpBin_DT_aggAll <- merge(tmpBin_DT_agg, tmpBin_DT_aggLen, by=c("peak_name", "peak_chromo", "peak_start", "peak_end", "status"))
      tmpBin_DT_aggAll$peakRatio <- tmpBin_DT_aggAll$sum_peaks/tmpBin_DT_aggAll$nSamples
      stopifnot(tmpBin_DT_aggAll$peakRatio >= 0 & tmpBin_DT_aggAll$peakRatio <= 1)
      
      write.table(tmpBin_DT_aggAll, append=T, file = logFile, col.names=T, row.names=F, sep="\t", quote=F)
      cat("\n", file = logFile, append=T)
      
      tmpBin_DT_aggAll$variable <- paste0("peak_ratio_", tmpBin_DT_aggAll$status)
      tmpBin_DT_aggAll <- tmpBin_DT_aggAll[, c("peak_name", "peak_chromo", "peak_start", "peak_end", "variable", "peakRatio")]
      # write.table(tmpBin_DT_agg, append=T, file = logFile, col.names=T, row.names=F, sep="\t", quote=F)
      tmpBin_DT_aggAll <- reshape(tmpBin_DT_aggAll, 
                                  idvar = c("peak_name", "peak_chromo", "peak_start", "peak_end"), 
                                  timevar = "variable",
                                  # varying = 5, 
                                  direction = "wide")
      colnames(tmpBin_DT_aggAll) <- gsub("^peakRatio\\.", "", colnames(tmpBin_DT_aggAll))
      write.table(tmpBin_DT_aggAll, append=T, file = logFile, col.names=T, row.names=F, sep="\t", quote=F)
      cat("\n", file = logFile, append=T)
      
      cat(paste0("... ", curr_tad, " - n peaks: ", nrow(tmpBin_DT_aggAll), "\n"), file = logFile, append=T)
      cat(paste0("... ", curr_tad, " - sum peak_ratio MSI: ", round(sum(tmpBin_DT_aggAll$peak_ratio_MSI, na.rm=T), 4), "\n"), file = logFile, append=T)
      cat(paste0("... ", curr_tad, " - sum peak_ratio MSS: ", round(sum(tmpBin_DT_aggAll$peak_ratio_MSS, na.rm=T), 4), "\n"), file = logFile, append=T)
      
      cat("\n", file = logFile, append=T)
      cat("\n", file = logFile, append=T)
      
      # NUMBER OF PEAKS
      data.frame(
        TAD = curr_tad,
        chromo = curr_chromo,
        start = curr_start,
        end = curr_end,
        histMark = histMark,
        nSamp_cond1 = length(cond1_samples),
        nSamp_cond2 = length(cond2_samples),
        binMat_totPeaks_cond1 = sum(tad_binMat[,paste0(cond1_samples, "_", cond1) ]),
        binMat_totPeaks_cond2 = sum(tad_binMat[,paste0(cond2_samples, "_", cond2) ]),
        countMat_sumPeaks_cond1 = sum(tad_countMat[,paste0(cond1_samples, "_", cond1, "_mean") ]),
        countMat_sumPeaks_cond2 = sum(tad_countMat[,paste0(cond2_samples, "_", cond2, "_mean") ]),
        stringsAsFactors = FALSE
      )
    }
    
    colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) == "nSamp_cond1"] <- paste0("nSamp_", cond1)
    colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) == "nSamp_cond2"] <- paste0("nSamp_", cond2)
    colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) == "binMat_totPeaks_cond1"] <- paste0("binMat_totPeaks_", cond1)
    colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) == "binMat_totPeaks_cond2"] <- paste0("binMat_totPeaks_", cond2)
    colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) == "countMat_sumPeaks_cond1"] <- paste0("countMat_sumPeaks_", cond1)
    colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) == "countMat_sumPeaks_cond2"] <- paste0("countMat_sumPeaks_", cond2)
    topTADs_histPeaks_DT
  }
  outFile <- file.path(outFold, "allHist_topTADs_histPeaks_DT.Rdata" )
  save(allHist_topTADs_histPeaks_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "allHist_topTADs_histPeaks_DT.Rdata" )
  stopifnot(file.exists(outFile))
  allHist_topTADs_histPeaks_DT <- eval(parse(text = load(outFile)))
}


#**************************************************
head(allHist_topTADs_histPeaks_DT)

allHist_topTADs_histPeaks_DT[, paste0("meanCount_", cond1)] <- allHist_topTADs_histPeaks_DT[, paste0("countMat_sumPeaks_", cond1)]/allHist_topTADs_histPeaks_DT[, paste0("nSamp_", cond1)]
allHist_topTADs_histPeaks_DT[, paste0("meanCount_", cond2)] <- allHist_topTADs_histPeaks_DT[, paste0("countMat_sumPeaks_", cond2)]/allHist_topTADs_histPeaks_DT[, paste0("nSamp_", cond2)]

allHist_topTADs_histPeaks_DT[, paste0("meanBin_", cond1)] <- allHist_topTADs_histPeaks_DT[, paste0("binMat_totPeaks_", cond1)]/allHist_topTADs_histPeaks_DT[, paste0("nSamp_", cond1)]
allHist_topTADs_histPeaks_DT[, paste0("meanBin_", cond2)] <- allHist_topTADs_histPeaks_DT[, paste0("binMat_totPeaks_", cond2)]/allHist_topTADs_histPeaks_DT[, paste0("nSamp_", cond2)]

allHist_topTADs_histPeaks_DT[, paste0("dir_meanCount")] <- ifelse( 
                            allHist_topTADs_histPeaks_DT[, paste0("meanCount_", cond1)] > allHist_topTADs_histPeaks_DT[, paste0("meanCount_", cond2)], cond1,
                            ifelse(allHist_topTADs_histPeaks_DT[, paste0("meanCount_", cond2)] > allHist_topTADs_histPeaks_DT[, paste0("meanCount_", cond1)], cond2, "none"))

allHist_topTADs_histPeaks_DT[, paste0("dir_meanBin")] <- ifelse(
   allHist_topTADs_histPeaks_DT[, paste0("meanBin_", cond1)] > allHist_topTADs_histPeaks_DT[, paste0("meanBin_", cond2)], cond1,
   ifelse(allHist_topTADs_histPeaks_DT[, paste0("meanBin_", cond2)] > allHist_topTADs_histPeaks_DT[, paste0("meanBin_", cond1)], cond2, "none"))

allHist_topTADs_histPeaks_DT[, paste0("dir_meanCount")] == allHist_topTADs_histPeaks_DT[, paste0("dir_meanBin")]

allHist_topTADs_histPeaks_withMean_DT <- allHist_topTADs_histPeaks_DT

outFile <- file.path(outFold, "allHist_topTADs_histPeaks_withMean_DT.Rdata" )
save(allHist_topTADs_histPeaks_withMean_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))


