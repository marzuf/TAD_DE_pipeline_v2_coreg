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
cond1 <- "MSI"
cond2 <- "MSS"
histMark <- "H3K27ac"

buildTable <- TRUE
cat(paste0("... buildTable = ", as.character(buildTable), "\n"))

outFold <- file.path(setDir, 
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg",
                     "HIST_PEAKS_LISTRATIO_v0",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "histpeaks_list_ratio_logFile.txt")
system(paste0("rm -f ", logFile))

# TOP TADs
cat("... load topTADs data\n")
topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region

topTADs <- topTADs[topTADs %in% c("chr1_TAD150", "chr6_TAD58", "chr12_TAD81")]

stopifnot(length(topTADs) > 0)

# HARD-CODED: FOCUS ON THE FOLLOWING TADs

txt <- paste0("!!! HARD-CODED !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... focus on following TADs:\t", paste0(topTADs, collapse = ", "), "\n")
printAndLog(txt, logFile)
txt <- paste0("... work with peak signal values aggregated by \"mean\"", "\n")
printAndLog(txt, logFile)

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
    binaryMatrix$peak_start <- as.numeric(as.character(binaryMatrix$peak_start))
    binaryMatrix$peak_end <- as.numeric(as.character(binaryMatrix$peak_end))
    
    binaryMatrix <- binaryMatrix[order(binaryMatrix$peak_chromo, binaryMatrix$peak_start, binaryMatrix$peak_end),]
    
    
    
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
    countMatrix$peak_start <- as.numeric(as.character(countMatrix$peak_start))
    countMatrix$peak_end <- as.numeric(as.character(countMatrix$peak_end))
    
    countMatrix <- countMatrix[order(countMatrix$peak_chromo, countMatrix$peak_start, countMatrix$peak_end),]
    
    
    
    topTADs_histPeaks_DT <- foreach(curr_tad = topTADs, .combine='rbind') %do% {
        
      txt <- paste0("> ", histMark, " - *** ", curr_tad, "***\n")
      printAndLog(txt, logFile)
      
      tad_idx <- which(topTADs_DT$region == curr_tad)
      stopifnot(length(tad_idx) == 1)
      
      curr_chromo <- topTADs_DT$chromo[tad_idx]
      curr_start <- topTADs_DT$start[tad_idx]
      curr_end <- topTADs_DT$end[tad_idx]

      stopifnot(is.numeric(curr_start))
      stopifnot(is.numeric(curr_end))
      stopifnot(is.numeric(countMatrix$peak_start))
      stopifnot(is.numeric(countMatrix$peak_end))
      stopifnot(is.numeric(binaryMatrix$peak_start))
      stopifnot(is.numeric(binaryMatrix$peak_end))
      
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
      

      
      ### >>> PREPARE BINARY FINAL TABLE
      tmpBin_DT <- tad_binMat
      tmpBin_DT <- melt(tmpBin_DT, 
                        id = c("peak_chromo", "peak_start", "peak_end", "peak_name"))
      tmpBin_DT$sample <- gsub(".+_(.+)_.+", "\\1", tmpBin_DT$variable)
      tmpBin_DT$status <- gsub(".+_.+_(.+)", "\\1", tmpBin_DT$variable)
      tmpBin_DT$variable <- NULL
      stopifnot(length(setdiff(tmpBin_DT$status, c(cond1,cond2))) == 0)
      
      binRatio_DT <- do.call(rbind, by(tmpBin_DT, tmpBin_DT$peak_name, function(subDT) {
        nSamplesCond1 <-  sum(subDT$status == cond1)
        nSamplesCond2 <-  sum(subDT$status == cond2)
        stopifnot(nSamplesCond1 + nSamplesCond2 == nrow(subDT))
        nPeaksCond1 <- sum(subDT$value[subDT$status == cond1])
        nPeaksCond2 <- sum(subDT$value[subDT$status == cond2])
        stopifnot(nPeaksCond1 + nPeaksCond2 == sum(subDT$value) )
        c(
          nSamplesCond1 = nSamplesCond1,
          nSamplesCond2 = nSamplesCond2,
          nPeaksCond1 = nPeaksCond1,
          nPeaksCond2 = nPeaksCond2
          )
      }))
      
      peakposDT <- tmpBin_DT[,c("peak_name", "peak_chromo", "peak_start", "peak_end")]
      peakposDT <- unique(peakposDT)
      rownames(peakposDT) <- peakposDT$peak_name
      peakposDT$peak_name <- NULL
      
      stopifnot(length(setdiff(rownames(peakposDT), rownames(binRatio_DT))) == 0)
      final_binRatio_DT <- merge(peakposDT, binRatio_DT, by="row.names")
      colnames(final_binRatio_DT)[colnames(final_binRatio_DT) == "Row.names"] <- "peak_name"

      final_binRatio_DT$ratioBinCond1 <- final_binRatio_DT$nPeaksCond1/final_binRatio_DT$nSamplesCond1
      final_binRatio_DT$ratioBinCond2 <- final_binRatio_DT$nPeaksCond2/final_binRatio_DT$nSamplesCond2
      
      
      ### >>> PREPARE COUNT FINAL TABLE
      tmpCount_DT <- tad_countMat
      tmpCount_DT <- melt(tmpCount_DT, 
                        id = c("peak_chromo", "peak_start", "peak_end", "peak_name"))
      tmpCount_DT$sample <- gsub(".+_(.+)_.+_.+", "\\1", tmpCount_DT$variable)
      tmpCount_DT$status <- gsub(".+_.+_(.+)_.+", "\\1", tmpCount_DT$variable)
      tmpCount_DT$variable <- NULL
      stopifnot(length(setdiff(tmpCount_DT$status, c(cond1,cond2))) == 0)
      
      
      countRatio_DT <- do.call(rbind, by(tmpCount_DT, tmpCount_DT$peak_name, function(subDT) {
        nSamplesCond1 <-  sum(subDT$status == cond1)
        nSamplesCond2 <-  sum(subDT$status == cond2)
        stopifnot(nSamplesCond1 + nSamplesCond2 == nrow(subDT))
        sumMeanPeaksCond1 <- sum(subDT$value[subDT$status == cond1])
        sumMeanPeaksCond2 <- sum(subDT$value[subDT$status == cond2])
        stopifnot( round(sumMeanPeaksCond1 + sumMeanPeaksCond2,6) == round(sum(subDT$value),6) )
        c(
          nSamplesCond1 = nSamplesCond1,
          nSamplesCond2 = nSamplesCond2,
          sumMeanPeaksCond1 = sumMeanPeaksCond1,
          sumMeanPeaksCond2 = sumMeanPeaksCond2
        )
      }))
      
      peakposDT <- tmpCount_DT[,c("peak_name", "peak_chromo", "peak_start", "peak_end")]
      peakposDT <- unique(peakposDT)
      rownames(peakposDT) <- peakposDT$peak_name
      peakposDT$peak_name <- NULL
      
      stopifnot(length(setdiff(rownames(peakposDT), rownames(countRatio_DT))) == 0)
      final_countRatio_DT <- merge(peakposDT, countRatio_DT, by="row.names")
      colnames(final_countRatio_DT)[colnames(final_countRatio_DT) == "Row.names"] <- "peak_name"

      final_countRatio_DT$meanCountCond1 <- final_countRatio_DT$sumMeanPeaksCond1/final_countRatio_DT$nSamplesCond1
      final_countRatio_DT$meanCountCond2 <- final_countRatio_DT$sumMeanPeaksCond2/final_countRatio_DT$nSamplesCond2
      
      
      
      all_DT <- merge(final_binRatio_DT, final_countRatio_DT,by=c("peak_name", "peak_chromo", "peak_start", "peak_end", "nSamplesCond1", "nSamplesCond2"))
      
      stopifnot(length(all_DT$peak_name) == length(unique(all_DT$peak_name)))
      
      all_DT$TAD <- curr_tad
      all_DT <- all_DT[,c("TAD", colnames(all_DT)[colnames(all_DT) != "TAD"])]
      
      colnames(all_DT) <- gsub("Cond1", paste0("_", cond1), colnames(all_DT))
      colnames(all_DT) <- gsub("Cond2", paste0("_", cond2), colnames(all_DT))
      
      head(all_DT,2)
      
      txt <- paste0("\n")
      printAndLog(txt, logFile)
      
      all_DT <- all_DT[order(all_DT$peak_chromo, all_DT$peak_start, all_DT$peak_end),]
      
      write_DT <- all_DT[,c("peak_name", "peak_chromo", "peak_start", "peak_end", "ratioBin_MSI", "ratioBin_MSS", "meanCount_MSI", "meanCount_MSS")]
      write_DT$ratioBin_MSI <- round(write_DT$ratioBin_MSI, 4)
      write_DT$ratioBin_MSS <- round(write_DT$ratioBin_MSS, 4)
      write_DT$meanCount_MSI <- round(write_DT$meanCount_MSI, 4)
      write_DT$meanCount_MSS <- round(write_DT$meanCount_MSS, 4)
      
      write_DT1 <- write_DT
      write_DT1 <- write_DT1[order(write_DT1$ratioBin_MSI, write_DT1$meanCount_MSI, decreasing=TRUE),]
      
      logFileDT <- file.path(outFold, paste0("histpeaks_list_ratio_logFile_", histMark, "_", curr_tad , ".txt"))
      write.table(write_DT1, file = logFileDT, col.names = TRUE, row.names=F, append=F, sep="\t", quote=F)
      cat(paste0("... written: ", logFileDT, "\n"))
      
      
      write.table(write_DT, file = logFile, col.names = TRUE, row.names=F, append=T, sep="\t", quote=F)
      
      txt <- paste0("\n")
      printAndLog(txt, logFile)
      
      all_DT
    } # end-iterate over the tads
    
    topTADs_histPeaks_DT$histMark <- histMark
    topTADs_histPeaks_DT <- topTADs_histPeaks_DT[,c("histMark", colnames(topTADs_histPeaks_DT)[colnames(topTADs_histPeaks_DT) != "histMark"])]
    

  } # end-iterate over histone marks
  outFile <- file.path(outFold, "allHist_topTADs_histPeaksRatio_DT.Rdata" )
  save(allHist_topTADs_histPeaks_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "allHist_topTADs_histPeaksRatio_DT.Rdata" )
  stopifnot(file.exists(outFile))
  allHist_topTADs_histPeaks_DT <- eval(parse(text = load(outFile)))
}





#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))

