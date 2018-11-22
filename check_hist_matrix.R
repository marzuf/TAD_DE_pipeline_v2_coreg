library(foreach)
library(doMC)
library(dplyr)
library(reshape2)

# Rscript expr_hist_peaks_topTADs.R

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


outFold <- file.path(setDir, 
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg",
                     "CHECK_HIST_MATRIX",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))


logFile <- file.path(outFold, "check_hist_matrix.txt")
system(paste0("rm -f ", logFile))

# TOP TADs
cat("... load topTADs data\n")
topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region


file1 <- file.path(setDir, "//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/EXPR_HIST_PEAKS_TOPTADS/TCGAcrc_msi_mss/allHist_topTADs_histPeaksExpr_DT.Rdata")
file2 <- file.path(setDir, "//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/HIST_PEAKS_TOPTADS/TCGAcrc_msi_mss/allHist_topTADs_histPeaks_DT.Rdata")
file3 <- file.path(setDir, "//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/HIST_PEAKS_TOPTADS/TCGAcrc_msi_mss/allHist_topTADs_histPeaks_withMean_DT.Rdata")

dt1 <- eval(parse(text = load(file1)))
dt2 <- eval(parse(text = load(file2)))
dt3 <- eval(parse(text = load(file3)))

tadBySample_DT <- dt1


nrow(tadBySample_DT)
tadBySample_DT$variable[tadBySample_DT$variable == "sum_expr_H3K27ac"] <- "sum_expr"
tadBySample_DT$variable[tadBySample_DT$variable == "sum_expr_H3K4me1"] <- "sum_expr"
tadBySample_DT <- unique(tadBySample_DT)
nrow(tadBySample_DT)

tadBySample_DT <-  reshape(tadBySample_DT, direction="wide", idvar=c("TAD", "samples", "status"), timevar=c("variable"))
colnames(tadBySample_DT) <- gsub("value.", "", colnames(tadBySample_DT))


tadBySample_DT <- tadBySample_DT[order(tadBySample_DT$TAD, tadBySample_DT$status, tadBySample_DT$samples),]

tadAgg_DT <- dt3


top_tads <- unique(as.character(topTADs_DT$region))

top_tads <- c("chr1_TAD150", "chr12_TAD81", "chr6_TAD58")


for(tad in top_tads) {
  
  cat(paste0("*** ", tad, "*** \n"), file = logFile, append=T)
  cat(paste0("*** ", tad, "*** \n"))
  
  sub_tadAgg_DT <- tadAgg_DT[tadAgg_DT$TAD == tad,]
  
  write.table(sub_tadAgg_DT, row.names = F, col.names = T, append = T, file = logFile, quote=F)
  write.table(sub_tadAgg_DT, row.names = F, col.names = T, append = T, file = "", quote=F)
  cat(paste0("\n"), file = logFile, append=T)
  
  sub_tadBySample_DT <- tadBySample_DT[tadBySample_DT$TAD == tad,]
  
  write.table(sub_tadBySample_DT, row.names = F, col.names = T, append = T, file = logFile, quote=F)
  cat(paste0("\n"), file = logFile, append=T)
  cat(paste0("\n"), file = logFile, append=T)
  
}


#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))

