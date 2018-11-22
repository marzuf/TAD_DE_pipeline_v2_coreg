library(foreach)
library(doMC)
library(dplyr)
library(reshape2)

# Rscript expr_hist_peaks_topTADs.R

# for each sample, retrieve the hist. and expression for each TAD
# too see if there are cell lines in which I have good correlation
# expression and histones and other that do not

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
                     "EXPR_HIST_PEAKS_TOPTADS",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))

# TOP TADs
cat("... load topTADs data\n")
topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region


# prepare expression data
cat("... load expression matrix\n")
exprF <- file.path(setDir, 
                   "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/FPKM",
                   curr_dataset, 
                   "gene_fpkm_with_region_DT.Rdata")
stopifnot(file.exists(exprF))
exprDT <- eval(parse(text = load(exprF)))

expr_samples <- colnames(exprDT)[grep("FPKM_", colnames(exprDT))]
expr_samples <- gsub(paste0(".+_(.+)_.+"), "\\1", expr_samples)
# expr_samples

all_hist_marks <- c("H3K27ac", "H3K4me1")

if(buildTable) {
  
  allHist_topTADs_histPeaksExpr_DT <- foreach(histMark = all_hist_marks, .combine='rbind') %do% {
    
    
    cat(paste0("> START: ", histMark, "\n"))
    
    
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
    
    countMat_samples <- colnames(countMatrix)[-c(1:4)]
    countMat_samples <- gsub(paste0(".+_(.+)_.+_.+"), "\\1", countMat_samples)
    countMat_samples
    
    countMat_status <- colnames(countMatrix)[-c(1:4)]
    countMat_status <- gsub(paste0(".+_.+_(.+)_.+"), "\\1", countMat_status)
    countMat_status
    
    samples_status <- setNames(countMat_status, countMat_samples)
    
    # select only those for which I have expression values
    # stopifnot(expr_samples %in% countMat_samples) # V9M has histones but not expression
    inter_samples <- intersect(countMat_samples, expr_samples)
    stopifnot(length(inter_samples) > 0)
    
    # for each sample: mean expression and mean peak counts
    countMat_sub <- countMatrix[,-c(1:4)]
    colnames(countMat_sub) <- gsub(paste0(".+_(.+)_.+_.+"), "\\1", colnames(countMat_sub))
    stopifnot(inter_samples %in% colnames(countMat_sub))
    countMat_sub <- countMat_sub[,inter_samples]
    countMat_sub <- cbind(countMatrix[,c(1:4)], countMat_sub)
    
    exprMat_sub <- exprDT
    exprMat_sub <- exprMat_sub[,grep("FPKM_", colnames(exprMat_sub))]
    colnames(exprMat_sub) <- gsub(paste0(".+_(.+)_.+"), "\\1", colnames(exprMat_sub))
    stopifnot(inter_samples %in% colnames(exprMat_sub))
    exprMat_sub <- exprMat_sub[, inter_samples]
    exprMat_sub <- cbind(exprDT$region, exprMat_sub)
    colnames(exprMat_sub)[1] <- "region"
    
    #### prepare binary matrix
    
    binaryMatrixF <-  file.path(setDir, 
                                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_BINARYMATRIX",
                                     curr_dataset, 
                                     histMark, 
                                     paste0(histMark, "_", "MSI_MSS", "_binaryMatrix.Rdata"))
    stopifnot(file.exists(binaryMatrixF))
    cat("... load binary matrix\n")
    binaryMat <- eval(parse(text = load(binaryMatrixF)))
    colnames(binaryMat) <- gsub(".+_(.+)_.+", "\\1", colnames(binaryMat))
    stopifnot(inter_samples %in% colnames(binaryMat))
    cat("... prepare binary matrix\n")
    binaryMat_sub <- binaryMat[,c("peak_name", "peak_chromo", "peak_start", "peak_end", inter_samples)] 
    
    # for each TAD, each cell line: correlation expression and hist marks
    
    topTADs_histPeaks_DT <- foreach(i_tad = 1:nrow(topTADs_DT), .combine='rbind') %dopar% {
      
      curr_tad <- topTADs_DT$region[i_tad]
      curr_chromo <- topTADs_DT$chromo[i_tad]
      curr_start <- topTADs_DT$start[i_tad]
      curr_end <- topTADs_DT$end[i_tad]
      
      # select peaks that are in curr_tad
      tad_countMat <- countMat_sub[countMat_sub$peak_chromo == curr_chromo &
                                     countMat_sub$peak_start >= curr_start & 
                                     countMat_sub$peak_end <= curr_end,]
      tad_countMat <- tad_countMat[, inter_samples]
      
      
      tad_binaryMat <- binaryMat_sub[binaryMat_sub$peak_chromo == curr_chromo &
                                       binaryMat_sub$peak_start >= curr_start & 
                                       binaryMat_sub$peak_end <= curr_end,]
      tad_binaryMat <- tad_binaryMat[, inter_samples]
      
      
      tad_exprMat <- exprMat_sub[exprMat_sub$region == curr_tad,]
      tad_exprMat <- tad_exprMat[, inter_samples]
      
      data.frame(
        TAD = curr_tad,
        samples = inter_samples,
        status = samples_status[inter_samples],
        sum_count = colSums(tad_countMat),
        sum_binary = colSums(tad_binaryMat),
        sum_expr = colSums(tad_exprMat),
        stringsAsFactors = FALSE
      )
      
    }
    
    rownames(topTADs_histPeaks_DT) <- NULL
    
    topTADs_histPeaks_DT_m <- melt(topTADs_histPeaks_DT, id = c("TAD", "samples", "status"))
    topTADs_histPeaks_DT_m$TAD <- as.character(topTADs_histPeaks_DT_m$TAD)
    topTADs_histPeaks_DT_m$samples <- as.character(topTADs_histPeaks_DT_m$samples)
    topTADs_histPeaks_DT_m$status <- as.character(topTADs_histPeaks_DT_m$status)
    topTADs_histPeaks_DT_m$variable <- as.character(topTADs_histPeaks_DT_m$variable)
    topTADs_histPeaks_DT_m$variable <- paste0(topTADs_histPeaks_DT_m$variable,  "_", histMark)
    topTADs_histPeaks_DT_m
  }
  
  outFile <- file.path(outFold,"allHist_topTADs_histPeaksExpr_DT.Rdata" )
  save(allHist_topTADs_histPeaksExpr_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold,"allHist_topTADs_histPeaksExpr_DT.Rdata" )
  allHist_topTADs_histPeaksExpr_DT <- eval(parse(text = load(outFile)))
}

# the expression values are duplicated
nrow(allHist_topTADs_histPeaksExpr_DT)
allHist_topTADs_histPeaksExpr_DT$variable[allHist_topTADs_histPeaksExpr_DT$variable == "sum_expr_H3K27ac"] <- "sum_expr"
allHist_topTADs_histPeaksExpr_DT$variable[allHist_topTADs_histPeaksExpr_DT$variable == "sum_expr_H3K4me1"] <- "sum_expr"
allHist_topTADs_histPeaksExpr_DT <- unique(allHist_topTADs_histPeaksExpr_DT)
nrow(allHist_topTADs_histPeaksExpr_DT)


histMark_expr_corr_DT <- do.call("rbind", 
  by(allHist_topTADs_histPeaksExpr_DT, allHist_topTADs_histPeaksExpr_DT$samples, function(x) {
  
  expr_tads <- x$TAD[x$variable == "sum_expr"]
  k27_tads <- x$TAD[x$variable == "sum_count_H3K27ac"]
  k4_tads <- x$TAD[x$variable == "sum_count_H3K4me1"]
  
  h3k27ac_counts <- setNames(x$value[x$variable == "sum_count_H3K27ac"], x$TAD[x$variable == "sum_count_H3K27ac"])
  h3k4me1_counts <- setNames(x$value[x$variable == "sum_count_H3K4me1"],  x$TAD[x$variable == "sum_count_H3K4me1"])
  expr <- setNames(x$value[x$variable == "sum_expr"],  x$TAD[x$variable == "sum_expr"])
  
  return(
    c(    sample = as.character(unique(x$samples)),
          h3k27ac_expr_corr = cor(h3k27ac_counts[intersect(expr_tads, k27_tads)], expr[intersect(expr_tads, k27_tads)]), 
          h3k4me1c_expr_corr = cor(h3k4me1_counts[intersect(expr_tads, k4_tads)], expr[intersect(expr_tads, k4_tads)])
    )
    )
})
)

rownames(histMark_expr_corr_DT) <- NULL

histMark_expr_corr_DT <- as.data.frame(histMark_expr_corr_DT)

outFile <- file.path(outFold,"histMark_expr_corr_DT.Rdata" )
save(histMark_expr_corr_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold,"histMark_expr_corr_DT.txt" )
write.table(histMark_expr_corr_DT, file = outFile, quote=F, sep="\t", row.names=F, col.names=T, append=F)
cat(paste0("... written: ", outFile, "\n"))

#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


