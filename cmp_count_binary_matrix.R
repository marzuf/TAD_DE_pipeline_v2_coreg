# /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss/H3K27ac/H3K27ac_MSI_MSS_binaryMatrix.Rdata

# /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss/H3K27ac/H3K27ac_MSI_countMatrix.Rdata

SSHFS <- TRUE
setDir <- ifelse(SSHFS, "/media/electron", "")

histMark <- "H3K27ac"
cond1 <- "MSI"
cond2 <- "MSS"

expType <- "mean"

binaryFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_BINARYMATRIX/TCGAcrc_msi_mss",
                                histMark, paste0(histMark, "_",cond1, "_", cond2, "_binaryMatrix.Rdata"))
stopifnot(file.exists(binaryFile))
binaryMat <- eval(parse(text = load(binaryFile)))
head(binaryMat,2)

countFile1 <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss",
                        histMark, paste0(histMark, "_",cond1, "_countMatrix.Rdata"))
stopifnot(file.exists(countFile1))
countMat_cond1 <- eval(parse(text = load(countFile1)))
head(countMat_cond1, 2)

countFile2 <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/HIST_COUNTMATRIX/TCGAcrc_msi_mss",
                        histMark, paste0(histMark, "_",cond2, "_countMatrix.Rdata"))
stopifnot(file.exists(countFile2))
countMat_cond2 <- eval(parse(text = load(countFile2)))
head(countMat_cond2, 2)


stopifnot(all.equal(countMat_cond1[,c("peak_name", "peak_chromo", "peak_start", "peak_end")], countMat_cond2[,c("peak_name", "peak_chromo", "peak_start", "peak_end")]))

peakRef_countDT <- countMat_cond1[,c("peak_name", "peak_chromo", "peak_start", "peak_end")]
countMat_cond1_expType_DT <- countMat_cond1[, grepl(paste0("_", cond1, "_", expType, "$"), colnames(countMat_cond1))]
countMat_cond2_expType_DT <- countMat_cond2[, grepl(paste0("_", cond2, "_", expType, "$"), colnames(countMat_cond2))]

countMat <- cbind(cbind(peakRef_countDT, countMat_cond1_expType_DT), countMat_cond2_expType_DT)

# binarize the countMat 
countMat_cond1_expType_DT_bin <- countMat_cond1_expType_DT
countMat_cond1_expType_DT_bin[countMat_cond1_expType_DT_bin > 0 ]<- 1

countMat_cond2_expType_DT_bin <- countMat_cond2_expType_DT
countMat_cond2_expType_DT_bin[countMat_cond2_expType_DT_bin > 0 ]<- 1

# split binary matrix for each condition
binaryMat_cond1 <- binaryMat[,grepl(paste0("_", cond1), colnames(binaryMat))]
binaryMat_cond2 <- binaryMat[,grepl(paste0("_", cond2), colnames(binaryMat))]

stopifnot(all.equal(binaryMat[,c("peak_name", "peak_chromo", "peak_start", "peak_end")], countMat_cond2[,c("peak_name", "peak_chromo", "peak_start", "peak_end")]))
stopifnot(all.equal(binaryMat[,c("peak_name", "peak_chromo", "peak_start", "peak_end")], countMat_cond1[,c("peak_name", "peak_chromo", "peak_start", "peak_end")]))

