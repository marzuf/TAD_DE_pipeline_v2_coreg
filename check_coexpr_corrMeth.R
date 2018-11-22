# question: how different would it be to take Spearman on raw rnaseq vs. PCC on qqnorm ?
# (just look at one dataset)


startTime <- Sys.time()
cat(paste0("> Rscript check_coexpr_corrMeth.R\n"))

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

## TOP 3
# Rscript check_coexpr_corrMeth.R TCGAcrc_msi_mss # running E
# Rscript check_coexpr_corrMeth.R GSE74927_neg_pos # running E
# Rscript check_coexpr_corrMeth.R GSE102073_stic_nostic # running E
# # LAST 3
# Rscript check_coexpr_corrMeth.R GSE65540_before_after # running P 
# Rscript check_coexpr_corrMeth.R GSE84231_lhb_rhb # running P
# Rscript check_coexpr_corrMeth.R GSE86356_tibMD1_tibNorm # running P

caller ="TopDom"
args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- args[1]
stopifnot(length(args) == 1)

corMethod_qqnorm <- "pearson"
corMethod_raw <- "spearman"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 400, 7)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CHECK_COEXPR_CORRMETH", paste0(curr_dataset, "_qqnorm_", corMethod_qqnorm, "_vs_raw_", corMethod_raw))
system(paste0("mkdir -p ", outFold))

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

script0_name <- "0_prepGeneData"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))

qqnormDT <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
stopifnot(names(pipeline_geneList) %in% rownames(qqnormDT))
qqnormDT <- qqnormDT[rownames(qqnormDT) %in% names(pipeline_geneList),]
stopifnot(nrow(qqnormDT) == length(pipeline_geneList))
rownames(qqnormDT) <- pipeline_geneList[rownames(qqnormDT)]
stopifnot(setequal(pipeline_geneList, rownames(qqnormDT)))

cor_qqnormMat <- cor(t(qqnormDT), method = corMethod_qqnorm)
stopifnot(nrow(cor_qqnormMat) == nrow(qqnormDT))
stopifnot(ncol(cor_qqnormMat) == nrow(qqnormDT))

# outFile <- file.path(outFold, paste0("cor_qqnormMat_", corMethod_qqnorm, ".Rdata"))
# save(cor_qqnormMat, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))
# cor_qqnormMat <- eval(parse(text = load(outFile)))
lowerTri_cor_qqnormMat <- cor_qqnormMat[lower.tri(cor_qqnormMat, diag=F)]

exprDT <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "rna_rnaseqDT.Rdata"))))
stopifnot(names(pipeline_geneList) %in% rownames(exprDT))
exprDT <- exprDT[rownames(exprDT) %in% names(pipeline_geneList),]
stopifnot(nrow(exprDT) == length(pipeline_geneList))
rownames(exprDT) <- pipeline_geneList[rownames(exprDT)]
stopifnot(setequal(pipeline_geneList, rownames(exprDT)))

cor_exprMat <- cor(t(exprDT), method = corMethod_raw)
stopifnot(nrow(cor_exprMat) == nrow(exprDT))
stopifnot(ncol(cor_exprMat) == nrow(exprDT))

# outFile <- file.path(outFold, paste0("cor_exprMat_", corMethod_raw, ".Rdata"))
# save(cor_exprMat, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))
# cor_exprMat <- eval(parse(text = load(outFile)))
lowerTri_cor_exprMat <- cor_exprMat[lower.tri(cor_exprMat, diag=F)]

stopifnot(rownames(cor_qqnormMat) == rownames(cor_exprMat))
stopifnot(colnames(cor_qqnormMat) == colnames(cor_exprMat))

outFile <- file.path(outFold, paste0("corrValues_qqnorm", corMethod_qqnorm, "_vs_raw", corMethod_raw, ".", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = lowerTri_cor_exprMat, # raw
     y = lowerTri_cor_qqnormMat, # qqnorm
     main = paste0(curr_dataset, " - corr. values qqnorm ", corMethod_qqnorm, " vs. raw ", corMethod_raw),
     xlab = paste0("raw ", corMethod_raw, " corr. values"),
     ylab = paste0("qqnorm ", corMethod_qqnorm, " corr. values"),
     pch=16, cex = 0.7)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

