# question: how different would it be to take Spearman on raw rnaseq vs. PCC on qqnorm ?
# (just look at one dataset)


startTime <- Sys.time()
cat(paste0("> Rscript check_correlationByCondition.R\n"))

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

## TOP 3
# Rscript check_correlationByCondition.R TCGAcrc_msi_mss # done
# Rscript check_correlationByCondition.R GSE74927_neg_pos # running E
# Rscript check_correlationByCondition.R GSE102073_stic_nostic # done
# # LAST 3
# Rscript check_correlationByCondition.R GSE65540_before_after # running E
# Rscript check_correlationByCondition.R GSE84231_lhb_rhb # running P
# Rscript check_correlationByCondition.R GSE86356_tibMD1_tibNorm # running P


caller ="TopDom"
curr_dataset <- "TCGAcrc_msi_mss"
SSFHS <- T

args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- args[1]
stopifnot(length(args) == 1)

myHeight <- 400
myWidth <- 400

corMethod_qqnorm <- "pearson"

buildTable <- TRUE

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CHECK_CORRELATIONBYCONDITION", paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))

settingFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline"), "SETTING_FILES_cleanInput", paste0("run_settings_", curr_dataset, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

script0_name <- "0_prepGeneData"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))

if(buildTable){
  qqnormDT <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
  stopifnot(names(pipeline_geneList) %in% rownames(qqnormDT))
  qqnormDT <- qqnormDT[rownames(qqnormDT) %in% names(pipeline_geneList),]
  stopifnot(nrow(qqnormDT) == length(pipeline_geneList))
  rownames(qqnormDT) <- pipeline_geneList[rownames(qqnormDT)]
  stopifnot(setequal(pipeline_geneList, rownames(qqnormDT)))
  
  samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
  stopifnot(samp1 %in% colnames(qqnormDT))
  
  samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))
  stopifnot(samp2 %in% colnames(qqnormDT))
  
  qqnormDT_all <- qqnormDT[,c(samp1, samp2)]
  cor_qqnormMat_all <- cor(t(qqnormDT_all), method = corMethod_qqnorm)
  stopifnot(nrow(cor_qqnormMat_all) == nrow(qqnormDT))
  stopifnot(ncol(cor_qqnormMat_all) == nrow(qqnormDT))
  
  outFile <- file.path(outFold, paste0("cor_qqnormMat_all.Rdata"))
  save(cor_qqnormMat_all, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  qqnormDT_samp1 <- qqnormDT[,c(samp1)]
  cor_qqnormMat_samp1 <- cor(t(qqnormDT_samp1), method = corMethod_qqnorm)
  stopifnot(nrow(cor_qqnormMat_samp1) == nrow(qqnormDT))
  stopifnot(ncol(cor_qqnormMat_samp1) == nrow(qqnormDT))
  
  outFile <- file.path(outFold, paste0("cor_qqnormMat_samp1.Rdata"))
  save(cor_qqnormMat_samp1, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  qqnormDT_samp2 <- qqnormDT[,c(samp2)]
  cor_qqnormMat_samp2 <- cor(t(qqnormDT_samp2), method = corMethod_qqnorm)
  stopifnot(nrow(cor_qqnormMat_samp2) == nrow(qqnormDT))
  stopifnot(ncol(cor_qqnormMat_samp2) == nrow(qqnormDT))
  
  outFile <- file.path(outFold, paste0("cor_qqnormMat_samp2.Rdata"))
  save(cor_qqnormMat_samp2, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  corDT_all <- cor_qqnormMat_all
  corDT_samp1 <- cor_qqnormMat_samp1
  corDT_samp2 <- cor_qqnormMat_samp2
  
} else{
  
  outFileAll <- file.path(outFold, paste0("cor_qqnormMat_all.Rdata"))
  outFile1 <- file.path(outFold, paste0("cor_qqnormMat_samp1.Rdata"))
  outFile2 <- file.path(outFold, paste0("cor_qqnormMat_samp2.Rdata"))
  cat("... load all samples data\n")
  corDT_all <- eval(parse(text = load(outFileAll)))
  cat("... load samp1 data\n")
  corDT_samp1 <- eval(parse(text = load(outFile1)))
  cat("... load samp2 data\n")
  corDT_samp2 <- eval(parse(text = load(outFile2)))
}


stopifnot(dim(corDT_all) == dim(corDT_samp1))
stopifnot(dim(corDT_all) == dim(corDT_samp2))

stopifnot(colnames(corDT_all) == colnames(corDT_samp1))
stopifnot(colnames(corDT_all) == colnames(corDT_samp2))

stopifnot(rownames(corDT_all) == rownames(corDT_samp1))
stopifnot(rownames(corDT_all) == rownames(corDT_samp2))

stopifnot(colnames(corDT_all) == rownames(corDT_all))
stopifnot(colnames(corDT_samp2) == rownames(corDT_samp2))
stopifnot(colnames(corDT_samp1) == rownames(corDT_samp1))

cat("... lower tri all data\n")
lowerTri_all <- corDT_all[lower.tri(corDT_all, diag=F)]
cat("... lower tri samp1\n")
lowerTri_samp1 <- corDT_samp1[lower.tri(corDT_samp1, diag=F)]
cat("... lower tri samp2\n")
lowerTri_samp2 <- corDT_samp2[lower.tri(corDT_samp2, diag=F)]

outFile <- file.path(outFold, "lowerTri_all_vs_samp1.png")
png(outFile, width = myWidth, height = myHeight)
plot(x = lowerTri_all,
     y = lowerTri_samp1,
     main = paste0(curr_dataset, " - ", "all vs. samp1"),
     pch = 16, cex = 0.7)

corTest <- cor.test(lowerTri_all, lowerTri_samp1, method="pearson")
legend("topleft", paste0("PCC = ", round(corTest$estimate, 2)), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "lowerTri_all_vs_samp2.png")
png(outFile, width = myWidth, height = myHeight)
plot(x = lowerTri_all,
     y = lowerTri_samp2,
     main = paste0(curr_dataset, " - ", "all vs. samp2"),
     pch = 16, cex = 0.7)

corTest <- cor.test(lowerTri_all, lowerTri_samp2, method="pearson")
legend("topleft", paste0("PCC = ", round(corTest$estimate, 2)), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, "samp1_vs_samp2.png")
png(outFile, width = myWidth, height = myHeight)
plot(x = lowerTri_samp1,
     y = lowerTri_samp2,
     main = paste0(curr_dataset, " - ", "samp1 vs. samp2"),
     pch = 16, cex = 0.7)

corTest <- cor.test(lowerTri_samp1, lowerTri_samp2, method="pearson")
legend("topleft", paste0("PCC = ", round(corTest$estimate, 2)), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

