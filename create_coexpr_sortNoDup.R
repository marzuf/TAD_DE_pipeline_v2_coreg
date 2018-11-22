startTime <- Sys.time()
cat(paste0("> Rscript create_coexpr_sortNoDup.R\n"))

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_coexpr_sortNoDup.R TCGAcrc_msi_mss
# Rscript create_coexpr_sortNoDup.R TCGAstad_EBVneg_EBVpos
# Rscript create_coexpr_sortNoDup.R GSE102073_stic_nostic
# Rscript create_coexpr_sortNoDup.R GSE74927_neg_pos
# Rscript create_coexpr_sortNoDup.R GSE65540_before_after
# Rscript create_coexpr_sortNoDup.R GSE84231_lhb_rhb
# Rscript create_coexpr_sortNoDup.R GSE86356_tibMD1_tibNorm

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the expression table to ensure alphabetical order of the genes
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> replace diag + upper.tri with NA, use melt with na.rm=TRUE
#    so that data saved are smaller !!!

caller ="TopDom"
args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- args[1]
corMethod <- "pearson"


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_COEXPR_SORTNODUP", paste0(curr_dataset, "_", corMethod))
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

################################################################ UPDATE 30.06.2018 -> sort the rownames to have gene1 < gene2
qqnormDT_sorted <- qqnormDT[sort(rownames(qqnormDT)),]
stopifnot(dim(qqnormDT_sorted) == dim(qqnormDT))
stopifnot(setequal(rownames(qqnormDT), rownames(qqnormDT_sorted)))
qqnormDT <- qqnormDT_sorted
qqnormDT_tmp <- qqnormDT

cor_qqnormMat <- cor(t(qqnormDT), method = corMethod)
stopifnot(nrow(cor_qqnormMat) == nrow(qqnormDT))
stopifnot(ncol(cor_qqnormMat) == nrow(qqnormDT))

###### UPDATE 30.06.2018
# coexprDT <- melt(cor_qqnormMat)
# colnames(coexprDT) <- c("gene1", "gene2", "coexpr")
# coexprDT <- coexprDT[coexprDT$gene1 != coexprDT$gene2,]

cor_qqnormMat_NA <- cor_qqnormMat
cor_qqnormMat_NA[lower.tri(cor_qqnormMat_NA, diag=T)] <- NA
coexprDT <- melt(cor_qqnormMat_NA, na.rm = T)
colnames(coexprDT) <- c("gene1", "gene2", "coexpr")
coexprDT$gene1 <- as.character(coexprDT$gene1)
coexprDT$gene2 <- as.character(coexprDT$gene2)
stopifnot(coexprDT$gene1 < coexprDT$gene2)

outFile <- file.path(outFold, "coexprDT.Rdata")
save(coexprDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
