startTime <- Sys.time()
cat(paste0("> Rscript create_coexpr.R\n"))

stop("!!! use _sortNoDup script !!!\n")

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_coexpr.R TCGAcrc_msi_mss
# Rscript create_coexpr.R GSE102073_stic_nostic
# Rscript create_coexpr.R GSE74927_neg_pos
# Rscript create_coexpr.R GSE65540_before_after -> ok
# Rscript create_coexpr.R GSE84231_lhb_rhb -> ok
# Rscript create_coexpr.R GSE86356_tibMD1_tibNorm -> running

caller ="TopDom"
args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- args[1]
corMethod <- "pearson"


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_COEXPR", paste0(curr_dataset, "_", corMethod))
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

cor_qqnormMat <- cor(t(qqnormDT), method = corMethod)
stopifnot(nrow(cor_qqnormMat) == nrow(qqnormDT))
stopifnot(ncol(cor_qqnormMat) == nrow(qqnormDT))

coexprDT <- melt(cor_qqnormMat)
colnames(coexprDT) <- c("gene1", "gene2", "coexpr")
coexprDT <- coexprDT[coexprDT$gene1 != coexprDT$gene2,]

outFile <- file.path(outFold, "coexprDT.Rdata")
save(coexprDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
