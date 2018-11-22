### LOCATION OF DATA FILES ###########################################################################################################################
# -> probe genomic coordinates (hg19): /mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/hm450_probe_coordinates.txt
# - contains all the 450k probes
# - do not use gene annotation from this file
# 
# -> /mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/HM450K_probe_FANTOM5_promoter_annotation_table.csv 
# - marco & sadeq derived mapping gene annotation for each probe
# - gene annotation good to use
# 
# -> /mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/meth/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv
# - beta values methylation file
# 
# /mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/gene_expression/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv  
# 
# Annotation file: /mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/pancanAtlas_pathways/pancanAtlas_pathway_caselist_2017.02.23.txt

# TCGA researchers identify four subtypes of stomach cancer
# https://www.cancer.gov/news-events/press-releases/2014/TCGAgastric
# 1) Tumors in the first group, which represented 9 percent of the tumors, were positive for Epstein-Barr virus (EBV) and had several other molecular commonalities.
# 2) Tumors in a second subgroup (22 percent of the tumors) had high microsatellite instability (MSI), which is the tendency for mutations to accumulate in repeated sequences of DNA. 
# The remaining subgroups differed in the level of somatic copy number alterations (SCNAs), which can result from duplication or deletion of sections of the genome. 
# 3) The tumors in the third subgroup, which comprised 20 percent of the tumors, were considered to have a low level of SCNAs and were called genomically stable. 
# 4) The remaining 50 percent of tumors were classified as chromosomally unstable, with a high level of SCNAs.


#######################################################################
printAndLog <- function(txt = NULL, logFile = NULL) {
  stopifnot(!is.null(txt))
  cat(txt)
  if(!is.null(logFile)) 
    cat(txt, file = logFile, append=T)
}


#######################################################################################################################################################

# Rscript 1_stomach_EBV_neg_pos.R
# everything hard-coded

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

cat(paste0("> START ", "prep_crc_methylation.R",  "\n"))

outFold <- "PREP_CRC_METHYLATION"
system(paste0("mkdir -p ", outFold))

sampType <- "COAD"
# msi_subtype <- "HM-SNV"
# mss_subtype <- "GS"

methylationFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/meth/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv")
annotFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/pancanAtlas_pathways/pancanAtlas_pathway_caselist_2017.02.23.txt")

logFile <- file.path(outFold, "prep_crc_methylation.txt")
system(paste0("rm -f ", logFile))


txt <- paste0("!!! HARD-CODED !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... sampType =\t", sampType, "\n")
printAndLog(txt, logFile)
# txt <- paste0("... msi_subtype =\t", msi_subtype, "\n")
# printAndLog(txt, logFile)
# txt <- paste0("... mss_subtype =\t", mss_subtype, "\n")
# printAndLog(txt, logFile)

# Retrieve the samples
annotDT <- read.delim(annotFile, stringsAsFactors = F, header=T)
annotDT <- annotDT[annotDT$type == sampType,]
unique(annotDT$subtype)
# "CIN"      "GS"       "HM-indel"      "HM-SNV" 
# CIN = chromosomal instability; 
# GS = genomic stable 
# HM-indel = hypermutated indels
# HM-SNV = hypermutated somatic number variation
# do not take samples with code >= 10

# do we have normal samples ???

sample_codes <- as.numeric(gsub("^.+-.+-.+-(.+)$", "\\1", annotDT$sample_id))
unique(sample_codes)
normal_tumor <- gsub("^.+-.+-.+-(.+)$", "\\1", annotDT$sample_id)
unique(normal_tumor)

annotDT <- annotDT[sample_codes <  10,]

stopifnot(annotDT$type == sampType)

load("../../other_datasets/TCGA_crc_msi_mss/28.07_prepData/msi_ID.Rdata")
load("../../other_datasets/TCGA_crc_msi_mss/28.07_prepData/mss_ID.Rdata")
msi_ID <- gsub("\\.", "-",msi_ID)
mss_ID <- gsub("\\.", "-",mss_ID)
unique(annotDT$subtype[annotDT$sample_id %in% msi_ID])  ## HM-indel, CIN ??????????????????????????????
unique(annotDT$type[annotDT$sample_id %in% msi_ID])  ## HM-indel, CIN ??????????????????????????????
table(annotDT$subtype[annotDT$sample_id %in% msi_ID])
# CIN HM-indel 
# 1       38 
# tot: 47
table(annotDT$subtype[annotDT$sample_id %in% mss_ID])
# CIN     GS HM-SNV 
# 125     24      2 
# tot: 227
unique(annotDT$subtype[annotDT$sample_id %in% mss_ID])  ## HM-indel, CIN, GS, HM-SNV  ??????????????????????????????
unique(annotDT$type[annotDT$sample_id %in% mss_ID])  ## HM-indel, CIN, GS, HM-SNV  ??????????????????????????????

msi_id <- msi_ID[msi_ID %in% annotDT$sample_id]
txt <- paste0("... found msi_id:\t", length(msi_id), "/", length(msi_ID), "\n")
printAndLog(txt, logFile)

mss_id <- mss_ID[mss_ID %in% annotDT$sample_id]
txt <- paste0("... found mss_id:\t", length(mss_id), "/", length(mss_ID), "\n")
printAndLog(txt, logFile)

# msi_id <- annotDT$sample_id[annotDT$subtype == msi_subtype & annotDT$type == sampType]
# mss_id <- annotDT$sample_id[annotDT$subtype == mss_subtype & annotDT$type == sampType]
save(msi_id, file = file.path(outFold, "msi_id.Rdata"))
save(mss_id, file = file.path(outFold, "mss_id.Rdata"))

if(file.exists(file.path(outFold, "msi_id.Rdata"))){
  cat(paste0("...written: ", file.path(outFold, "msi_id.Rdata"), "\n"))
}else{
  stop("error\n")
}

if(file.exists(file.path(outFold, "mss_id.Rdata"))){
  cat(paste0("...written: ", file.path(outFold, "mss_id.Rdata"), "\n"))
}else{
  stop("error\n")
}

#cat("... load TCGA rnaseq \t")
#TCGA_rnaseqFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/gene_expression/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv")
#TCGA_DT <- read.delim(TCGA_rnaseqFile, header=T, stringsAsFactors=F)
#cat("done\n")
#colnames(TCGA_DT) <- substr(x = colnames(TCGA_DT), start = 1, stop = 15)
#colnames(TCGA_DT) <- gsub("\\.", "-", colnames(TCGA_DT))
#stopifnot(any(c(msi_id, mss_id) %in% colnames(TCGA_DT)))
#inSTAD_id = colnames(TCGA_DT)[colnames(TCGA_DT) %in% c(msi_id, mss_id)]
#STAD_RNAseq_DT <- TCGA_DT[,c("gene_id", inSTAD_id)]
#save(STAD_RNAseq_DT, file = file.path(outFold, "STAD_RNAseq_DT.Rdata"))
#if(file.exists(file.path(outFold, "STAD_RNAseq_DT.Rdata"))){
#cat(paste0("...written: ", file.path(outFold, "STAD_RNAseq_DT.Rdata"), "\n"))
#}else{
#stop("error\n")
#}

cat("... load TCGA methylation \t")
methylation_DT <- read.delim(methylationFile, header=T, row.names = 1, stringsAsFactors = F)
colnames(methylation_DT) <- substr(x = colnames(methylation_DT), start = 1, stop = 15)
colnames(methylation_DT) <- gsub("\\.", "-", colnames(methylation_DT))
stopifnot(any(c(msi_id, mss_id) %in% colnames(methylation_DT)))
inMet_id = colnames(methylation_DT)[colnames(methylation_DT) %in% c(msi_id, mss_id)]
COAD_methylation_DT <- methylation_DT[,inMet_id]
cat("done\n")
save(COAD_methylation_DT, file = file.path(outFold, "COAD_methylation_DT.Rdata"))

if(file.exists(file.path(outFold, "COAD_methylation_DT.Rdata"))){
  cat(paste0("...written: ", file.path(outFold, "COAD_methylation_DT.Rdata"), "\n"))
}else{
  stop("error\n")
}

cat("*** DONE\n")
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)


# methylProbeAnnotFile <- file.path(setDir, "/mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/hm450_probe_coordinates.txt")
# methylGeneAnnotFile <- file.path(setDir, "/mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/HM450K_probe_FANTOM5_promoter_annotation_table.csv")


