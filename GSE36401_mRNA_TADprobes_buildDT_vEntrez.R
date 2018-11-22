library(foreach)
library(doMC)
library(dplyr)

# old script name: tad_genes_mean_FPKM.R
# Rscript GSE36401_mRNA_TADprobes_buildDT_vEntrez.R

# in this version assign mRNA to TAD based on entrezID

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 20))

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE36401_Akhtar-Zaidi2012/mRNA_DATA_entrezID", curr_dataset)
system(paste0("mkdir -p ", outFold))

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")
TADpos_DT <- read.delim(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = F)

topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data", paste0(curr_dataset, "_topTADs10.bed"))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2tDT$entrezID <- as.character(g2tDT$entrezID)


# "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE36401_Akhtar-Zaidi2012/mRNA_DATA/TCGAcrc_msi_mss"

exprDT <- eval(parse(text = load(file.path(outFold, "entrez_exprMatrix_madnorm.Rdata"))))
# ????????????????????????????????????????????????????????????????????????? 10.07.2018 -> how to deal with NA ??????????????????????????????????????
na_rows <- apply(exprDT,1,function(x)any(is.na(x)))
cat("# of rows with NA: ", sum(na_rows), "/", length(na_rows), "\n")
dim(exprDT)
exprDT <- exprDT[!na_rows,]
dim(exprDT)
exprDT_mat <- exprDT
exprDT_mat[1:5,1:5]
# ?????????????????????????????????????????????????????????????????????????

# NUMBER OF GENES WITH G-2-T ASSIGNMENT

cat("# of genes with g2t assignment: ", sum(rownames(exprDT_mat) %in% g2tDT$entrezID), "/", nrow(exprDT_mat), "\n")
# 16054/16106

exprDT_mat <- exprDT_mat[rownames(exprDT_mat) %in% g2tDT$entrezID, ]
nrow(exprDT_mat)

# MAP THE PROBES TO THE TADs
exprDT <- as.data.frame(exprDT_mat)
sampleNames <- colnames(exprDT)
exprDT$entrezID <- rownames(exprDT)

exprDT <- left_join(exprDT, g2tDT, by = "entrezID")
stopifnot(!any(is.na(exprDT)))

cond1_columns <- sampleNames[grepl(paste0("_", cond1), sampleNames)]
cond2_columns <- sampleNames[grepl(paste0("_", cond2), sampleNames)]
stopifnot(length(cond1_columns) > 0)
stopifnot(length(cond2_columns) > 0)

TAD_DT <- TADpos_DT[grepl("_TAD", TADpos_DT$region),]

exprDT$region <- foreach(i = seq_len(nrow(exprDT)), .combine='c') %dopar% {
  curr_entrez_chromo <- exprDT$chromo[i]
  curr_entrez_start <- exprDT$start[i]
  curr_entrez_end <- exprDT$end[i]
  found_idx <- which(TADpos_DT$chromo == curr_entrez_chromo &
          curr_entrez_start >= TADpos_DT$start & 
          curr_entrez_end <= TADpos_DT$end)
  stopifnot(length(found_idx) <= 1)
  if(length(found_idx) == 0) return("ambiguous")
  TADpos_DT$region[found_idx]
}

mRNA_with_region_DT_vEntrez <- exprDT

outFile <- file.path(outFold, "mRNA_with_region_DT_vEntrez.Rdata")
save(mRNA_with_region_DT_vEntrez, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

mean_mRNA_status_DT_vEntrez <- foreach(i = seq_len(nrow(TAD_DT)), .combine='rbind') %dopar% {
  tad_chr <- TAD_DT$chromo[i]
  tad_start <- TAD_DT$start[i]
  tad_end <- TAD_DT$end[i]
  tad <- TAD_DT$region[i]
    
  tad_mRNA_DT <- exprDT[
    (as.character(exprDT$chromo) == as.character(tad_chr)) & 
      (exprDT$start >= tad_start & exprDT$end <= tad_end),
  ]
  if(nrow(tad_mRNA_DT) == 0) {
    meanSamplesGenes_cond1 <- NA
    meanSamplesGenes_cond2 <- NA
    cond1_cond2_ratio <- NA
  } else {
    meanSamplesGenes_cond1 <- mean(colMeans(tad_mRNA_DT[,cond1_columns], na.rm=T), na.rm=T)
    meanSamplesGenes_cond2 <- mean(colMeans(tad_mRNA_DT[,cond2_columns], na.rm=T), na.rm=T)
	cond1_cond2_ratio <- round(meanSamplesGenes_cond1/meanSamplesGenes_cond2,4)
  }
  data.frame(
    region = tad,
    start = tad_start,
    end = tad_end,
    nGenes_mRNA = nrow(tad_mRNA_DT),
    mean_cond1 = meanSamplesGenes_cond1,
    mean_cond2 = meanSamplesGenes_cond2,
    cond1_cond2_ratio = cond1_cond2_ratio,
    stringsAsFactors = FALSE
  )
}

colnames(mean_mRNA_status_DT_vEntrez)[colnames(mean_mRNA_status_DT_vEntrez) == "mean_cond1"] <- paste0("mean_", cond1)
colnames(mean_mRNA_status_DT_vEntrez)[colnames(mean_mRNA_status_DT_vEntrez) == "mean_cond2"] <- paste0("mean_", cond2)
colnames(mean_mRNA_status_DT_vEntrez)[colnames(mean_mRNA_status_DT_vEntrez) == "cond1_cond2_ratio"] <- paste0(cond1, "_", cond2, "_ratio")

outFile <- file.path(outFold, "mean_mRNA_status_DT_vEntrez.txt")
write.table(mean_mRNA_status_DT_vEntrez, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

tmp <- mean_mRNA_status_DT_vEntrez
head(tmp)
head(topTADdt)
tmp <- tmp[match(topTADdt$region, tmp$region),]
head(tmp)
tmp <- tmp[tmp$region %in% topTADdt$region,]
write.table(tmp, col.names=T, row.names=F, sep="\t", quote=F)

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



