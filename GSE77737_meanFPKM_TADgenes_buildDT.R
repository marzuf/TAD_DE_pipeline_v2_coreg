library(foreach)
library(doMC)

# old script name: tad_genes_mean_FPKM.R
# Rscript GSE77737_meanFPKM_TADgenes_buildDT.R

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 20))

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"

# outFold <- "."
# update 10.07.2018: the outfiles are now in /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/FPKM/TCGAcrc_msi_mss
outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/FPKM", curr_dataset)
system(paste0("mkdir -p ", outFold))

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")
TADpos_DT <- read.delim(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = F)

# topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/bed_topTADs.bed")
# topTADfile <- file.path("bed_topTADs.bed")
# updated file location 10.07.2018
topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data", paste0(curr_dataset, "_topTADs10.bed"))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)

fpkmFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/FPKM")
fpkmFiles <- list.files(fpkmFold, pattern="FPKMs.txt", full.names=T, recursive = TRUE)

# HARD-CODED: check I should have 6 files: V457and V503 for MSI and V389, V410, V576, V9M for MSS
stopifnot(length(fpkmFiles) == 6)

all_DF <- foreach(exprFile = fpkmFiles) %dopar% {
  currDT <- read.delim(exprFile, header=T, stringsAsFactors = FALSE)
  currDT <- currDT[,c("gene_id","locus", "FPKM")]
  ms_status <- basename(dirname(exprFile))
  sample <- gsub("GSM.+_(.+?)_.+", "\\1", basename(exprFile))
  colnames(currDT)[colnames(currDT) == "FPKM"] <- paste0("FPKM_", sample, "_", ms_status)
  currDT
}

merged_all_DF <- Reduce(function(...) merge(..., by=c("gene_id", "locus")), all_DF)

stopifnot( length(unique(c(sapply(all_DF, nrow), nrow(merged_all_DF)))) ==1 )

merged_all_DF$gene_chr <- gsub("(chr.+):(.+)-(.+)", "\\1", merged_all_DF$locus)
merged_all_DF$gene_start <- gsub("(chr.+):(.+)-(.+)", "\\2", merged_all_DF$locus)
merged_all_DF$gene_end <- gsub("(chr.+):(.+)-(.+)", "\\3", merged_all_DF$locus)

merged_all_DF$gene_start <- as.numeric(as.character(merged_all_DF$gene_start))
merged_all_DF$gene_end <- as.numeric(as.character(merged_all_DF$gene_end))
stopifnot(!is.na(merged_all_DF$gene_start))
stopifnot(!is.na(merged_all_DF$gene_end))

cond1_columns <- colnames(merged_all_DF)[grepl(paste0("_", cond1), colnames(merged_all_DF))]
cond2_columns <- colnames(merged_all_DF)[grepl(paste0("_", cond2), colnames(merged_all_DF))]

TAD_DT <- TADpos_DT[grepl("_TAD", TADpos_DT$region),]

merged_all_DF$region <- foreach(i = seq_len(nrow(merged_all_DF)), .combine='c') %dopar% {
  curr_gene_chromo <- merged_all_DF$gene_chr[i]
  curr_gene_start <- merged_all_DF$gene_start[i]
  curr_gene_end <- merged_all_DF$gene_end[i]
  
  found_idx <- which(TADpos_DT$chromo == curr_gene_chromo &
          curr_gene_start >= TADpos_DT$start & 
          curr_gene_end <= TADpos_DT$end)
  stopifnot(length(found_idx) <= 1)
  if(length(found_idx) == 0) return("ambiguous")
  TADpos_DT$region[found_idx]
}

gene_fpkm_with_region_DT <- merged_all_DF
outFile <- file.path(outFold, "gene_fpkm_with_region_DT.Rdata")
save(gene_fpkm_with_region_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


mean_fpkm_status_DT <- foreach(i = seq_len(nrow(TAD_DT)), .combine='rbind') %dopar% {
  tad_chr <- TAD_DT$chromo[i]
  tad_start <- TAD_DT$start[i]
  tad_end <- TAD_DT$end[i]
  tad <- TAD_DT$region[i]
    
  tad_fpkmDT <- merged_all_DF[
    (as.character(merged_all_DF$gene_chr) == as.character(tad_chr)) & 
      (merged_all_DF$gene_start >= tad_start & merged_all_DF$gene_end <= tad_end),
  ]
  if(nrow(tad_fpkmDT) == 0) {
    meanSamplesGenes_cond1 <- NA
    meanSamplesGenes_cond2 <- NA
    cond1_cond2_ratio <- NA
  } else {
    meanSamplesGenes_cond1 <- mean(colMeans(tad_fpkmDT[,cond1_columns], na.rm=T), na.rm=T)
    meanSamplesGenes_cond2 <- mean(colMeans(tad_fpkmDT[,cond2_columns], na.rm=T), na.rm=T)
	cond1_cond2_ratio <- round(meanSamplesGenes_cond1/meanSamplesGenes_cond2,4)
  }
  data.frame(
    region = tad,
    start = tad_start,
    end = tad_end,
    nGenes_fpkm = nrow(tad_fpkmDT),
    mean_cond1 = meanSamplesGenes_cond1,
    mean_cond2 = meanSamplesGenes_cond2,
    cond1_cond2_ratio = cond1_cond2_ratio,
    stringsAsFactors = FALSE
  )
}

colnames(mean_fpkm_status_DT)[colnames(mean_fpkm_status_DT) == "mean_cond1"] <- paste0("mean_", cond1)
colnames(mean_fpkm_status_DT)[colnames(mean_fpkm_status_DT) == "mean_cond2"] <- paste0("mean_", cond2)
colnames(mean_fpkm_status_DT)[colnames(mean_fpkm_status_DT) == "cond1_cond2_ratio"] <- paste0(cond1, "_", cond2, "_ratio")

outFile <- file.path(outFold, "mean_fpkm_status_DT.txt")
write.table(mean_fpkm_status_DT, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

tmp <- mean_fpkm_status_DT
head(tmp)
head(topTADdt)
tmp <- tmp[match(topTADdt$region, tmp$region),]
head(tmp)
tmp <- tmp[tmp$region %in% topTADdt$region,]
write.table(tmp, col.names=T, row.names=F, sep="\t", quote=F)

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



