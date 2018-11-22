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
outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/CHECK_FPKM_MAPPING")
system(paste0("mkdir -p ", outFold))

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")
TADpos_DT <- read.delim(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = F)

# topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/bed_topTADs.bed")
# topTADfile <- file.path("bed_topTADs.bed")
# updated file location 10.07.2018
topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data", paste0(curr_dataset, "_topTADs10.bed"))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)

gene2tadDT_file <- file.path(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2tDT$entrezID <- as.character(g2tDT$entrezID)

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

nrow(merged_all_DF)
# 25207

cond1_columns <- colnames(merged_all_DF)[grepl(paste0("_", cond1), colnames(merged_all_DF))]
cond2_columns <- colnames(merged_all_DF)[grepl(paste0("_", cond2), colnames(merged_all_DF))]

TAD_DT <- TADpos_DT[grepl("_TAD", TADpos_DT$region),]



### MAP BY POSITION:
merged_all_DF$region <- foreach(i = seq_len(nrow(merged_all_DF)), .combine='c') %dopar% {
  curr_gene_chromo <- merged_all_DF$gene_chr[i]
  curr_gene_start <- merged_all_DF$gene_start[i]
  curr_gene_end <- merged_all_DF$gene_end[i]
  
  # found_idx <- which(TADpos_DT$chromo == curr_gene_chromo &
  #                      curr_gene_start >= TADpos_DT$start & 
  #                      curr_gene_end <= TADpos_DT$end)
  # less stringent:
  found_idx <- which(TADpos_DT$chromo == curr_gene_chromo &
                       curr_gene_start >= TADpos_DT$start & 
                       curr_gene_start <= TADpos_DT$end)
  
  stopifnot(length(found_idx) <= 1)
  if(length(found_idx) == 0) return("ambiguous")
  TADpos_DT$region[found_idx]
}

gene_fpkm_with_region_DT <- merged_all_DF

gene_fpkm_with_region_DT_posMap <- gene_fpkm_with_region_DT

outFile <- file.path(outFold, "gene_fpkm_with_region_DT_posMap.Rdata")
save(gene_fpkm_with_region_DT_posMap, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

### MAP BY GENE SYMBOL:
gene2tad_DT <- g2tDT
entrezDT$entrezID <- as.character(entrezDT$entrezID)
entrezDT$symbol <- as.character(entrezDT$symbol)
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

merged_all_DF$region <- foreach(i = seq_len(nrow(merged_all_DF)), .combine='c') %dopar% {
  
  curr_gene <- as.character(merged_all_DF$gene_id[i])

  curr_entrezID <- entrezDT$entrezID[entrezDT$symbol == curr_gene]
  
  if(length(curr_entrezID) == 0){
    return(NA)  
  }  
  
  curr_region <- gene2tad_DT$region[gene2tad_DT$entrezID == curr_entrezID]
  
  if(length(curr_region) == 0){
    return(NA)  
  } else {
    return(curr_region)
  }
}

gene_fpkm_with_region_DT_symbolMap <- merged_all_DF
outFile <- file.path(outFold, "gene_fpkm_with_region_DT_symbolMap.Rdata")
save(gene_fpkm_with_region_DT_symbolMap, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/CHECK_FPKM_MAPPING/gene_fpkm_with_region_DT_posMap.Rdata")
load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/CHECK_FPKM_MAPPING/gene_fpkm_with_region_DT_symbolMap.Rdata")
ncol(gene_fpkm_with_region_DT_posMap)

posDT <- gene_fpkm_with_region_DT_posMap
colnames(posDT) <- paste0(colnames(posDT), "_posMap")

symbolDT <- gene_fpkm_with_region_DT_symbolMap
colnames(symbolDT) <- paste0(colnames(symbolDT), "_symbolMap")

cmpDT <- cbind(posDT[, c("gene_id_posMap", "region_posMap")], symbolDT[,c("gene_id_symbolMap", "region_symbolMap")])
stopifnot(cmpDT$gene_id_posMap == cmpDT$gene_id_symbolMap)

any(is.na(cmpDT$region_posMap))

sum(grepl("_TAD", cmpDT$region_posMap))
# 20148
sum(grepl("_TAD", cmpDT$region_symbolMap))
# 20254
sum(is.na(cmpDT$region_symbolMap))
#1324
sum(is.na(cmpDT$region_symbolMap))

bothDT <- na.omit(cmpDT)
bothDT <- bothDT[grepl("_TAD", bothDT$region_posMap) & grepl("_TAD", bothDT$region_symbolMap),]

nrow(bothDT)
# 19193
sum(bothDT$region_posMap == bothDT$region_symbolMap)
# 18394

gcolnames(gene_fpkm_with_region_DT_posMap)[colnames(gene_fpkm_with_region_DT_posMap) == "region"] <- "region_posMap"
colnames(gene_fpkm_with_region_DT_symbolMap)[colnames(gene_fpkm_with_region_DT_symbolMap) == "region"] <- "region_symbolMap"
head(gene_fpkm_with_region_DT_posMap)




cmpDT <- merge(gene_fpkm_with_region_DT_posMap[, c("gene_id", "region_posMap")], 
               gene_fpkm_with_region_DT_symbolMap[, c("gene_id", "region_symbolMap")],
               by = "gene_id")

require(reshape2)
cmpDT <- merge(gene_fpkm_with_region_DT_posMap[, c("gene_id", "region_posMap")], 
               gene_fpkm_with_region_DT_symbolMap[, c("gene_id", "region_symbolMap")],
               by = "gene_id", all.x=TRUE, all.y=TRUE)


head(cmpDT)
nrow(cmpDT)
nrow(gene_fpkm_with_region_DT_posMap)
nrow(gene_fpkm_with_region_DT_symbolMap)


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


