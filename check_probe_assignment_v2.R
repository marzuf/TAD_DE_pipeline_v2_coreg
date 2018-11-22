SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggpubr)

registerDoMC(ifelse(SSHFS, 20, 2))

# Rscript check_probe_assignment_v2.R

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"

# HARD CODED, they should be 2 MSI and 4 MSS !
check_nCond1 <- 4
check_nCond2 <- 5

outFold <- "check_probe_assignment_v2"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "check_probe_topTADs.txt")
system(paste0("rm -f ", logFile))

topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data", paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADfile))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)

mRNAfile <- file.path(setDir, 
                        paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data"),
                        paste0("GSE36401_Akhtar-Zaidi2012/mRNA_DATA/", curr_dataset, "/probe_mRNA_with_region_DT.Rdata"))
stopifnot(file.exists(mRNAfile))
mrnaDT <- eval(parse(text = load(mRNAfile)))
mrnaDT$region <- as.character(mrnaDT$region)

mRNAfile_vEntrez <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data",
                      "GSE36401_Akhtar-Zaidi2012/mRNA_DATA_entrezID/TCGAcrc_msi_mss/mRNA_with_region_DT_vEntrez.Rdata")
stopifnot(file.exists(mRNAfile_vEntrez))
mrnaDT_vEntrez <- eval(parse(text = load(mRNAfile_vEntrez)))
mrnaDT_vEntrez$region <- as.character(mrnaDT_vEntrez$region)

gene2tadDT_file <- file.path(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2tDT$entrezID <- as.character(g2tDT$entrezID)

dataset_geneListFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER",
                                  curr_dataset, "0_prepGeneData", "pipeline_geneList.Rdata")
stopifnot(file.exists(dataset_geneListFile))
geneList <- eval(parse(text = load(dataset_geneListFile)))

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)

probAssign_file <- file.path(setDir, 
                             "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE36401_Akhtar-Zaidi2012/mRNA_DATA_entrezID/TCGAcrc_msi_mss_v1/probeIDentrezID_DT.Rdata")
probAssign_DT <- eval(parse(text = load(probAssign_file)))

cond1_columns <- colnames(mrnaDT)[grepl(paste0("_", cond1), colnames(mrnaDT))]
cond2_columns <- colnames(mrnaDT)[grepl(paste0("_", cond2), colnames(mrnaDT))]

if(curr_dataset == "TCGAcrc_msi_mss") {
  stopifnot(length(cond1_columns) == check_nCond1)
  stopifnot(length(cond2_columns) == check_nCond2)
} else{
  stop("not available\n")
}

get_symbol <- function(x_entrez, entrez_table) {
  entrez_table$symbol[as.character(entrez_table$entrezID) == as.character(x_entrez)]
}


top_tads <- unique(as.character(topTADdt$region))

# my_comparisons <- list( c("MSI", "MSS"))
# p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 2)                   # Add global p-value

top_tads <- c("chr1_TAD150", "chr12_TAD81", "chr6_TAD58")

for(tad in top_tads) {
  
  cat(paste0("\n*** ", tad, "\t", 
             as.character(topTADdt$chromo[topTADdt$region == tad]), "\t", 
             as.character(topTADdt$start[topTADdt$region == tad]), "\t", 
             as.character(topTADdt$end[topTADdt$region == tad]),
             "***\n"), file = logFile, append = TRUE)
  
  curr_df <- mrnaDT[mrnaDT$region == tad, c("probe", "probe_chromo", "probe_start", "probe_end")]
  head(curr_df)
  
  curr_df_withEntrez <- merge(curr_df, probAssign_DT[,c("probeID", "entrezID")], by.x="probe", by.y="probeID")
  curr_df_withEntrez <- merge(curr_df_withEntrez, entrezDT[,c("entrezID", "symbol")], by="entrezID")
  curr_df_withEntrez$TAD <- tad
  curr_df_withEntrez <- curr_df_withEntrez[,c("TAD", "probe", "probe_chromo", "probe_start", "probe_end", "entrezID", "symbol")]
  
  cat(paste0("> TABLE LIFTOVER PROBE MAPPING\n"), file = logFile, append = TRUE)
  write.table(curr_df_withEntrez, file = logFile, append=T, quote=F, row.names = F, col.names = T, sep="\t")
  
  curr_df_vE <- mrnaDT_vEntrez[mrnaDT_vEntrez$region == tad, c("entrezID", "chromo", "start", "end")]
  head(curr_df_vE)
  
  # curr_df_vE_withEntrez <- merge(curr_df_vE, probAssign_DT[,c("probeID", "entrezID")], by.x="probe", by.y="probeID")
  curr_df_vE_withEntrez <- merge(curr_df_vE, entrezDT[,c("entrezID", "symbol")], by="entrezID")
  curr_df_vE_withEntrez$TAD <- tad
  curr_df_vE_withEntrez <- curr_df_vE_withEntrez[,c("TAD", "chromo", "start", "end", "entrezID", "symbol")]
  
  cat(paste0("> TABLE ENTREZ ASSIGNMENT\n"), file = logFile, append = TRUE)
  write.table(curr_df_vE_withEntrez, file = logFile, append=T, quote=F, row.names = F, col.names = T, sep="\t")
  
  
  g2t_genes <- g2tDT$entrezID[g2tDT$region == tad]
  dataset_genes <- pipeline_geneList[pipeline_geneList %in% g2t_genes]
  curr_df_dataset <- data.frame(
    TAD = tad,
    entrezID = as.character(dataset_genes),
    symbol = as.character(sapply(dataset_genes, get_symbol, entrezDT)),
    stringsAsFactors = FALSE
    
  )
  cat(paste0("> TABLE DATASET PIPELINE GENES\n"), file = logFile, append = TRUE)
  write.table(curr_df_dataset, file = logFile, append=T, quote=F, row.names = F, col.names = T, sep="\t")
  
  
}
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))
