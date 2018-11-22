SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggpubr)

registerDoMC(ifelse(SSHFS, 20, 2))

# Rscript GSE77737_FPKM_topTADs_plot.R

plot_countType <- c("", "_log10")

plotType <- "png"
myHeight <- 7
myWidth <- 10

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"

# HARD CODED, they should be 2 MSI and 4 MSS !
check_nCond1 <- 2
check_nCond2 <- 4

outFold <- "check_probe_assignment_v2"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "check_probe_topTADs_FPKM.txt")
system(paste0("rm -f ", logFile))

topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data", paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADfile))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)

# load the file created with FPKM data for all samples (output file from: GSE77737_meanFPKM_TADgenes_buildDT.R)
fpkmDTfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/FPKM", curr_dataset,  "gene_fpkm_with_region_DT.Rdata")
stopifnot(file.exists(fpkmDTfile))
fpkmDT <- eval(parse(text = load(fpkmDTfile)))
fpkmDT$region <- as.character(fpkmDT$region)

cond1_columns <- colnames(fpkmDT)[grepl(paste0("_", cond1), colnames(fpkmDT))]
cond2_columns <- colnames(fpkmDT)[grepl(paste0("_", cond2), colnames(fpkmDT))]


entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)


# HARD CODED, they should be 2 MSI and 4 MSS !
if(curr_dataset == "TCGAcrc_msi_mss") {
  stopifnot(length(cond1_columns) == check_nCond1)
  stopifnot(length(cond2_columns) == check_nCond2)
} else{
  stop("not available\n")
}

top_tads <- unique(as.character(topTADdt$region))

# my_comparisons <- list( c("MSI", "MSS"))
# p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 2)                   # Add global p-value



get_symbol <- function(x_gene, entrez_table) {
  x <- entrez_table$entrezID[as.character(entrez_table$symbol) == as.character(x_gene)]
  if(length(x) == 0) {
    x <- NA
  }
  return(x)
}

top_tads <- c("chr1_TAD150", "chr12_TAD81", "chr6_TAD58")


for(tad in top_tads) {
  
  cat(paste0("\n*** ", tad, "\t", 
             as.character(topTADdt$chromo[topTADdt$region == tad]), "\t", 
             as.character(topTADdt$start[topTADdt$region == tad]), "\t", 
             as.character(topTADdt$end[topTADdt$region == tad]),
             "***\n"), file = logFile, append = TRUE)
  
  
  curr_df <- fpkmDT[fpkmDT$region == tad, c("gene_id", "locus", cond1_columns, cond2_columns)]
  
  
  if(nrow(curr_df) == 0) {
    cat("! No probe mapping to topTAD: ", tad, "!\n")
    next
  }
  
  curr_df <- merge(curr_df, entrezDT[,c("symbol", "entrezID", "chromo", "start", "end")], by.x = "gene_id", by.y="symbol", all.x = T, all.y = F)
  #curr_df$entrezID <- sapply(curr_df$gene_id, function(x) get_symbol(x, entrezDT))
  curr_df <- curr_df[,c("entrezID", "gene_id", "locus","chromo", "start", "end")]
  colnames(curr_df)[2] <- "symbol"
  
  cat(paste0("> TABLE FPKM\n"), file = logFile, append = TRUE)
  write.table(curr_df, file = logFile, append=T, quote=F, row.names = F, col.names = T, sep="\t")
  
  
}
  
  
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))

