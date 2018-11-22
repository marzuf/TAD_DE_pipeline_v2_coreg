SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

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


# outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/topTADs_average_GSE_data/BOXPLOT_FPKMs")
# this was the folder for the older version fpkm_topTADs.R
outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/GSE77737_TOPTADs_FPKM_BOXPLOT", curr_dataset)
system(paste0("mkdir -p ", outFold))

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


for(tad in top_tads) {
  curr_df <- fpkmDT[fpkmDT$region == tad, c("gene_id", cond1_columns, cond2_columns)]
  
  if(nrow(curr_df) == 0) {
    cat("! No probe mapping to topTAD: ", tad, "!\n")
    next
  }
  
  varying_cols <- c(which(colnames(curr_df) %in% cond1_columns), which(colnames(curr_df) %in% cond2_columns))
  curr_df_m <- reshape(curr_df, idvar="gene_id", direction = "long", varying=2:7,  timevar="sample", sep="")
  curr_df_m$status <- gsub(".+_(.+)", "\\1", curr_df_m$sample)
  colnames(curr_df_m)[3] <- "FPKM"
  curr_df_m$FPKM_log10 <- log10(curr_df_m$FPKM+0.001)
  
  subTit <- paste0("nGenes in TAD = ", length(unique(curr_df_m$gene_id)), "/", nrow(fpkmDT), 
                   "\n# ", cond1, " = ", length(cond1_columns), 
                   " - # ", cond2, " = ", length(cond2_columns))
  
  for(countType in plot_countType) {
    
    p <- ggboxplot(curr_df_m, 
                   x = "status",
                   xlab = paste0("CRC status"), 
                   y = paste0("FPKM", countType),
                   ylab = paste0("FPKM", countType),
                   legend.title="",
                   # legend="right",
                   legend="none",
                   title = paste0(tad, ": FPKM"),
                   subtitle=subTit,
                   fill = "status", palette = c("steelblue3", "tan2"),
                   add = "jitter")
    p <- p + font("legend.text", size = 14)
    p <- p + font("subtitle", face="bold")
    p <- p + font("title", size = 18, face="bold")
    p <- p + font("xlab", size = 16, face="bold")
    p <- p + font("ylab", size = 16)
    p <- p + font("xy.text", size = 14)
    p <- p + theme(plot.title = element_text(hjust=0.5))
    p <- p + stat_compare_means(label.x=0.5)
    if(SSHFS) p
    
    outFile <- file.path(outFold, paste0(tad, "_FPKM", countType, ".", plotType))
    ggsave(plot=p, filename = outFile, height=myHeight, width=myWidth )
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}
  
