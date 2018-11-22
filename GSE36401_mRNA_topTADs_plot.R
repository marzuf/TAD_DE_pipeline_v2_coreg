SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

library(foreach)
library(doMC)
library(ggpubr)

registerDoMC(ifelse(SSHFS, 20, 2))

# Rscript GSE36401_mRNA_topTADs_plot.R

plot_countType <- c("", "_log10")


plotType <- "png"
myHeight <- 7
myWidth <- 10

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"

# HARD CODED, they should be 2 MSI and 4 MSS !
check_nCond1 <- 4
check_nCond2 <- 5

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE36401_Akhtar-Zaidi2012/GSE36401_mRNA_TOPTADs_BOXPLOT_v2", curr_dataset)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "mRNA_topTADs_boxplot_logFile.txt")
system(paste0("rm -f ", logFile))

topTADfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data", paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADfile))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)

mRNAfile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE36401_Akhtar-Zaidi2012/mRNA_DATA/", curr_dataset, "/probe_mRNA_with_region_DT.Rdata"))
stopifnot(file.exists(mRNAfile))
mrnaDT <- eval(parse(text = load(mRNAfile)))
mrnaDT$region <- as.character(mrnaDT$region)

cond1_columns <- colnames(mrnaDT)[grepl(paste0("_", cond1), colnames(mrnaDT))]
cond2_columns <- colnames(mrnaDT)[grepl(paste0("_", cond2), colnames(mrnaDT))]

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
  
  cat(paste0("*** ", tad, "*** \n"), file = logFile)
  
  curr_df <- mrnaDT[mrnaDT$region == tad, c("probe", cond1_columns, cond2_columns)]
  if(nrow(curr_df) == 0) {
    cat("! No probe mapping to topTAD: ", tad, "!\n")
    next
  }
  varying_cols <- c(which(colnames(curr_df) %in% cond1_columns), which(colnames(curr_df) %in% cond2_columns))
  curr_df_m <- reshape(curr_df, idvar="probe", direction = "long", varying=varying_cols,  timevar="sample", sep="")
  curr_df_m$status <- gsub(".+_(.+)", "\\1", curr_df_m$sample)
  colnames(curr_df_m)[3] <- "mRNA"
  curr_df_m$mRNA_log10 <- log10(curr_df_m$mRNA+0.001)
  
  write.table(curr_df_m, row.names = F, col.names = T, append = T, file = logFile)

  subTit <- paste0("nProbes in TAD = ", length(unique(curr_df_m$probe)), "/", nrow(mrnaDT), 
                   "\n# ", cond1, " = ", length(cond1_columns), 
                   " - # ", cond2, " = ", length(cond2_columns))
  
  for(countType in plot_countType) {
    
    p <- ggboxplot(curr_df_m, 
                   x = "status",
                   xlab = paste0("CRC status"), 
                   y = paste0("mRNA", countType),
                   ylab = paste0("mRNA", countType),
                   legend.title="",
                   # legend="right",
                   legend="none",
                   title = paste0(tad, ": mRNA"),
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
    
    outFile <- file.path(outFold, paste0(tad, "_mRNA", countType, ".", plotType))
    ggsave(plot=p, filename = outFile, height=myHeight, width=myWidth )
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}
  
