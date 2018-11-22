SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

library(foreach)
library(doMC)
library(ggpubr)

# Rscript GSE77737_histone_marks_peaksChr6TAD58.R

registerDoMC(ifelse(SSHFS, 20, 2))

curr_dataset <- "TCGAcrc_msi_mss"
cond1 <- "MSI"
cond2 <- "MSS"

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/GSE77737_PEAKSCHR6TAD58_aggregHISTMARKS_BOXPLOT", curr_dataset)
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- 7
myWidth <- 10

# averageFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/peaksChr6TAD58_average_GSE_data")
# averageFiles <- list.files(averageFolder, pattern="_peaksChr6TAD58.tab$", full.names = TRUE)
# MSI_files <- averageFiles[grep(paste0(cond1, "_"), basename(averageFiles))]
# MSS_files <- averageFiles[grep(paste0(cond2, "_"), basename(averageFiles))]

# UPDATE 10.07.2018
mainInFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE77737_peaksChr6TAD58_avgHistMarks")
# # /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE77737_peaksChr6TAD58_avgHistMarks/MSS/H3K4me1/GSM2058089_V968_MSS_H3K4me1_peaksChr6TAD58.tab
cond1_files <- list.files(file.path(mainInFold, cond1), recursive = TRUE, pattern = "_peaksChr6TAD58.tab", full.names = T)
cond2_files <- list.files(file.path(mainInFold, cond2), recursive = TRUE, pattern = "_peaksChr6TAD58.tab", full.names = T)

cond1_DT <- foreach(i_file = cond1_files, .combine='rbind') %dopar% {
  curr_dt <- read.delim(i_file, header=F, col.names=c("name", "size", "covered", "sum", "mean0", "mean"), stringsAsFactors = F)
  curr_dt$sample <- gsub(paste0("(.+)_", cond1, "_.+_peaksChr6TAD58.tab"), "\\1", basename(i_file))
  curr_dt$mark <- gsub(paste0(".+_", cond1, "_(.+)_peaksChr6TAD58.tab"), "\\1", basename(i_file))
  curr_dt$status <- cond1
  curr_dt
}

cond2_DT <- foreach(i_file = cond2_files, .combine='rbind') %dopar% {
  curr_dt <- read.delim(i_file, header=F, col.names=c("name", "size", "covered", "sum", "mean0", "mean"), stringsAsFactors = F)
  curr_dt$sample <- gsub(paste0("(.+)_", cond2, "_.+_peaksChr6TAD58.tab"), "\\1", basename(i_file))
  curr_dt$mark <- gsub(paste0(".+_", cond2, "_(.+)_peaksChr6TAD58.tab"), "\\1", basename(i_file))
  curr_dt$status <- cond2
  curr_dt
}

# name - name field from bed, which should be unique
# size - size of bed (sum of exon sizes)
# covered - # bases within exons covered by bigWig
# sum - sum of values over all bases covered
# mean0 - average over bases with non-covered bases counting as zeroes
# mean - average over just covered bases

all_DT <- rbind(cond1_DT, cond2_DT)

top_tads <- c(unique(all_DT$name), "chr6_TAD58_allPeaks")

for(tad in top_tads) {
  
  if(tad == "chr6_TAD58_allPeaks") {
	  curr_df <- all_DT
    tadTit <- "chr6_TAD58 - circled peaks"
  } else {
    curr_df <- all_DT[all_DT$name == tad,]
    tadTit <- tad
  }

  
  for(yvar in c("sum", "mean0", "mean")) {
    
    # p <- ggboxplot(curr_df, 
    #           x = "status",
    #           xlab = paste0("CRC status"), 
    #           y = yvar,
    #           ylab = paste0("Signal ", yvar),
    #           legend.title="",
    #           legend="right",
    #           title = tad,
    #           # fill = "mark", palette = c("#00AFBB", "#E7B800"),
    #           fill = "mark", palette = c("steelblue3", "tan2"),
    #           add = "jitter")
    # p <- p + font("legend.text", size = 14)
    # p <- p + font("title", size = 18, face="bold")
    # p <- p + font("xlab", size = 16, face="bold")
    # p <- p + font("ylab", size = 16)
    # p <- p + font("xy.text", size = 14)
    # p <- p + theme(plot.title = element_text(hjust=0.5))
    # if(SSHFS) p

	### UPDATE: faceted version
    tmpDT <- aggregate(sample ~ status + mark, FUN=length, data = curr_df)
    subTit <- paste0(apply(tmpDT, 1, function(x) paste0("# ", x[1], " - ", x[2], " = ", x[3], collapse=" ")), collapse="; ")
    
    p <- ggboxplot(curr_df, x = "status", y = yvar,
                   fill = "status", palette = c("steelblue3", "tan2"),
                   add = "jitter",
                   subtitle = subTit, 
                   title = tadTit,
                   legend="none",
                   facet.by = "mark", short.panel.labs = T)
    # Use only p.format as label. Remove method name.
    p <- p + stat_compare_means(
      aes(label = paste0("p = ", ..p.format..)), label.x=0.5)
    p <- p + font("legend.text", size = 14)
    p <- p + font("title", size = 18, face="bold")
    p <- p + font("subtitle", face="bold")
    p <- p + font("xlab", size = 16, face="bold")
    p <- p + font("ylab", size = 16)
    p <- p + font("xy.text", size = 14)
    p <- p + theme(plot.title = element_text(hjust=0.5))
    if(SSHFS) p
    outFile <- file.path(outFold, paste0(tad, "_", yvar, "_histMarks.", plotType))
    ggsave(plot=p, filename = outFile, height=myHeight, width=myWidth )
    cat(paste0("... written: ", outFile, "\n"))
  }
}
  
