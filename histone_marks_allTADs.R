library(foreach)
library(doMC)


SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 20, 2))

curr_dataset <- "TCGAcrc_msi_mss"

plotType <- "svg"
myHeight <- 7
myWidth <- 7

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/allTADs_average_GSE_data/HIST_MARKS_ALL_TADs")
system(paste0("mkdir -p ", outFold))

pvalCombFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/11_runEmpPvalCombined"), "emp_pval_combined.Rdata")
pval_comb <- eval(parse(text = load(pvalCombFile)))
pval_comb <- p.adjust(pval_comb, method="BH")

meanLogfcFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/3_runMeanTADLogFC"), "all_meanLogFC_TAD.Rdata")
meanLogFC <- eval(parse(text = load(meanLogfcFile)))

averageFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/allTADs_average_GSE_data")
averageFiles <- list.files(averageFolder, pattern="_allTADs.tab$", full.names = TRUE)

MSI_files <- averageFiles[grep("MSI_", basename(averageFiles))]
MSS_files <- averageFiles[grep("MSS_", basename(averageFiles))]

msi_DT <- foreach(i_file = MSI_files, .combine='rbind') %dopar% {
  curr_dt <- read.delim(i_file, header=F, col.names=c("name", "size", "covered", "sum", "mean0", "mean"), stringsAsFactors = F)
  curr_dt$sample <- gsub("MSI_(.+)_H3K.+", "\\1", basename(i_file))
  curr_dt$mark <- gsub("MSI_.+_.+_(H3K.+)_allTADs.tab", "\\1", basename(i_file))
  curr_dt$status <- "MSI"
  curr_dt
}

sum_mean_msi_DT <- aggregate(sum ~ name+mark, FUN=mean, data=msi_DT)
colnames(sum_mean_msi_DT)[colnames(sum_mean_msi_DT) == "sum"] <- "sum_MSI"

mean_mean_msi_DT <- aggregate(mean ~ name+mark, FUN=mean, data=msi_DT)
colnames(mean_mean_msi_DT)[colnames(mean_mean_msi_DT) == "mean"] <- "mean_MSI"

mean0_mean_msi_DT <- aggregate(mean0 ~ name+mark, FUN=mean, data=msi_DT)
colnames(mean0_mean_msi_DT)[colnames(mean0_mean_msi_DT) == "mean0"] <- "mean0_MSI"


mss_DT <- foreach(i_file = MSS_files, .combine='rbind') %dopar% {
  curr_dt <- read.delim(i_file, header=F, col.names=c("name", "size", "covered", "sum", "mean0", "mean"), stringsAsFactors = F)
  curr_dt$sample <- gsub("MSS_(.+)_H3K.+", "\\1", basename(i_file))
  curr_dt$mark <- gsub("MSS_.+_.+_(H3K.+)_allTADs.tab", "\\1", basename(i_file))
  curr_dt$status <- "MSS"
  curr_dt
}

# negative = more expressed in MSI

sum_mean_mss_DT <- aggregate(sum ~ name+mark, FUN=mean, data=mss_DT)
colnames(sum_mean_mss_DT)[colnames(sum_mean_mss_DT) == "sum"] <- "sum_MSS"

mean_mean_mss_DT <- aggregate(mean ~ name+mark, FUN=mean, data=mss_DT)
colnames(mean_mean_mss_DT)[colnames(mean_mean_mss_DT) == "mean"] <- "mean_MSS"

mean0_mean_mss_DT <- aggregate(mean0 ~ name+mark, FUN=mean, data=mss_DT)
colnames(mean0_mean_mss_DT)[colnames(mean0_mean_mss_DT) == "mean0"] <- "mean0_MSS"


sum_mean_DT <- merge(sum_mean_msi_DT, sum_mean_mss_DT, by =c("name", "mark"))
sum_mean_DT$sum_MSS_MSI_log2Ratio <- log2(sum_mean_DT$sum_MSS/sum_mean_DT$sum_MSI)
  
mean_mean_DT <- merge(mean_mean_msi_DT, mean_mean_mss_DT, by =c("name", "mark"))
mean_mean_DT$mean_MSS_MSI_log2Ratio <- log2(mean_mean_DT$mean_MSS/mean_mean_DT$mean_MSI)

mean0_mean_DT <- merge(mean0_mean_msi_DT, mean0_mean_mss_DT, by =c("name", "mark"))
mean0_mean_DT$mean0_MSS_MSI_log2Ratio <- log2(mean0_mean_DT$mean0_MSS/mean0_mean_DT$mean0_MSI)

marks_ratio_DT <- merge(merge(sum_mean_DT[,c("name", "mark", "sum_MSS_MSI_log2Ratio")],
                              mean_mean_DT[,c("name", "mark", "mean_MSS_MSI_log2Ratio")], by=c("name", "mark")), 
                        mean0_mean_DT[,c("name", "mark", "mean0_MSS_MSI_log2Ratio")], by=c("name", "mark"))

activatedMSI <- names(meanLogFC)[meanLogFC < 0]
activatedMSS <- names(meanLogFC)[meanLogFC > 0]

marks <- c("H3K27ac", "H3K4me1")
marksFun <- c("sum", "mean",  "mean0")

addLeg <- function(xvect, yvect, mypos="topright", corMeth="Pearson") {
  corTest <- cor.test(curr_x, curr_y)
  # legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%.2f", corTest$p.value), ")")
  legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%1.2e", corTest$p.value), ")")
  legend(mypos, legTxt, bty="n")
}

for(h_mark in marks) {
  for(h_fun in marksFun) {

      curr_histValues <- setNames(marks_ratio_DT[marks_ratio_DT$mark == h_mark, paste0(h_fun, "_MSS_MSI_log2Ratio")], 
                                  marks_ratio_DT$name[marks_ratio_DT$mark == h_mark])
      
      ### PLOT RELATIONSHIP WITH ADJ. COMBINED P-VAL FOR TADS ACTIVATED IN MSI
      # curr_x <- eval(parse(text = paste0(tolower(h_mark), "_", h_fun)))[activatedMSI] 
      curr_x <- curr_histValues[activatedMSI]
      curr_y <- -log10(pval_comb[activatedMSI])
      subTit <- paste0(curr_dataset, " - ", length(activatedMSI), " TADs")      
      outFile <- file.path(outFold, paste0(h_mark, "_", h_fun, "_", "pvalComb_MSIactiv.svg"))
      svg(outFile, height = myHeight, width = myWidth)
      plot(
        x = curr_x,
        xlab = paste0("log2 MSS/MSI - ", h_mark, ", ", h_fun, " signal"),
        y = curr_y,
        ylab = "-log10 p-val combined",
        main = "MSI activated TADs",
        pch=16,
        cex=0.7
      )
      mtext(text = subTit, side=3, font=3)
      addLeg(curr_x, curr_y)
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))

      ### PLOT RELATIONSHIP WITH ADJ. COMBINED P-VAL FOR TADS ACTIVATED IN MSI
      # curr_x <- eval(parse(text = paste0(tolower(h_mark), "_", h_fun)))[activatedMSS] 
      curr_x <- curr_histValues[activatedMSS]
      curr_y <- -log10(pval_comb[activatedMSS])
      subTit <- paste0(curr_dataset, " - ", length(activatedMSS), " TADs")
      outFile <- file.path(outFold, paste0(h_mark, "_", h_fun, "_", "pvalComb_MSSactiv.svg"))
      svg(outFile, height = myHeight, width = myWidth)
      plot(
        x = curr_x,
        xlab = paste0("log2 MSS/MSI - ", h_mark, ", ", h_fun, " signal"),
        y = curr_y,
        ylab = "-log10 p-val combined",
        main = "MSS activated TADs",
        pch=16,
        cex=0.7
      )
      mtext(text = subTit, side=3, font=3)
      addLeg(curr_x, curr_y)
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      ### PLOT RELATIONSHIP WITH LOG FC AND SIGNAL
      # curr_x <- eval(parse(text = paste0(tolower(h_mark), "_", h_fun)))[c(activatedMSI, activatedMSS)] 
      curr_x <- curr_histValues[c(activatedMSI, activatedMSS)]
      curr_y <- meanLogFC[c(activatedMSI, activatedMSS)] 
      subTit <- paste0(curr_dataset, " - ", length(c(activatedMSI, activatedMSS)), " TADs")
      outFile <- file.path(outFold, paste0(h_mark, "_", h_fun, "_", "logFC_MSI_MSS.svg"))
      svg(outFile, height = myHeight, width = myWidth)
      plot(
        x = curr_x,
        xlab = paste0("log2 MSS/MSI - ", h_mark, ", ", h_fun, " signal"),
        y = curr_y,
        ylab = "mean LogFC",
        main = "MSI + MSS activated TADs",
        pch=16,
        cex=0.7
      )
      mtext(text = subTit, side=3, font=3)
      addLeg(curr_x, curr_y)
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
  }
}

  
