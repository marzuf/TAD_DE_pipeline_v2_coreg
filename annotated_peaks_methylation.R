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

#######################################################################################################################################################

# Rscript annotated_peaks_methylation.R
# everything hard-coded

library(foreach)

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS)  setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")

source("coreg_utils.R")

cat(paste0("> START ", "annotated_peaks_methylation.R",  "\n"))

inFold <- "PREP_CRC_METHYLATION"
outFold <- "ANNOTATED_PEAKS_METHYLATION"
system(paste0("mkdir -p ", outFold))


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 300, 7)
myWidth <- myHeight

# load probe coordinates
cat("... load probe coordinates DT\n")
probeCoord_DT <- read.delim(file.path(setDir, "/mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/hm450_probe_coordinates.txt"), header=T, stringsAsFactors = F)
stopifnot(!grepl("chr", probeCoord_DT$Chromosome))
probeCoord_DT$Chromosome <- paste0("chr", probeCoord_DT$Chromosome)
nrow(probeCoord_DT)
# 485577
stopifnot(is.numeric(probeCoord_DT$Genomic_Coordinate[1]))

logFile <- file.path(outFold, "annotated_peaks_crc_methylation.txt")
system(paste0("rm -f ", logFile))

methDT <- eval(parse(text = load(file.path(inFold, "COAD_methylation_DT.Rdata" ))))

mss_ID <- eval(parse(text = load(file.path(inFold, "mss_id.Rdata" ))))

msi_ID <- eval(parse(text = load(file.path(inFold, "msi_id.Rdata" ))))

all(colnames(methDT) %in% c(msi_ID,mss_ID))

probeMatchTol <- 10000
args <- commandArgs(trailingOnly = TRUE)
probeMatchTol <- as.numeric(args[1])
stopifnot(!is.na(probeMatchTol))

txt <- paste0("!!! HARD-CODED !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... probeMatchTol =\t", probeMatchTol, "\n")
printAndLog(txt, logFile)

mss_id <- mss_ID
mss_ID <- mss_ID[mss_ID %in% colnames(methDT)]
msi_id <- msi_ID
msi_ID <- msi_ID[msi_ID %in% colnames(methDT)]

# msi_id <- msi_ID[msi_ID %in% annotDT$sample_id]
txt <- paste0("... found msi_id with methylation data:\t", length(msi_ID), "/", length(msi_id), "\n")
printAndLog(txt, logFile)

# mss_id <- mss_ID[mss_ID %in% annotDT$sample_id]
txt <- paste0("... found mss_id with methylation data:\t", length(mss_ID), "/", length(mss_id), "\n")
printAndLog(txt, logFile)


# load info about the peaks
####################################################

# HARD CODED MAIN FOLDER
mainFold <- file.path(setDir, 
                      "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017")
h3k27ac_refPeakFile <- file.path(mainFold, "H3K27ac", paste0("H3K27ac", "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(h3k27ac_refPeakFile))
h3k27ac_refPeak_DT <- read.delim(h3k27ac_refPeakFile, header=F, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"), stringsAsFactors = FALSE)
stopifnot(is.numeric(h3k27ac_refPeak_DT$peak_start))
stopifnot(is.numeric(h3k27ac_refPeak_DT$peak_end))

h3k4me1_refPeakFile <- file.path(mainFold, "H3K4me1", paste0("H3K4me1", "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(h3k4me1_refPeakFile))
h3k4me1_refPeak_DT <- read.delim(h3k4me1_refPeakFile, header=F, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"), stringsAsFactors = FALSE)
stopifnot(is.numeric(h3k4me1_refPeak_DT$peak_start))
stopifnot(is.numeric(h3k4me1_refPeak_DT$peak_end))


### identified manually with threshold >= 2 ratioBinFC, annotated as regulatory regions for genes in curren TAD
all_peaks <- list(
  c(TAD_name = "chr6_TAD58", hist_mark = "H3K27ac", peak_name = "peak132954"),
  c(TAD_name = "chr6_TAD58", hist_mark = "H3K27ac", peak_name = "peak132947"),
  c(TAD_name = "chr6_TAD58", hist_mark = "H3K27ac", peak_name = "peak132953"),
  c(TAD_name = "chr6_TAD58", hist_mark = "H3K4me1", peak_name = "peak202346"),
  c(TAD_name = "chr6_TAD58", hist_mark = "H3K4me1", peak_name = "peak202349"),
  
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K27ac", peak_name = "peak35663"),
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K27ac", peak_name = "peak35647"),
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K27ac", peak_name = "peak35671"),
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K27ac", peak_name = "peak35655"),
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K4me1", peak_name = "peak52951"),  
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K4me1", peak_name = "peak52958"),  
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K4me1", peak_name = "peak52955"),  
  c(TAD_name = "chr12_TAD81", hist_mark = "H3K4me1", peak_name = "peak52952")
)  


txt <- paste0("> METHYLATION - ANNOTATED PEAKS\n")
printAndLog(txt, logFile)


peak_meth_DT <- foreach(i = 1:length(all_peaks), .combine='rbind') %do% {
  
  curr_tad <- all_peaks[[i]][["TAD_name"]]
  curr_peak <- all_peaks[[i]][["peak_name"]]
  curr_hist <- all_peaks[[i]][["hist_mark"]]

  txt <- paste0("... *** ", curr_tad, " - ", curr_hist, " - ", curr_peak, "\n")
  printAndLog(txt, logFile)
  
  # retrieve the position of the peak
  
  histDT <- eval(parse(text = paste0(tolower(curr_hist), "_refPeak_DT")))
  
  # retrieve the position
  curr_peak_start <- histDT$peak_start[histDT$peak_name == curr_peak]
  stopifnot(length(curr_peak_start) == 1)
  stopifnot(is.numeric(curr_peak_start))
  
  curr_peak_end <- histDT$peak_end[histDT$peak_name == curr_peak]
  stopifnot(length(curr_peak_end) == 1)
  stopifnot(is.numeric(curr_peak_end))
  
  curr_peak_chromo <- histDT$peak_chromo[histDT$peak_name == curr_peak]
  stopifnot(length(curr_peak_chromo) == 1)
  
  stopifnot(is.numeric(probeCoord_DT$Genomic_Coordinate))
  
  nestedProbe <- which(probeCoord_DT$Chromosome == curr_peak_chromo  &
                         probeCoord_DT$Genomic_Coordinate >= curr_peak_start & 
                         probeCoord_DT$Genomic_Coordinate <= curr_peak_end)
  
  tolLeftProbe <- which(probeCoord_DT$Chromosome == curr_peak_chromo  &
                           abs(probeCoord_DT$Genomic_Coordinate - curr_peak_start) <= probeMatchTol)
  
  tolRightProbe <- which(probeCoord_DT$Chromosome == curr_peak_chromo  &
                           abs(probeCoord_DT$Genomic_Coordinate - curr_peak_end) <= probeMatchTol)
  
  matchProbePeaks <- intersect(intersect(nestedProbe, tolLeftProbe), tolRightProbe)

  txt <- paste0("...... found # matching probes: ", length(matchProbePeaks), "\n")
  printAndLog(txt, logFile)
                         
  if(length(matchProbePeaks) > 0) {
    matchingProbes_id <- probeCoord_DT$Probe_id[matchProbePeaks]
    matchingProbes_gene <- probeCoord_DT$Gene_Symbol[matchProbePeaks]
    matchingProbes_chromo <- probeCoord_DT$Chromosome[matchProbePeaks]
    matchingProbes_coord <- probeCoord_DT$Genomic_Coordinate[matchProbePeaks]
    
    if(! any(matchingProbes_id %in% rownames(methDT)) ) {
      dt <- NULL
    } else {
      stopifnot(is.numeric(matchingProbes_coord))
      stopifnot(matchingProbes_chromo == curr_peak_chromo)
      stopifnot(any(rownames(methDT) %in% probeCoord_DT$Probe_id))
      
      stopifnot(any(rownames(methDT) %in% matchingProbes_id))
      
      subMeth_DT <- methDT[rownames(methDT) %in% matchingProbes_id,]
      
      mean_msi <- rowMeans(subMeth_DT[, msi_ID], na.rm=T)
      mean_mss <- rowMeans(subMeth_DT[, mss_ID], na.rm=T)
      stopifnot(names(mean_msi) == names(mean_mss) )
      av_probes <- names(mean_msi)
      av_probes_coord <- sapply(av_probes, function(x) probeCoord_DT$Genomic_Coordinate[probeCoord_DT$Probe_id == x]   )
      av_probes_gene <- sapply(av_probes, function(x) probeCoord_DT$Gene_Symbol[probeCoord_DT$Probe_id == x]   )
      
      dt <- data.frame(
        TAD = curr_tad,
        hist_mark = curr_hist,
        peak_name = curr_peak,
        peak_chromo = curr_peak_chromo,
        peak_start = curr_peak_start,
        peak_end = curr_peak_end,
        matchProbe_ID = av_probes,
        matchProbe_coordD = av_probes_coord,
        matchProbe_gene = av_probes_gene,
        meanMeth_MSI = as.numeric(mean_msi),
        meanMeth_MSS = as.numeric(mean_mss),
        stringsAsFactors = FALSE
      )
      dt$meanMeth_FC <- dt$meanMeth_MSI/dt$meanMeth_MSS
      dt <- dt[order(dt$hist_mark, dt$meanMeth_FC, decreasing=T),]
      
      dt_out <- dt
      dt_out$meanMeth_FC <- round(dt_out$meanMeth_FC, 4)
      dt_out$meanMeth_MSI <- round(dt_out$meanMeth_MSI, 4)
      dt_out$meanMeth_MSS <- round(dt_out$meanMeth_MSS, 4)
      
      outFile <- file.path(outFold, paste0(curr_tad, "_tol",probeMatchTol, "_peak_meth_DT.txt"))
      write.table(dt_out, file = outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
      cat(paste0("... written: ", outFile, "\n"))
    }
    
  } else {
    dt <- NULL
    
  }
  dt
}

if(!is.null(peak_meth_DT) > 0) {
  
  
  
  peak_meth_DT <- peak_meth_DT[order(peak_meth_DT$TAD, peak_meth_DT$hist_mark, peak_meth_DT$meanMeth_FC, decreasing = T),]  
  
  peak_meth_DT$meanMeth_FC <- round(peak_meth_DT$meanMeth_FC, 4)
  peak_meth_DT$meanMeth_MSI <- round(peak_meth_DT$meanMeth_MSI, 4)
  peak_meth_DT$meanMeth_MSS <- round(peak_meth_DT$meanMeth_MSS, 4)
  
  outFile <- file.path(outFold, paste0("peak_meth_DT_tol", probeMatchTol, ".Rdata"))
  save(peak_meth_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFold, paste0("peak_meth_DT_tol", probeMatchTol, ".txt"))
  write.table(peak_meth_DT, file = outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
  cat(paste0("... written: ", outFile, "\n"))
  
}



cat("\n*** DONE\n")
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)

