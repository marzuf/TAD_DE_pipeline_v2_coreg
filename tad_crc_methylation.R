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

# Rscript tad_crc_methylation.R
# everything hard-coded

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source("coreg_utils.R")

cat(paste0("> START ", "tad_crc_methylation.R",  "\n"))

inFold <- "PREP_CRC_METHYLATION"
outFold <- "TAD_CRC_METHYLATION"
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

logFile <- file.path(outFold, "tad_crc_methylation.txt")
system(paste0("rm -f ", logFile))

methDT <- eval(parse(text = load(file.path(inFold, "COAD_methylation_DT.Rdata" ))))

mss_ID <- eval(parse(text = load(file.path(inFold, "mss_id.Rdata" ))))

msi_ID <- eval(parse(text = load(file.path(inFold, "msi_id.Rdata" ))))

all(colnames(methDT) %in% c(msi_ID,mss_ID))

# txt <- paste0("!!! HARD-CODED !!!\n")
# printAndLog(txt, logFile)
# txt <- paste0("... sampType =\t", sampType, "\n")
# printAndLog(txt, logFile)

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


# chr1_TAD150	chr1	89440001	89920000
# chr12_TAD81	chr12	54160001	54600000
# chr6_TAD58	chr6	32520001	32840000
  
all_TADs <- list(
  c(TAD_name = "chr1_TAD150", chromo = "chr1", start = 89440001, end = 89920000),
  c(TAD_name = "chr12_TAD81", chromo = "chr12", start = 54160001, end = 54600000),
  c(TAD_name = "chr6_TAD58", chromo = "chr6", start = 32520001, end = 32840000)
)  

txt <- paste0("> METHYLATION - WHOLE TAD\n")
printAndLog(txt, logFile)


for(i in 1:length(all_TADs)) {
  
  curr_tad <- all_TADs[[i]][["TAD_name"]]
  curr_chromo <- all_TADs[[i]][["chromo"]]
  curr_start <- as.numeric(all_TADs[[i]][["start"]])
  curr_end <- as.numeric(all_TADs[[i]][["end"]])
  
  txt <- paste0("*** ", curr_tad, "\t", curr_start, " - ", curr_end, " ***\n")
  printAndLog(txt, logFile)
  
  # find if any probe is matching
  sub_probes <- probeCoord_DT[probeCoord_DT$Chromosome == curr_chromo  &
    probeCoord_DT$Genomic_Coordinate >= curr_start & 
      probeCoord_DT$Genomic_Coordinate <= curr_end, 
      ]
  tad_probesID <- as.character(sub_probes$Probe_id)
  stopifnot(any(tad_probesID %in% rownames(methDT)))
  tad_probes_id <- tad_probesID
  tad_probesID <- tad_probesID[tad_probesID %in% rownames(methDT)]
  
  txt <- paste0("... found # of probes in TAD:\t", length(tad_probesID), "\n")  
  printAndLog(txt, logFile)
  
  methDT_msi <- methDT[tad_probesID, msi_ID]
  methDT_mss <- methDT[tad_probesID, mss_ID]
  
  msi_values <- as.numeric(unlist(methDT_msi))
  mss_values <- as.numeric(unlist(methDT_mss))
  
  meanMSI <- mean(msi_values, na.rm=T)
  meanMSS <- mean(mss_values, na.rm=T)
  
  txt <- paste0("... mean methylation MSI:\t", round(meanMSI, 4), "\n")  
  printAndLog(txt, logFile)
  
  txt <- paste0("... mean methylation MSS:\t", round(meanMSS, 4), "\n")  
  printAndLog(txt, logFile)
  
  tmpDT <- data.frame(status = c(rep("MSI", length(msi_values)), 
                                     rep("MSS", length(mss_values))) ,
                      values =c(msi_values, mss_values), stringsAsFactors = FALSE )
  

  outFile <- file.path(outFold, paste0(curr_tad, "_cmp_meth.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))                     
  boxplot(values ~ status, data =tmpDT,
          main = curr_tad)
  foo <- dev.off()
  cat(paste0("... written: ", outFile,  "\n"))
  
  wt <- wilcox.test(x = msi_values,
              y = mss_values)
  
  txt <- paste0("... two-sided Wilcoxon's test MSI vs. MSS:\t", sprintf("%2.2e", wt$p.value) , "\n")   #round(wt$p.value, 4)
  printAndLog(txt, logFile)
  
}


txt <- paste0("\n> METHYLATION - HISTONE PEAKS IN TAD\n")
printAndLog(txt, logFile)


# load histone peaks to investigate methylation in peak regions

all_hist_marks <- c("H3K27ac", "H3K4me1")

mainFold <- file.path(setDir, 
                      "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017")
histMark <- "H3K27ac"
for(histMark in all_hist_marks) {
  refPeakFile <- file.path(mainFold, histMark, paste0(histMark, "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
  stopifnot(file.exists(refPeakFile))
  refPeaks_DT <- read.delim(refPeakFile, stringsAsFactors = FALSE, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"), header=FALSE)
  # check the peaks are not overlapping
  peak_ends <- refPeaks_DT$end[-length(refPeaks_DT$end)]
  peak_starts <- refPeaks_DT$start[-1]
  stopifnot(peak_ends <= peak_starts)
  
  
  for(i in 1:length(all_TADs)) {
    
    curr_tad <- all_TADs[[i]][["TAD_name"]]
    curr_chromo <- all_TADs[[i]][["chromo"]]
    curr_start <- as.numeric(all_TADs[[i]][["start"]])
    curr_end <- as.numeric(all_TADs[[i]][["end"]])
    
    txt <- paste0("*** ", histMark, ":\t ", curr_tad, "\t", curr_start, " - ", curr_end, " ***\n")
    printAndLog(txt, logFile)
    
    
    chromo_refPeaks_DT <- refPeaks_DT[refPeaks_DT$peak_chromo == curr_chromo,]
    stopifnot(nrow(chromo_refPeaks_DT) > 0)
    # peaks overlapping with the TADs
    
    overlap_idxs <- which( 
      (chromo_refPeaks_DT$peak_start <= curr_end & chromo_refPeaks_DT$peak_end >= curr_start) |
                             # or curr_tad is nested:
                             (curr_start > chromo_refPeaks_DT$peak_start & curr_start < chromo_refPeaks_DT$peak_end & curr_end > chromo_refPeaks_DT$peak_start & curr_end < chromo_refPeaks_DT$peak_end)
    )
    
    tadPeaks_DT <- chromo_refPeaks_DT[overlap_idxs,]
    
    txt <- paste0("... # of peaks in TAD:\t", nrow(tadPeaks_DT), "\n")
    printAndLog(txt, logFile)    

    if(nrow(tadPeaks_DT) == 0) next
        
    # find probes that overlap with peaks
    
    matchProbes <- sapply(1:nrow(tadPeaks_DT), function(x) {
      peak_chromo <- tadPeaks_DT$peak_chromo[x]
      peak_start <- tadPeaks_DT$peak_start[x]
      peak_end <- tadPeaks_DT$peak_end[x]
      stopifnot(peak_chromo == curr_chromo)
      stopifnot(as.character(peak_chromo) == as.character(curr_chromo))
      stopifnot(is.numeric(peak_start))
      stopifnot(is.numeric(peak_end))
      stopifnot(is.numeric(probeCoord_DT$Genomic_Coordinate))
      # find peaks that overlap with probes
      which(probeCoord_DT$Chromosome == curr_chromo  &
              probeCoord_DT$Genomic_Coordinate >= peak_start & 
              probeCoord_DT$Genomic_Coordinate <= peak_end)
      
      
    })
    matchProbes <- unique(unlist(matchProbes))
    
    txt <- paste0("... # of probes with coordinates in peaks:\t", length(matchProbes), "\n")
    printAndLog(txt, logFile)    
    
    if(length(matchProbes) == 0) next
    
    matchProbes_DT <- probeCoord_DT[matchProbes,]
    
    stopifnot(is.numeric(curr_start))
    stopifnot(is.numeric(curr_end))
    stopifnot(is.numeric(matchProbes_DT$Genomic_Coordinate))

    stopifnot(matchProbes_DT$Genomic_Coordinate >= curr_start & matchProbes_DT$Genomic_Coordinate <= curr_end )
    stopifnot(matchProbes_DT$Chromosome == curr_chromo)
    
    tad_probesID <- as.character(matchProbes_DT$Probe_id)
    stopifnot(any(tad_probesID %in% rownames(methDT)))
    
    tad_probesID <- tad_probesID[tad_probesID %in% rownames(methDT)]
    txt <- paste0("... found # of probes in ", histMark, " peaks in TAD:\t", length(tad_probesID), "\n")  
    printAndLog(txt, logFile)
    
    methDT_msi <- methDT[tad_probesID, msi_ID]
    methDT_mss <- methDT[tad_probesID, mss_ID]
    
    msi_values <- as.numeric(unlist(methDT_msi))
    mss_values <- as.numeric(unlist(methDT_mss))
    
    meanMSI <- mean(msi_values, na.rm=T)
    meanMSS <- mean(mss_values, na.rm=T)
    
    txt <- paste0("... mean methylation MSI:\t", round(meanMSI, 4), "\n")  
    printAndLog(txt, logFile)
    
    txt <- paste0("... mean methylation MSS:\t", round(meanMSS, 4), "\n")  
    printAndLog(txt, logFile)
    
    tmpDT <- data.frame(status = c(rep("MSI", length(msi_values)), 
                                   rep("MSS", length(mss_values))) ,
                        values =c(msi_values, mss_values), stringsAsFactors = FALSE )
    
    outFile <- file.path(outFold, paste0(histMark, "peaks_", curr_tad, "_cmp_meth.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))                     
    boxplot(values ~ status, data =tmpDT,
            main = paste0("Methylation in ", histMark, " peaks - ", curr_tad))
    foo <- dev.off()
    cat(paste0("... written: ", outFile,  "\n"))
    
    wt <- wilcox.test(x = msi_values,
                      y = mss_values)
    
    txt <- paste0("... two-sided Wilcoxon's test MSI vs. MSS:\t", sprintf("%2.2e", wt$p.value) , "\n")   #round(wt$p.value, 4)
    printAndLog(txt, logFile)
    
    
    
    
  }
  
  
  
  
}











# # msi_id <- annotDT$sample_id[annotDT$subtype == msi_subtype & annotDT$type == sampType]
# # mss_id <- annotDT$sample_id[annotDT$subtype == mss_subtype & annotDT$type == sampType]
# save(msi_id, file = file.path(outFold, "msi_id.Rdata"))
# save(mss_id, file = file.path(outFold, "mss_id.Rdata"))

# if(file.exists(file.path(outFold, "msi_id.Rdata"))){
#   cat(paste0("...written: ", file.path(outFold, "msi_id.Rdata"), "\n"))
# }else{
#   stop("error\n")
# }
# 
# if(file.exists(file.path(outFold, "mss_id.Rdata"))){
#   cat(paste0("...written: ", file.path(outFold, "mss_id.Rdata"), "\n"))
# }else{
#   stop("error\n")
# }

#cat("... load TCGA rnaseq \t")
#TCGA_rnaseqFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/gene_expression/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv")
#TCGA_DT <- read.delim(TCGA_rnaseqFile, header=T, stringsAsFactors=F)
#cat("done\n")
#colnames(TCGA_DT) <- substr(x = colnames(TCGA_DT), start = 1, stop = 15)
#colnames(TCGA_DT) <- gsub("\\.", "-", colnames(TCGA_DT))
#stopifnot(any(c(msi_id, mss_id) %in% colnames(TCGA_DT)))
#inSTAD_id = colnames(TCGA_DT)[colnames(TCGA_DT) %in% c(msi_id, mss_id)]
#STAD_RNAseq_DT <- TCGA_DT[,c("gene_id", inSTAD_id)]
#save(STAD_RNAseq_DT, file = file.path(outFold, "STAD_RNAseq_DT.Rdata"))
#if(file.exists(file.path(outFold, "STAD_RNAseq_DT.Rdata"))){
#cat(paste0("...written: ", file.path(outFold, "STAD_RNAseq_DT.Rdata"), "\n"))
#}else{
#stop("error\n")
#}



cat("\n*** DONE\n")
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)


# methylProbeAnnotFile <- file.path(setDir, "/mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/hm450_probe_coordinates.txt")
# methylGeneAnnotFile <- file.path(setDir, "/mnt/pd2/marco/CNS_PNET/CNS_PNET_v2/data/annotation/Methlation_450k/HM450K_probe_FANTOM5_promoter_annotation_table.csv")


