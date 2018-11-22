
source("coreg_utils.R")

library(GenomicRanges)

# Rscript regElements_histPeaks.R
cat(paste0("> START ", "regElements_histPeaks.R",  "\n"))

# in this version assign mRNA to TAD based on entrezID

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"))

annotFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg", "reg_elements_data")

curr_dataset <- "TCGAcrc_msi_mss"
outFold <- file.path(setDir, 
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg",
                     "REGELEMENTS_HISTPEAKS",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "regElements_histPeaks_logFile.txt")
system(paste0("rm -f ", logFile))

all_chromo <- paste0("chr", c(1:22,"X"))


########################################
########## LOAD REGULATORY ELEMENTS
########################################

annot_datasets <- c(
  GeneHancer="genehancer.csv",
  vistaEnhancer="vistaEnhancers.txt",
  ensemblOtherRegRegions="grch37_ensembl_other_reg_regions.csv",
  ensemblRegFeatures="grch37_ensembl_other_reg_regions.csv",
  # ensemblRegEvidence="grch37_ensembl_other_reg_regions.csv",
  ensemblRegFeaturesGff="homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff",
  ensemblRegMotifGff="homo_sapiens.GRCh37.motiffeatures.20161117.gff"
)

annotFiles <- setNames(file.path(annotFold, annot_datasets), names(annot_datasets))

cat(paste0("... load ", "GeneHancer", " dataset\n"))
genehancerDT <- read.delim(annotFiles["GeneHancer"], header = TRUE, sep="," , stringsAsFactors=F)
# chrom,source,feature name,start,end,score,strand,frame,attributes
# chr6,GeneHancer,Enhancer,133518045,133518161,0.22,.,.,genehancer_id=GH06I133518;connected_gene=EYA4;score=7.10;connected_gene=HSPE1P21;score=0.44;connected_gene=ENSG00000272428;score=0.34;connected_gene=TARID;score=0.05
colnames(genehancerDT)[colnames(genehancerDT) == "chrom"] <- "chromo"
stopifnot(is.numeric(genehancerDT$start))
stopifnot(is.numeric(genehancerDT$end))
genehancerDT <- genehancerDT[,c("chromo", "source", "feature.name", "start", "end", "attributes")]
genehancerDT <- genehancerDT[genehancerDT$chromo %in% all_chromo,]
stopifnot(nrow(genehancerDT) > 0)
genehancerDT$my_desc <- genehancerDT$attributes
genehancerDT$DS <- "GeneHancer"

cat(paste0("... load ", "vistaEnhancer", " dataset\n"))
vistaEnhancerDT <- read.delim(annotFiles["vistaEnhancer"], header = FALSE, stringsAsFactors=F)
# 609     chr1    3190581 3191428 element_705     900
# 647     chr1    8130439 8131887 element_1833    900
colnames(vistaEnhancerDT)[2] <- "chromo"
colnames(vistaEnhancerDT)[3] <- "start"
colnames(vistaEnhancerDT)[4] <- "end"
stopifnot(is.numeric(vistaEnhancerDT$start))
stopifnot(is.numeric(vistaEnhancerDT$end))
vistaEnhancerDT <- vistaEnhancerDT[vistaEnhancerDT$chromo %in% all_chromo,]
stopifnot(nrow(vistaEnhancerDT) > 0)
vistaEnhancerDT$my_desc <- paste0("vistaEnhancer_", vistaEnhancerDT$V5)
vistaEnhancerDT$DS <- "vistaEnhancer"

cat(paste0("... load ", "ensemblOtherRegRegions", " dataset\n"))
otherRegDT<- read.delim(annotFiles["ensemblOtherRegRegions"], header = TRUE, stringsAsFactors=F)
# Chromosome/scaffold Name        Start (bp)      End (bp)        Feature type    Feature type class      Feature type description
# 19      2634884 2634884 FANTOM predictions      Transcription Start Site        FANTOM TSS, relaxed
# 19      2634532 2634552 FANTOM predictions      Transcription Start Site        FANTOM TSS, relaxed
colnames(otherRegDT)[1] <- "chromo"
colnames(otherRegDT)[2] <- "start"
colnames(otherRegDT)[3] <- "end"
stopifnot(is.numeric(otherRegDT$start))
stopifnot(is.numeric(otherRegDT$end))
otherRegDT$chromo <- paste0("chr", otherRegDT$chromo)
otherRegDT <- otherRegDT[otherRegDT$chromo %in% all_chromo,]
stopifnot(nrow(otherRegDT) > 0)
otherRegDT$my_desc <- otherRegDT$Feature.type.description
otherRegDT$DS <- "ensemblOtherReg"

cat(paste0("... load ", "ensemblRegFeatures", " dataset\n"))
regFeatDT <- read.delim(annotFiles["ensemblRegFeatures"], header = T, stringsAsFactors=F)
colnames(regFeatDT)[1:5] <- c("chromo", "start", "end", "bound_start", "bound_end")
stopifnot(is.numeric(regFeatDT$start))
stopifnot(is.numeric(regFeatDT$end))
regFeatDT$chromo <- paste0("chr", regFeatDT$chromo)
regFeatDT <- regFeatDT[regFeatDT$chromo %in% all_chromo,]
stopifnot(nrow(regFeatDT) > 0)
regFeatDT$my_desc <- paste0(regFeatDT$Feature.type.description)
regFeatDT$DS <- "ensemblRegFeatures"
#   Chromosome/scaffold name,Start (bp),End (bp),Bound start (bp),Bound end (bp),Feature type,Feature type description,Activity
# 15,102118789,102119129,102118695,102119230,TF binding site,Transcription factor binding site,INACTIVE

cat(paste0("... load ", "ensemblRegEvidence", " dataset\n"))
# regEvidDT <- read.delim(annotFiles["ensemblRegEvidence"], header = F, stringsAsFactors=F)
# stopifnot(is.numeric(regEvidDT$start))
# stopifnot(is.numeric(regEvidDT$end))
# regEvidDT <- regEvidDT[regEvidDT$chromo %in% all_chromo,]
# stopifnot(nrow(regEvidDT) > 0)
# regEvidDT$my_desc <- 

cat(paste0("... load ", "ensemblRegFeaturesGff", " dataset\n"))
regFeaturesGffDT <- read.delim(annotFiles["ensemblRegFeaturesGff"], header = F, stringsAsFactors=F)
# 15      Regulatory_Build        regulatory_region       102118789       102119129       .       .       .       ID=ENSR00000368862;bound_end=102119230;bound_start=102118695;description=Transcription factor binding site;feature_type=TF binding site
# 7       Regulatory_Build        regulatory_region       1619070 1619637 .       .       .       ID=ENSR00000408425;bound_end=1619683;bound_start=1619063;description=Transcription factor binding site;feature_type=TF binding site
colnames(regFeaturesGffDT)[1] <- "chromo"
colnames(regFeaturesGffDT)[4] <- "start"
colnames(regFeaturesGffDT)[5] <- "end"
stopifnot(is.numeric(regFeaturesGffDT$start))
stopifnot(is.numeric(regFeaturesGffDT$end))
regFeaturesGffDT$chromo <- paste0("chr", regFeaturesGffDT$chromo)
regFeaturesGffDT$V6 <- NULL
regFeaturesGffDT$V7 <- NULL
regFeaturesGffDT$V8 <- NULL
regFeaturesGffDT <- regFeaturesGffDT[regFeaturesGffDT$chromo %in% all_chromo,]
stopifnot(nrow(regFeaturesGffDT) > 0)
regFeaturesGffDT$my_desc <- paste0(regFeaturesGffDT$V3, ";", regFeaturesGffDT$V9)
regFeaturesGffDT$DS <- "regFeaturesGFF"

cat(paste0("... load ", "ensemblRegMotifGff", " dataset\n"))
motifFeaturesGffDT <- read.delim(annotFiles["ensemblRegMotifGff"], header = F, stringsAsFactors=F)
# 1       .       TF_binding_site 11866228        11866234        9.25    +       .       binding_matrix=MA0409.1;motif_feature_type=SRebp1
# 5       .       TF_binding_site 168006637       168006643       9.25    +       .       binding_matrix=MA0409.1;motif_feature_type=SRebp1
colnames(motifFeaturesGffDT)[1] <- "chromo"
colnames(motifFeaturesGffDT)[4] <- "start"
colnames(motifFeaturesGffDT)[5] <- "end"
stopifnot(is.numeric(motifFeaturesGffDT$start))
stopifnot(is.numeric(motifFeaturesGffDT$end))
motifFeaturesGffDT$chromo <- paste0("chr", motifFeaturesGffDT$chromo)
motifFeaturesGffDT$V6 <- NULL
motifFeaturesGffDT$V7 <- NULL
motifFeaturesGffDT$V8 <- NULL
motifFeaturesGffDT <- motifFeaturesGffDT[motifFeaturesGffDT$chromo %in% all_chromo,]
stopifnot(nrow(motifFeaturesGffDT) > 0)
motifFeaturesGffDT$my_desc <- paste0(motifFeaturesGffDT$V3, ";", motifFeaturesGffDT$V9)
motifFeaturesGffDT$DS <- "motifFeaturesGFF"

all_annot_DT <- list(
  GeneHancer=genehancerDT[,c("DS", "chromo", "start","end", "my_desc")],
  vistaEnhancer=vistaEnhancerDT[,c("DS", "chromo", "start","end", "my_desc")],
  ensemblOtherRegRegions=otherRegDT[,c("DS", "chromo", "start","end", "my_desc")],
  ensemblRegFeatures=regFeatDT[,c("DS", "chromo", "start","end", "my_desc")],
  # ensemblRegEvidence=regEvidDT[,c("DS", "chromo", "start","end", "my_desc")],
  ensemblRegFeaturesGff=regFeaturesGffDT[,c("DS", "chromo", "start","end", "my_desc")],
  ensemblRegMotifGff=regFeaturesGffDT[,c("DS", "chromo", "start","end", "my_desc")]
)

regFeatures_DT <- do.call(rbind, all_annot_DT)
rownames(regFeatures_DT) <- NULL
colnames(regFeatures_DT) <- c("DS","feature_chromo","feature_start","feature_end","my_desc")
regFeatures_DT <- unique(regFeatures_DT)

########################################
########## LOAD TOP TADs
########################################
cat("... load topTADs data\n")
topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region
topTADs <- topTADs[topTADs %in% c("chr1_TAD150", "chr6_TAD58", "chr12_TAD81")]
stopifnot(length(topTADs) > 0)

########################################
########## LOAD REFERENCE PEAKS
########################################

cat("... prepare histone peaks data\n")

mainFold <- file.path(setDir, 
                      "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017")
h3k27ac_DTfile <- file.path(mainFold, "H3K27ac", paste0("H3K27ac", "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(h3k27ac_DTfile))

h3k4me1_DTfile <- file.path(mainFold, "H3K4me1", paste0("H3K4me1", "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(h3k4me1_DTfile))

h3k27acPeaksDT <- read.delim(h3k27ac_DTfile, stringsAsFactors = FALSE, header=F, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"))
head(h3k27acPeaksDT)
h3k27acPeaksDT$histMark <- "H3K27ac"
h3k4me1PeaksDT <- read.delim(h3k4me1_DTfile, stringsAsFactors = FALSE, header=F, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"))
head(h3k4me1PeaksDT)
h3k4me1PeaksDT$histMark <- "H3K4me1"
hist_DT <- rbind(h3k27acPeaksDT,h3k4me1PeaksDT)


# HARD-CODED: FOCUS ON THE FOLLOWING TADs

txt <- paste0("!!! HARD-CODED !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... focus on following TADs:\t", paste0(topTADs, collapse = ", "), "\n")
printAndLog(txt, logFile)

txt <- paste0("... use following annotation data:\t", paste0(annot_datasets, collapse = ", "), "\n")
printAndLog(txt, logFile)

topTADs <- c("chr1_TAD150", "chr6_TAD58", "chr12_TAD81")

curr_tad="chr1_TAD150"

for(curr_tad in topTADs) {
  
  txt <- paste0(">>>", curr_tad, "\n")
  printAndLog(txt, logFile)
  
  tad_idx <- which(topTADs_DT$region == curr_tad)
  stopifnot(length(tad_idx) == 1)
  
  curr_chromo <- topTADs_DT$chromo[tad_idx]
  curr_start <- topTADs_DT$start[tad_idx]
  curr_end <- topTADs_DT$end[tad_idx]
  
  stopifnot(is.numeric(curr_start))
  stopifnot(is.numeric(curr_end))
  
  subHist_DT <- hist_DT[hist_DT$peak_chromo == curr_chromo &
                          hist_DT$peak_start >= curr_start &
                          hist_DT$peak_end <= curr_end,]
  subHist_DT <- subHist_DT[order(subHist_DT$peak_chromo, subHist_DT$peak_start, subHist_DT$peak_end ),] 
  
  sub_histPeaks_GR <- GRanges(seqnames = subHist_DT$peak_chromo,
                          ranges = IRanges(start=subHist_DT$peak_start, end = subHist_DT$peak_end))
  
  subFeatures_DT <- regFeatures_DT[regFeatures_DT$feature_chromo == curr_chromo &
                          regFeatures_DT$feature_start >= curr_start &
                          regFeatures_DT$feature_end <= curr_end,]
  subFeatures_DT <- subFeatures_DT[order(subFeatures_DT$feature_chromo, subFeatures_DT$feature_start, subFeatures_DT$feature_end),]
  
  sub_regFeatures_GR <- GRanges(seqnames = subFeatures_DT$feature_chromo,
                            ranges = IRanges(start=subFeatures_DT$feature_start, end = subFeatures_DT$feature_end))
  
  # findOverlaps uses the triplet (sequence name, range, strand) 
  # to determine which features FROM THE QUERY overlap which features IN THE SUBJECT
  overlap_histPeaks_regFeatures_DT <- findOverlaps(query = sub_histPeaks_GR,
               subject = sub_regFeatures_GR)
  
  histPeaksOverlap_DT <- subHist_DT[overlap_histPeaks_regFeatures_DT@from,]
  regFeaturesOverlap_DT <- subFeatures_DT[overlap_histPeaks_regFeatures_DT@to,]
  
  overlap_DT <- cbind(histPeaksOverlap_DT, regFeaturesOverlap_DT)
  rownames(overlap_DT) <- NULL
    
  stopifnot(pmin(overlap_DT$feature_end, overlap_DT$peak_end) >= pmax(overlap_DT$feature_start, overlap_DT$peak_start) )
  
  outFile <- file.path(outFold, paste0(curr_tad, "_overlap_histPeaks_regFeatures.txt"))
  write.table(overlap_DT, file = outFile, sep="\t", quote=F, col.names=TRUE, row.names=F, append=F)
  cat(paste0("... written: ", outFile, "\n"))
  
  summary_overlap_DT <- overlap_DT[,c("histMark", "peak_chromo"	,"peak_start",	"peak_end"	,"peak_name")]
  summary_overlap_DT <- unique(summary_overlap_DT)
  summary_overlap_DT <- summary_overlap_DT[order(summary_overlap_DT$histMark,
                                                 summary_overlap_DT$peak_chromo,
                                                 summary_overlap_DT$peak_start,
                                                 summary_overlap_DT$peak_end),]
  outFile <- file.path(outFold, paste0(curr_tad, "_annotated_histPeaks_regFeatures_summaryDT.txt"))
  write.table(summary_overlap_DT, file = outFile, sep="\t", quote=F, col.names=TRUE, row.names=F, append=F)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}






#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))










