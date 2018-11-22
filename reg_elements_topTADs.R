
source("coreg_utils.R")

# Rscript reg_elements_topTADs.R
cat(paste0("> START ", "reg_elements_topTADs.R",  "\n"))

# in this version assign mRNA to TAD based on entrezID

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"))

annotFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg", "reg_elements_data")

curr_dataset <- "TCGAcrc_msi_mss"
outFold <- file.path(setDir, 
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg",
                     "REG_ELEMENTS_TOPTADs",
                     curr_dataset)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "reg_elements_topTADs_logFile.txt")
system(paste0("rm -f ", logFile))


annot_datasets <- c(
  GeneHancer="genehancer.csv",
  vistaEnhancer="vistaEnhancers.txt",
  ensemblOtherRegRegions="grch37_ensembl_other_reg_regions.csv",
  # ensemblRegFeatures="grch37_ensembl_other_reg_regions.csv",
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

cat(paste0("... load ", "vistaEnhancer", " dataset\n"))
vistaEnhancerDT <- read.delim(annotFiles["vistaEnhancer"], header = FALSE, stringsAsFactors=F)
# 609     chr1    3190581 3191428 element_705     900
# 647     chr1    8130439 8131887 element_1833    900
colnames(vistaEnhancerDT)[2] <- "chromo"
colnames(vistaEnhancerDT)[3] <- "start"
colnames(vistaEnhancerDT)[4] <- "end"
stopifnot(is.numeric(vistaEnhancerDT$start))
stopifnot(is.numeric(vistaEnhancerDT$end))

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

cat(paste0("... load ", "ensemblRegFeatures", " dataset\n"))
# regFeatDT <- read.delim(annotFiles["ensemblRegFeatures"], header = F, stringsAsFactors=F)

cat(paste0("... load ", "ensemblRegEvidence", " dataset\n"))
# regEvidDT <- read.delim(annotFiles["ensemblRegEvidence"], header = F, stringsAsFactors=F)

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


all_annot_DT <- list(
  GeneHancer=genehancerDT,
  vistaEnhancer=vistaEnhancerDT,
  ensemblOtherRegRegions=otherRegDT,
  # ensemblRegFeatures=regFeatDT,
  # ensemblRegEvidence=regEvidDT,
  ensemblRegFeaturesGff=regFeaturesGffDT,
  ensemblRegMotifGff=regFeaturesGffDT
)

# all_annot_DT <- list(
#   GeneHancer=genehancerDT,
#   vistaEnhancer=vistaEnhancerDT,
#   ensemblOtherRegRegions=otherRegDT
# )


# TOP TADs
cat("... load topTADs data\n")
topTADsFile <- file.path(setDir,
                         "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data",
                         paste0(curr_dataset, "_topTADs10.bed"))
stopifnot(file.exists(topTADsFile))
topTADs_DT <- read.delim(topTADsFile, col.names = c("chromo", "start", "end", "region"), stringsAsFactors = F)
topTADs <- topTADs_DT$region
topTADs <- topTADs[topTADs %in% c("chr1_TAD150", "chr6_TAD58", "chr12_TAD81")]
stopifnot(length(topTADs) > 0)

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
  
  stopifnot(length(all_annot_DT) > 0)
  for(i in 1:length(all_annot_DT)) {
    
    curr_dt <- all_annot_DT[[i]]
    
    txt <- "\n"
    printAndLog(txt, logFile)
    txt <- paste0("... *** dataset: ", names(all_annot_DT)[i], "\n")
    printAndLog(txt, logFile)
    
    inTAD_DT <- curr_dt[curr_dt$chromo == curr_chromo &
                          curr_dt$start >= curr_start &
                          curr_dt$end <= curr_end,]
    
    write.table(inTAD_DT, file = logFile, sep="\t", quote=F, col.names=TRUE, row.names=F, append=T)
    txt <- "\n"
    printAndLog(txt, logFile)
    
    outFile <- file.path(outFold, paste0(curr_tad, "_", names(all_annot_DT)[i], ".txt"))
    write.table(inTAD_DT, file = outFile, sep="\t", quote=F, col.names=FALSE, row.names=F, append=F)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}






#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
cat(paste0("... written: ", logFile, "\n"))










