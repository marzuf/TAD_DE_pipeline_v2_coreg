# Rscript interaction_in_hic.R

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source("coreg_utils.R")

cat(paste0("> START ", "find_match_tissue_TAD.R",  "\n"))

inFold <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom/input_TopDom")
inFold <- file.path(setDir, "/mnt/nas_marie/TAD_DA_pipeline/TAD_call_pipeline_TopDom/input_TopDom")

tissue <- "SB"
binSize <- 40000
binSizeKb <- binSize/1000

nTop <- 5

outFold <- "INTERACTION_IN_HIC"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0(tissue, "_HiC_TAD_interactions.txt"))
system(paste0("rm -f ", logFile))

txt <- paste0("!!! HARD-CODED !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... binSize\t=\t", binSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... tissue\t=\t", tissue,"\n")
printAndLog(txt, logFile)
txt <- paste0("... nTop interactions\t=\t", nTop,"\n")
printAndLog(txt, logFile)

all_TADs <- list(
  c(TAD_name = "chr1_TAD150", chromo = "chr1", start = 89440001, end = 89920000),
  c(TAD_name = "chr12_TAD81", chromo = "chr12", start = 54160001, end = 54600000),
  c(TAD_name = "chr6_TAD58", chromo = "chr6", start = 32520001, end = 32840000)
)  

for(i in 1:length(all_TADs)) {
  
  curr_tad <- all_TADs[[i]][["TAD_name"]]
  curr_chromo <- all_TADs[[i]][["chromo"]]
  curr_start <- as.numeric(all_TADs[[i]][["start"]])
  curr_end <- as.numeric(all_TADs[[i]][["end"]])
  
  txt <- paste0("\n*** CONSENSUS TAD:\t", curr_tad, "\t", curr_start, " - ", curr_end, " ***\n")
  printAndLog(txt, logFile)
  
  tissueFile <- file.path(inFold, paste0(tissue, "_", curr_chromo, "_", binSizeKb, "k_matrix_pos_zero.txt"))
#  cat(paste0("tissueFile = ", tissueFile, "\n"))
  stopifnot(file.exists(tissueFile))  
  
  cat("... load Hi-C data\n")
  hic_DT <- read.delim(tissueFile, header=F, stringsAsFactors = FALSE)
  hic_DT_s <- hic_DT
  hic_DT[1:3,1:10]
  stopifnot( ncol(hic_DT) == nrow(hic_DT) + 3)
  
  hic_DT <- hic_DT[,-c(1:3)]
  stopifnot( ncol(hic_DT) == nrow(hic_DT) )
  
  # index of the start (1 -> 1, 40001 -> 2)
  startBin <- (curr_start-1)/binSize + 1
  # index of the start (40000 -> 1, 80000 -> 2)
  endBin <- curr_end/binSize

  stopifnot(endBin >= startBin)  
  
  sub_hic_DT <- hic_DT[startBin:endBin, startBin:endBin]
  cat("... subset Hi-C matrix\n")
  dim(sub_hic_DT)  
  
  # replace the diago with 0 
  diag(sub_hic_DT) <- 0
  
  j <- 0
  
  while(j < nTop) {
    
    ### MAX ONLY:
    max_interaction_row <- which.max(as.numeric(unlist(sub_hic_DT)))%/%ncol(sub_hic_DT) + 1
    max_interaction_col <- which.max(as.numeric(unlist(sub_hic_DT)))%%ncol(sub_hic_DT)
    if(max_interaction_col == 0) {
      max_interaction_col <- ncol(sub_hic_DT)
      max_interaction_row <- max_interaction_row - 1
    }
    stopifnot( sub_hic_DT[max_interaction_row, max_interaction_col] == max(sub_hic_DT) )
    stopifnot( sub_hic_DT[max_interaction_col, max_interaction_row] == max(sub_hic_DT) )
    
    binA_start <- curr_start + (max_interaction_row-1)*binSize
    binA_end <- binA_start + binSize - 1 
    binA_bin <- (binA_start-1)/binSize + 1
    
    binB_start <- curr_start + (max_interaction_col-1)*binSize
    binB_end <- binB_start + binSize - 1
    binB_bin <- (binB_start-1)/binSize + 1
    
    
    txt <- paste0("... ", curr_tad, " - interaction top ", j+1, "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... ", "interaction value = \t", max(sub_hic_DT), "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... ", "interaction binA = \t", binA_start , " - ", binA_end , "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("...... ", "interaction binB = \t", binB_start , " - ", binB_end , "\n")
    printAndLog(txt, logFile)
    
    
    stopifnot(hic_DT[binA_bin,binB_bin] ==  max(sub_hic_DT) )
    
    sub_hic_DT[max_interaction_row, max_interaction_col] <- 0
    sub_hic_DT[max_interaction_col, max_interaction_row] <- 0
    
    j <- j+1
    
  }
  
  
}
