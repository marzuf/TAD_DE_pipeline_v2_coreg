
# Rscript find_match_tissue_TAD.R

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source("coreg_utils.R")

cat(paste0("> START ", "find_match_tissue_TAD.R",  "\n"))

inFold <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom/domains_TopDom")
#inFold <- file.path(setDir, "/mnt/nas_marie/TAD_DA_pipeline/TAD_call_pipeline_TopDom/output_TopDom")

tissue <- "SB"

outFold <- "MATCH_TISSUE_TADs"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0(tissue, "_matchingTADs.txt"))
system(paste0("rm -f ", logFile))


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

  tissueFile <- file.path(inFold, paste0(tissue, "_", curr_chromo, "_domains.bed"))
  stopifnot(file.exists(tissueFile))
  
  tissueTAD_dt <- read.delim(tissueFile, col.names=c("chromo", "start", "end"), stringsAsFactors = FALSE, header=F)
  
  tissue_midPos <- 0.5*(tissueTAD_dt$end + tissueTAD_dt$start)
  
  curr_midPos <- 0.5*(curr_start + curr_end)
  
  match_tad_idx <- which.min( abs(tissue_midPos - curr_midPos) )
  
  # tissueTAD_dt[match_tad_idx,]
  
  txt <- paste0("*** CONSENSUS TAD:\t", curr_tad, "\t", curr_start, " - ", curr_end, " ***\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... matching TAD in ", tissue, ":\t", tissueTAD_dt$chromo[match_tad_idx], "\t", tissueTAD_dt$start[match_tad_idx], " - ", tissueTAD_dt$end[match_tad_idx], "\n")
  printAndLog(txt, logFile)
  
  
  overlap_idxs <- which( (tissueTAD_dt$start <= curr_end & tissueTAD_dt$end >= curr_start) |
                           # or curr_tad is nested:
                           (curr_start > tissueTAD_dt$start & curr_start < tissueTAD_dt$end & curr_end > tissueTAD_dt$start & curr_end < tissueTAD_dt$end)
  )
  
  txt <- paste0("... overlapping TAD(s) in ", tissue, ":\n")
  printAndLog(txt, logFile)
  write.table(tissueTAD_dt[overlap_idxs,], file = logFile, col.names=F, row.names=F, sep="\t", quote=F, append=TRUE)
  write.table(tissueTAD_dt[overlap_idxs,], file = "", col.names=F, row.names=F, sep="\t", quote=F)
}
