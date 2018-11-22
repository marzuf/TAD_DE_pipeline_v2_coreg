options(scipen=100)

# args <- commandArgs(trailingOnly = TRUE)
# all_tads <- args

setDir <- "/media/electron"

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/zeroBasedTADcoord")
system(paste0("mkdir -p ", outFold))

tadcoordFile <- file.path(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")

tadcoordDT <- read.delim(tadcoordFile, header=F, stringsAsFactors=F, col.names=c("chromo", "tad", "start", "end"))

matchDT <- tadcoordDT

matchDT <- matchDT[grepl("_TAD", matchDT$tad),]

matchDT$start <- matchDT$start-1

# outFile <- file.path(outFold, paste0("bed_allTADs_zeroBased.bed"))
write.table(matchDT[,c("chromo", "start", "end", "tad")], col.names=F, row.names=F, quote=F, sep="\t", file = outFile)
cat(paste0("... written: ", outFile, "\n"))
