args <- commandArgs(trailingOnly = TRUE)

all_tads <- args

tadcoordFile <- "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt"

tadcoordDT <- read.delim(tadcoordFile, header=F, stringsAsFactors=F, col.names=c("chromo", "tad", "start", "end"))


matchDT <- tadcoordDT[tadcoordDT$tad %in% all_tads,]

matchDT <- matchDT[match(all_tads, matchDT$tad),]


write.table(matchDT, col.names=F, row.names=F, quote=F, sep="\t")

cat("\n")

write.table(matchDT[,c("chromo", "start", "end")], col.names=F, row.names=F, quote=F, sep="\t")

cat("\n")

write.table(matchDT[,c("chromo", "start", "end", "tad")], col.names=F, row.names=F, quote=F, sep="\t")
