startTime <- Sys.time()
cat(paste0("> Rscript create_dist.R\n"))

stop("!!! use _sortNoDup script !!!\n")

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

caller ="TopDom"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_DIST")
system(paste0("mkdir -p ", outFold))
 
gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT$midPos <- (gene2tadDT$start + gene2tadDT$end)/2
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]

all_chromo <- unique(gene2tadDT$chromo)

all_dist_pairs <- foreach(chromo = all_chromo, .combine='rbind') %dopar% {
  chromo_g2t_dt <- gene2tadDT[gene2tadDT$chromo == chromo,]
  if(nrow(chromo_g2t_dt) == 1) return(NULL)
  chromoDT <- as.data.frame(t(combn(chromo_g2t_dt$entrezID, m=2)))
  colnames(chromoDT) <- c("gene1", "gene2")
  genePairDist <- apply(chromoDT, 1, function(x) 
    abs(chromo_g2t_dt$midPos[chromo_g2t_dt$entrezID == x[1]] - chromo_g2t_dt$midPos[chromo_g2t_dt$entrezID == x[2]]))
  stopifnot(length(genePairDist) == nrow(chromoDT))
  chromoDT$chromo <- chromo
  chromoDT$dist <- genePairDist
  chromoDT
}

outFile <- file.path(outFold, "all_dist_pairs.Rdata")
save(all_dist_pairs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
