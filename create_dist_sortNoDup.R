startTime <- Sys.time()
cat(paste0("> Rscript create_dist_sortNoDup.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

caller ="TopDom"

# Rscript create_dist_sortNoDup.R

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> add as.character() in apply !!! use "gene1" etc. instead of 1 index in apply 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_DIST_SORTNODUP")
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
  #### UPDATE 30.06.2018 -> SORT TO ENSURE ALPHABETICAL ORDER
  chromo_g2t_dt <- chromo_g2t_dt[order(as.character(chromo_g2t_dt$entrezID)),]
  chromoDT <- as.data.frame(t(combn(chromo_g2t_dt$entrezID, m=2)))
  colnames(chromoDT) <- c("gene1", "gene2")
  #### UPDATE 30.06.2018 -> ENSURE CHARACTER 
  chromoDT$gene1 <- as.character(chromoDT$gene1)
  chromoDT$gene2 <- as.character(chromoDT$gene2)
  #### UPDATE 30.06.2018 -> ENSURE CHARACTER IN APPLY !!! AND NAME "GENE1" AND "GENE2" INSTEAD OF INDICES TO ENSURE
  genePairDist <- apply(chromoDT, 1, function(x) 
    abs(chromo_g2t_dt$midPos[as.character(chromo_g2t_dt$entrezID) == as.character(x["gene1"])] - 
          chromo_g2t_dt$midPos[as.character(chromo_g2t_dt$entrezID) == as.character(x["gene2"])]))
  stopifnot(length(genePairDist) == nrow(chromoDT))
  chromoDT$chromo <- chromo
  chromoDT$dist <- genePairDist
  #### UPDATE 30.06.2018 -> CHECK TO ENSURE ALPHABETICAL ORDER
  stopifnot(chromoDT$gene1 < chromoDT$gene2)
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
