startTime <- Sys.time()
cat(paste0("> Rscript create_combn.R\n"))

#suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
#suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

caller ="TopDom"
familyData ="hgnc"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
#registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_COMBN")
system(paste0("mkdir -p ", outFold))
 
gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT$midPos <- (gene2tadDT$start + gene2tadDT$end)/2
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]


# 
# gene2tadDT <- gene2tadDT[1:10,]

# all_gene_pairs <- combn(gene2tadDT$entrezID,m=2)

outFile = file.path(outFold, "all_gene_pairs.Rdata")
# save(all_gene_pairs, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))

cat(paste0("... start loading: ", Sys.time(), "\n"))
load(outFile)
cat(paste0("... end loading: ", Sys.time(), "\n"))

# cat(paste0("... start convert DT: ", Sys.time(), "\n"))
# all_gene_pairs_DT <- as.data.frame(t(all_gene_pairs))
# colnames(all_gene_pairs_DT) <- c("gene1", "gene2")
# outFile = file.path(outFold, "all_gene_pairs_DT.Rdata")
# save(all_gene_pairs_DT, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))




######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
