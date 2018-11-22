startTime <- Sys.time()
cat(paste0("> Rscript prep_all_genes_data_v2.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

caller ="TopDom"
familyData ="hgnc"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "PREP_ALL_GENES_v2")
system(paste0("mkdir -p ", outFold))
 
gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT$midPos <- (gene2tadDT$start + gene2tadDT$end)/2
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TAD_DT <- read.table(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"))
TAD_DT <- TAD_DT[grepl("_TAD", TAD_DT$region),]

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)

# 
# gene2tadDT <- gene2tadDT[1:10,]

all_pairs <- combn(gene2tadDT$entrezID,m=2)

all_familyData <- paste0(familyData, c("_family", "_family_short"))

for(i_fam in all_familyData) {
  

  pairsDT <- foreach(i = 1:ncol(all_pairs), .combine="rbind") %dopar% {
    
    cat(paste0("... ", i, "/", ncol(all_pairs), "\n"))
    gene1 <- all_pairs[1,i]
    gene2<- all_pairs[2,i]
    i_1 <- which(gene2tadDT$entrezID == gene1)
    i_2 <- which(gene2tadDT$entrezID == gene2)
    
    
    # foreach(i_1 = 1:(nrow(gene2tadDT)-1), .combine='rbind') %:%
    #      foreach(i_2 = (i_1+1):nrow(gene2tadDT), .combine='rbind') %dopar% {
    #        
    #        cat(paste0("... ", i_1, "/", nrow(gene2tadDT), " - ", i_2, "/", nrow(gene2tadDT), "\n"))
           
           # gene1 <- gene2tadDT$entrezID[i_1]
           # gene2 <- gene2tadDT$entrezID[i_2]
           
           distGenes <- ifelse(gene2tadDT$chromo[i_1] == gene2tadDT$chromo[i_2], 
                               abs(gene2tadDT$midPos[i_1]-gene2tadDT$midPos[i_2]), NA)
           
           sameTAD <- as.numeric(gene2tadDT$region[i_1] == gene2tadDT$region[i_2])
           
           sameFamily <- ifelse(gene1 %in% familyDT$entrezID & gene2 %in% familyDT$entrezID, 
                                as.numeric(familyDT[familyDT$entrezID == gene1, i_fam] == familyDT[familyDT$entrezID == gene2, i_fam]),
                                NA)
           data.frame(
             gene1=gene1,
             gene2=gene2,
             distGenes=distGenes,
             sameTAD=sameTAD,
             sameFamily=sameFamily,
             stringsAsFactors = FALSE
           )
         }
  stopifnot(nrow(pairsDT) == (0.5*nrow(gene2tadDT)*(nrow(gene2tadDT)-1)))
  
  outFile <- file.path(outFold, "/", paste0(i_fam, "_pairsDT.Rdata"))
  save(pairsDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
}



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))