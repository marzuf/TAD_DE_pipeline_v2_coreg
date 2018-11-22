suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

startTime <- Sys.time()

cat(paste0("> Rscript pair_counts_Fisher.R\n"))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}


### HARD CODED
buildTable <- TRUE
caller <- "TopDom"

### RETRIEVE FROM COMMAND LINE
# Rscript pair_counts_Fisher.R
# Rscript pair_counts_Fisher.R TCGAcrc_msi_mss 50 hgnc
#  Rscript pair_counts_Fisher.R <dataset> <nTADs> <family>
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 3)
  curr_dataset <- args[1]
  nTopTADs <- as.numeric(args[2])
  stopifnot(!is.na(nTopTADs))
  familyData <- args[3]
} else{
  txt <- paste0("> Default parameters:\n")
  nTopTADs <- 50
  curr_dataset <- "TCGAcrc_msi_mss"
  familyData <- "hgnc"
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "PAIR_COUNTS_FISHER",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("pair_counts_fisher.txt"))  
system(paste0("rm -f ", logFile))

printAndLog(txt, logFile)
txt <- paste0("... nTopTADs = ",  nTopTADs, "\n")
printAndLog(txt, logFile)
txt <- paste0("... curr_dataset = ",  curr_dataset, "\n")
printAndLog(txt, logFile)
txt <- paste0("... familyData = ",  familyData, "\n")
printAndLog(txt, logFile)

txt <- paste0("> ! Hard-coded parameters:\n")
printAndLog(txt, logFile)
txt <- paste0("... caller = ",  caller, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

toprankingScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking")
source(paste0(toprankingScriptDir, "/", "get_topTADs.R"))

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)

topTADs <- get_topTADs(mydataset = curr_dataset, nTop = nTopTADs, TADonly=TRUE)

#============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
script0_name <- "0_prepGeneData"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline: ", length(pipeline_geneList), "\n")
printAndLog(txt, logFile)


curr_geneList <- as.character(pipeline_geneList[pipeline_geneList %in% familyDT$entrezID])
txt <- paste0("...... number of genes with family data available: ", length(curr_geneList), "\n")
printAndLog(txt, logFile)


curr_gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList, ]
rm(gene2tadDT)

all_familyData <- paste0(familyData, c("_family", "_family_short"))

cat("... form pair of genes\n")

stopifnot(curr_geneList %in% curr_gene2tadDT$entrezID)
all_gene_pairs <- combn(x = curr_geneList, m = 2)
all_gene_pairs_DT_s <- data.frame(gene1 = all_gene_pairs[1,],
                                gene2 = all_gene_pairs[2,],
                                stringsAsFactors = FALSE)

for(i_fam in all_familyData) {
  cat(paste0("> START: ", i_fam, "\n"))
  
  if(buildTable) {
    
    all_gene_pairs_DT <- all_gene_pairs_DT_s
    
    # retrieve TAD for gene1
    cat("...... retrieve TAD for gene1\n")
    all_gene_pairs_DT <- left_join(all_gene_pairs_DT, curr_gene2tadDT[,c("entrezID", "region")], by=c("gene1" = "entrezID"))
    colnames(all_gene_pairs_DT)[colnames(all_gene_pairs_DT) == "region"] <- "tad1"
    
    cat("...... retrieve TAD for gene2v")
    all_gene_pairs_DT <- left_join(all_gene_pairs_DT, curr_gene2tadDT[,c("entrezID", "region")], by=c("gene2" = "entrezID"))
    colnames(all_gene_pairs_DT)[colnames(all_gene_pairs_DT) == "region"] <- "tad2"
    
    cat("... retrieve same TAD gene1-gene2\n")
    all_gene_pairs_DT$sameTAD <- as.numeric(all_gene_pairs_DT$tad1 == all_gene_pairs_DT$tad2 )
    
    cat("... retrieve topTAD gene1-gene2\n")
    all_gene_pairs_DT$topTAD <- as.numeric(all_gene_pairs_DT$tad1 %in% topTADs & 
                                             all_gene_pairs_DT$tad2  %in% topTADs &
                                             all_gene_pairs_DT$sameTAD)
    
    cat("...... retrieve family for gene1\n")
    all_gene_pairs_DT <- left_join(all_gene_pairs_DT, familyDT[,c("entrezID", i_fam)], by=c("gene1" = "entrezID"))
    colnames(all_gene_pairs_DT)[colnames(all_gene_pairs_DT) == i_fam] <- paste0(i_fam, "1")
    
    cat("...... retrieve family for gene2\n")
    all_gene_pairs_DT <- left_join(all_gene_pairs_DT, familyDT[,c("entrezID", i_fam)], by=c("gene2" = "entrezID"))
    colnames(all_gene_pairs_DT)[colnames(all_gene_pairs_DT) == i_fam] <- paste0(i_fam, "2")
    
    cat("... retrieve same family gene1-gene2\n")
    all_gene_pairs_DT$sameFamily <- as.numeric(all_gene_pairs_DT[,paste0(i_fam, "1")] == all_gene_pairs_DT[,paste0(i_fam, "2")])
    
    outFile <- file.path(outFold, paste0(i_fam, "_all_gene_pairs_DT.Rdata"))
    save(all_gene_pairs_DT, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
  } else {
    outFile <- file.path(outFold, paste0(i_fam, "_all_gene_pairs_DT.Rdata"))
    all_gene_pairs_DT <- eval(parse(text= load(outFile)))
    }
}
  

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


