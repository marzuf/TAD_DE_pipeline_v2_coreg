options(scipen = 100)
# !!! need to change that gene2family also stores NA !!!
suppressPackageStartupMessages(library(WGCNA, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

cat(paste0("> Rscript pair_counts_Fisher_heatmap.R\n"))

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

caller <- "TopDom"

startTime <- Sys.time()

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
if(!SSHFS) allowWGCNAThreads()
registerDoMC(ifelse(SSHFS, 2, 40))

source(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking/get_topTADs.R")))
source(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking/get_topGenes.R")))

### HARD CODED
caller <- "TopDom"
buildTable <- TRUE

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

### RETRIEVE FROM COMMAND LINE
# Rscript pair_counts_Fisher_heatmap.R
# Rscript pair_counts_Fisher_heatmap.R TCGAcrc_msi_mss 50 10 hgnc
#  Rscript pair_counts_Fisher_heatmap.R <dataset> <nTADs> <nTADs_plot <family>
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 4)
  curr_dataset <- args[1]
  nTopTADs <- as.numeric(args[2])
  stopifnot(!is.na(nTopTADs))
  nTopTADs_plot <- as.numeric(args[3])
  stopifnot(!is.na(nTopTADs_plot))
  familyData <- args[4]
} else{
  txt <- paste0("> Default parameters:\n")
  nTopTADs <- 50
  nTopTADs_plot <- 10
  curr_dataset <- "TCGAcrc_msi_mss"
  familyData <- "hgnc"
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "PAIR_COUNTS_FISHER_HEATMAP",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("pair_counts_fisher_heatmap.txt"))  
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

inFoldAssignment <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "PAIR_COUNTS_FISHER",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))

# for the heatmap Fisher tests
topTADs_plot <- get_topTADs(mydataset = curr_dataset, nTop = nTopTADs_plot, TADonly = TRUE)

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)

gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))


#============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
script0_name <- "0_prepGeneData"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline: ", length(pipeline_geneList), "\n")
printAndLog(txt, logFile)

all_familyData <- paste0(familyData, c("_family", "_family_short"))

for(i_fam in all_familyData) {
  txt <- paste0("> START: ", i_fam, "\n")
  printAndLog(txt, logFile)
  cat("...... load assignment data \n")  
  assignmentDT <- eval(parse(text = load(file.path(inFoldAssignment, paste0(i_fam, "_all_gene_pairs_DT.Rdata")))))
  stopifnot(topTADs_plot %in% assignmentDT$tad1, topTADs_plot %in% assignmentDT$tad2)
  stopifnot(assignmentDT$gene1 %in% familyDT$entrezID, assignmentDT$gene2 %in% familyDT$entrezID)
  stopifnot(assignmentDT$gene1 %in% gene2tadDT$entrezID, assignmentDT$gene2 %in% gene2tadDT$entrezID)
  
  txt <- paste0("...... # gene pairs retrieved: ", nrow(assignmentDT), "\n")
  printAndLog(txt, logFile)
  
  #===== FISHER TESTS SAME/DIFF TADs VS. SAME/DIFF FAM <================================ ALL DATA
  
  txt <- paste0("...... Fisher test: sameTADs/diffTADs vs. sameFam/diffFam \n")
  printAndLog(txt, logFile)
  
  # do Fisher TADs for all the data
  n_sameTAD_sameFam <- sum(assignmentDT$sameTAD == 1 & assignmentDT$sameFamily == 1)
  n_sameTAD_diffFam <- sum(assignmentDT$sameTAD == 1 & assignmentDT$sameFamily == 0)
  
  n_diffTAD_sameFam <- sum(assignmentDT$sameTAD == 0 & assignmentDT$sameFamily == 1)
  n_diffTAD_diffFam <- sum(assignmentDT$sameTAD == 0 & assignmentDT$sameFamily == 0)
  
  txt <- paste0("...... # gene pairs sameTAD + sameFam: ", n_sameTAD_sameFam, "/", nrow(assignmentDT), " (", round(n_sameTAD_sameFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("...... # gene pairs sameTAD + diffFam: ", n_sameTAD_diffFam, "/", nrow(assignmentDT), " (", round(n_sameTAD_diffFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("...... # gene pairs diffTAD + sameFam: ", n_diffTAD_sameFam, "/", nrow(assignmentDT), " (", round(n_diffTAD_sameFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("...... # gene pairs diffTAD + diffFam: ", n_diffTAD_diffFam, "/", nrow(assignmentDT), " (", round(n_diffTAD_diffFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  testMat <- matrix(c(n_sameTAD_sameFam, n_sameTAD_diffFam, n_diffTAD_sameFam, n_diffTAD_diffFam), nrow=2,
         byrow = FALSE, dimnames = list(c("sameTAD","diffTAD"), c("sameFam","diffFam")))
  
  stopifnot(sum(testMat) == nrow(assignmentDT))
  
  sink(logFile, append = TRUE)
  print(testMat)
  sink()
  
  sameDiffTADs_sameDiffFam_FT <- fisher.test(testMat, alternative = "greater")$p.value
  txt <- paste0("... Fisher's test p-val = ", sprintf("%.4f", sameDiffTADs_sameDiffFam_FT), "\n")
  printAndLog(txt, logFile)
  
  #===== FISHER TESTS top/other TADs VS. SAME/DIFF FAM <================================ ALL DATA
  
  txt <- paste0("...... Fisher test: topTADs/otherTADs vs. sameFam/diffFam \n")
  printAndLog(txt, logFile)
  
  # do Fisher TADs for all the data
  n_topTAD_sameFam <- sum(assignmentDT$topTAD == 1 & assignmentDT$sameFamily == 1)
  n_topTAD_diffFam <- sum(assignmentDT$topTAD == 1 & assignmentDT$sameFamily == 0)
  
  n_otherTAD_sameFam <- sum(assignmentDT$topTAD == 0 & assignmentDT$sameFamily == 1)
  n_otherTAD_diffFam <- sum(assignmentDT$topTAD == 0 & assignmentDT$sameFamily == 0)
  
  txt <- paste0("...... # gene pairs topTAD + sameFam: ", n_topTAD_sameFam, "/", nrow(assignmentDT), " (", round(n_topTAD_sameFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("...... # gene pairs topTAD + diffFam: ", n_topTAD_diffFam, "/", nrow(assignmentDT), " (", round(n_topTAD_diffFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("...... # gene pairs otherTAD + sameFam: ", n_otherTAD_sameFam, "/", nrow(assignmentDT), " (", round(n_otherTAD_sameFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("...... # gene pairs otherTAD + diffFam: ", n_otherTAD_diffFam, "/", nrow(assignmentDT), " (", round(n_otherTAD_diffFam/nrow(assignmentDT)*100, 2), "%)\n")
  printAndLog(txt, logFile)
  
  testMat <- matrix(c(n_topTAD_sameFam, n_topTAD_diffFam, n_otherTAD_sameFam, n_otherTAD_diffFam), nrow=2,
                    byrow = FALSE, dimnames = list(c("topTAD","otherTAD"), c("sameFam","diffFam")))
  
  stopifnot(sum(testMat) == nrow(assignmentDT))
  
  sink(logFile, append = TRUE)
  print(testMat)
  sink()
  
  topOtherTADs_sameDiffFam_FT <- fisher.test(testMat, alternative = "greater")$p.value
  txt <- paste0("... Fisher's test p-val = ", sprintf("%.4f", topOtherTADs_sameDiffFam_FT), "\n")
  printAndLog(txt, logFile)
  
  ### plot heatmap for the nTopTADs_plot
  # mod1 = TADs
  # mod2 = families
  
  topGenes_plot <- get_topGenes(mydataset = curr_dataset, nTop = nTopTADs_plot, TADonly = TRUE) 
  mod1Labels <- unique(as.character(sapply(topGenes_plot, function(x) gene2tadDT$region[gene2tadDT$entrezID == x])))
  mod1Order <- get_topTADs(mydataset = curr_dataset, nTop = nTopTADs_plot, TADonly = TRUE) 
  
  # gene2mod1_topTADs <- sapply(topGenes_plot, function(x) gene2tadDT$region[gene2tadDT$entrezID == x])
  # gene2mod2_topTADs <- sapply(topGenes_plot, function(x) ifelse(x %in% familyDT$entrezID, familyDT[,i_fam][familyDT$entrezID == x], "unassigned"))
  # universe should be all families and all genes from the pipeline  
  gene2mod1_topTADs <- sapply(pipeline_geneList, function(x) gene2tadDT$region[gene2tadDT$entrezID == x])
  gene2mod2_topTADs <- sapply(pipeline_geneList, function(x) 
    ifelse(x %in% familyDT$entrezID, familyDT[,i_fam][familyDT$entrezID == x], "unassigned"))
  
  mod1Labels_s <- mod1Labels[match(mod1Order[mod1Order %in% mod1Labels], mod1Labels)]
  stopifnot(setequal(mod1Labels, mod1Labels_s))
  stopifnot(length(mod1Labels) == length(mod1Labels_s))
  stopifnot(!is.na(mod1Labels_s))
  mod1Labels <- mod1Labels_s
  
  
  mod2Labels <- unique(as.character(sapply(topGenes_plot, function(x) ifelse(x %in% familyDT$entrezID, familyDT[,i_fam][familyDT$entrezID == x], "unassigned"))))
  mod2Order <- c(names(sort(table(familyDT[,i_fam]), decreasing = T)), "unassigned")
  
  mod2Labels_s <- mod2Labels[match(mod2Order[mod2Order %in% mod2Labels], mod2Labels)]
  stopifnot(setequal(mod2Labels, mod2Labels_s))
  stopifnot(length(mod2Labels) == length(mod2Labels_s))
  stopifnot(!is.na(mod2Labels_s))
  mod2Labels <- mod2Labels_s
  
  type_mod2 <- paste0(i_fam)
  type_mod1 <- paste0("topTADs", nTopTADs_plot)
  cond1_name <- type_mod1
  cond2_name <- type_mod2
  
  
  # Numbers of female and consensus modules
  nMod1  <- length(mod1Labels) # the TADS
  nMod2 <- length(mod2Labels) # the families
  
  # Initialize tables of p-values and of the corresponding counts
  pTable <- matrix(0, nrow = nMod1, ncol = nMod2);
  CountTbl <- matrix(0, nrow = nMod1, ncol = nMod2);
  
  stopifnot(names(gene2mod1_topTADs) == names(gene2mod2_topTADs))
  
  # Execute all pairwaise comparisons
  # outerloop -> y axis heatmap
  cat("... start computing pairwise mod1-mod2 module comparisons (Fisher tests) for topTADs - can take a while \n")
  cat("start:", paste0(Sys.time()), "\tend:")
  
  mod1_mod2_testDT_topTADs <-   foreach(i_mod1=1:nMod1, .combine='rbind') %:%
    foreach(i_mod2=1:nMod2, .combine='rbind') %dopar% {
      mod1Members <- ( as.character(gene2mod1_topTADs) == mod1Labels[i_mod1] )
      mod2Members <- ( as.character(gene2mod2_topTADs)  == mod2Labels[i_mod2])
      
      # universe should be full set of genes
      stopifnot(length(mod1Members) == length(mod2Members))
      stopifnot(length(mod1Members) == length(pipeline_geneList))
      
      pval <- -log10(fisher.test(mod1Members, mod2Members, alternative = "greater")$p.value);
      count <- sum( ( as.character(gene2mod1_topTADs) == mod1Labels[i_mod1] ) & 
                      ( as.character(gene2mod2_topTADs)  == mod2Labels[i_mod2]))
      data.frame(i_mod1 = i_mod1, 
                 i_mod2 = i_mod2,
                 pval = pval,
                 count = count)
    }
  pTable <- as.matrix(sparseMatrix( i = mod1_mod2_testDT_topTADs$i_mod1, j = mod1_mod2_testDT_topTADs$i_mod2, x = mod1_mod2_testDT_topTADs$pval  ))
  CountTbl <- as.matrix(sparseMatrix( i = mod1_mod2_testDT_topTADs$i_mod1, j = mod1_mod2_testDT_topTADs$i_mod2, x = mod1_mod2_testDT_topTADs$count  ))
  
  cat(paste0(Sys.time()), "\n")
  
  
  ##### VISUALIZE THE COMPARISON
  # To display the p-value and count tables in an informative way, create a color-coded table of the intersection
  # counts (colors indicate the p-value significance)
  # Truncate p values smaller than 10^{-50} to 10^{-50} 
  pTable[is.infinite(pTable)] <- 1.3*max(pTable[is.finite(pTable)]);
  pTable[pTable>50 ] <- 50 ;
  # Marginal counts (really module sizes)
  Cond1ModTotals <- apply(CountTbl, 1, sum)
  Cond2ModTotals <- apply(CountTbl, 2, sum)
  
  # # Actual plotting
  outFile <- file.path(outFold, paste0(i_fam, "_fisherTest_", cond1_name,"_", cond2_name, "_heatmap.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width=myWidth))
  par(mfrow=c(1,1))
  par(cex = 1.0)
  par(mar=c(10, 10.4, 2.7, 1)+0.3)
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  # outerloop = y axis
  
  empty_CountTbl <- matrix(rep("", length(CountTbl)), nrow = dim(CountTbl)[1])
  
  labeledHeatmap(Matrix = pTable,
                 xLabels = paste("", mod2Labels),
                 yLabels = paste("", mod1Labels),
                 colorLabels = TRUE,
                 xSymbols = paste(mod2Labels, ": ", Cond2ModTotals, sep=""),
                 ySymbols = paste(mod1Labels, ": ", Cond1ModTotals, sep=""),
                 textMatrix = CountTbl,
                 # textMatrix = empty_CountTbl,
                 colors = greenWhiteRed(100)[50:100],
                 main = paste0("Correspondence of ", type_mod2, " - ", type_mod1, " assignment"),
                 verticalSeparator.x = 1:ncol(pTable),
                 cex.text = 1.0, 
                 # cex.lab = 1.0, 
                 cex.lab.x = 0.7,
                 cex.lab.y = 0.7,
                 setStdMargins = FALSE)
  mtext(type_mod2, side=1, line=7)
  mtext(type_mod1, side=2, line=5)
  mtext(paste0(curr_dataset, " - ", i_fam), side=3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}



######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


