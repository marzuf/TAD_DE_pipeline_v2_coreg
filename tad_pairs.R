startTime <- Sys.time()
cat(paste0("> Rscript tad_pairs.R\n"))

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

### HARD_CODED
caller <- "TopDom"
buildTable <- FALSE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 2 * 10^6
nbrBins <- 10

### RETRIEVE FROM COMMAND LINE
# Rscript tad_pairs.R
# Rscript tad_pairs.R TCGAcrc_msi_mss 50 hgnc
#  Rscript tad_pairs.R <dataset> <nTADs> <family>
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

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "TAD_PAIRS", paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("tad_pairs_logFile_boxplotUpdate.txt"))  
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
txt <- paste0("... distLimit = ",  distLimit, "\n")
printAndLog(txt, logFile)
txt <- paste0("... nbrBins = ",  nbrBins, "\n")
printAndLog(txt, logFile)


col1 <- "red"
col2 <- "green"
col3 <- "blue"

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

toprankingScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking")
source(paste0(toprankingScriptDir, "/", "get_topTADs.R"))

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT$midPos <- (gene2tadDT$start + gene2tadDT$end)/2

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))

familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)
# for all families
# -> for all gene pairs
#  --> distance
#  --> coexpression
#  --> same TAD ?

topTADs <- get_topTADs(mydataset = curr_dataset, nTop = nTopTADs, TADonly = TRUE)

#============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline: ", length(pipeline_geneList), "\n")
printAndLog(txt, logFile)
stopifnot(pipeline_geneList %in% gene2tadDT$entrezID)

pipeline_regionList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_regionList.Rdata"))))
pipeline_regionList <- pipeline_regionList[grepl("_TAD", pipeline_regionList)]
txt <- paste0("... number of TADs used in the ", curr_dataset, " pipeline: ", length(pipeline_regionList), "\n")
printAndLog(txt, logFile)
stopifnot(pipeline_regionList %in% gene2tadDT$region)

DE_dt <- eval(parse(text = load(file.path(dataset_pipDir, script1_name, "DE_topTable.Rdata"))))
stopifnot(DE_dt$genes %in% rownames(DE_dt))
stopifnot(names(pipeline_geneList) %in% rownames(DE_dt))
stopifnot(names(pipeline_geneList) %in% DE_dt$genes)
# but not pipeline_geneList
# all(rownames(DE_dt) %in% gene2tadDT$entrezID) => FALSE
DE_dt <- DE_dt[DE_dt$genes %in% names(pipeline_geneList),]
rownames(DE_dt) <- unlist(sapply(DE_dt$genes, function(x) pipeline_geneList[names(pipeline_geneList)==x] ))
stopifnot(all(pipeline_geneList %in% rownames(DE_dt)))
DE_dt$genes <- rownames(DE_dt)

# take only genes from the pipeline
ds_familyDT <- familyDT[familyDT$entrezID %in% pipeline_geneList,]
# just ensure that I kept only the "TRUE" family data (not the single "unassigned" family)
stopifnot(!any(grepl("assigned", ds_familyDT$hgnc_family_short)))

ds_gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList,]

rm(gene2tadDT)

txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline with family info.: ", nrow(ds_familyDT), "\n")
printAndLog(txt, logFile)

all_familyData <- paste0(familyData, c("_family", "_family_short"))

for(i_fam in all_familyData) {
  
  cat(paste0("> START: ", i_fam, "\n"))
  
  all_tads <- unique(pipeline_regionList)
  cat(paste0("... # TADs: ", length(all_tads), "\n"))
  
  if(buildTable) {
    tad_dist_logFC_DT <- foreach(curr_tad = all_tads, .combine='rbind') %dopar% {
      tad_genes <- ds_gene2tadDT$entrezID[ds_gene2tadDT$region == curr_tad]  
      stopifnot(tad_genes %in% pipeline_geneList)
      if(length(tad_genes) == 1) {
        return(data.frame(
          gene1 = tad_genes,
          gene2 = tad_genes,
          dist = NA,
          prodLogFC = NA,
          sameFamily = NA,
          TAD = curr_tad,
          stringsAsFactors = FALSE
        ))
      }
      # create combn of gene pairs
      genePairs <- combn(x = tad_genes, m = 2)
      stopifnot(ncol(genePairs) == (length(tad_genes) * (length(tad_genes) -1) *0.5))
      
      genepairs_DT <- foreach(i_pair = 1:ncol(genePairs), .combine='rbind') %do% {
        gene1 <- genePairs[1,i_pair]
        gene2 <- genePairs[2,i_pair]
        chromo1 <- ds_gene2tadDT$chromo[ds_gene2tadDT$entrezID == gene1]
        chromo2 <- ds_gene2tadDT$chromo[ds_gene2tadDT$entrezID == gene2]
        stopifnot(length(chromo1) == length(chromo2))
        stopifnot(length(chromo1) == 1)
        pair_dist <- ifelse(chromo1 == chromo2, 
                            abs(ds_gene2tadDT$midPos[ds_gene2tadDT$entrezID == gene1] - ds_gene2tadDT$midPos[ds_gene2tadDT$entrezID == gene2]), NA)
        
        
        pair_sameFamily <- ifelse( gene1 %in% ds_familyDT$entrezID & gene2 %in% ds_familyDT$entrezID,
                                   as.numeric( ds_familyDT[ds_familyDT$entrezID == gene1,i_fam] == ds_familyDT[ds_familyDT$entrezID == gene2,i_fam]),
                                   NA)
        pair_prodLogFC <- DE_dt[gene1,]$logFC * DE_dt[gene2,]$logFC 
        data.frame(
          gene1 = gene1,
          gene2 = gene2,
          dist = pair_dist,
          prodLogFC = pair_prodLogFC,
          sameFamily = pair_sameFamily,
          TAD = curr_tad,
          stringsAsFactors = FALSE
        )
      } # end iterating over gene pairs
      genepairs_DT
    } # end iterating over the tads
    outFile <- file.path(outFold, paste0(i_fam, "_tad_dist_logFC_DT.Rdata"))
    save(tad_dist_logFC_DT, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
  } else {
    outFile <- file.path(outFold, paste0(i_fam, "_tad_dist_logFC_DT.Rdata"))
    # outFile <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/TAD_PAIRS/TCGAcrc_msi_mss/hgnc_family_tad_dist_logFC_DT.Rdata"
    # outFile <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/TAD_PAIRS/TCGAcrc_msi_mss_nTop50_hgnc/hgnc_family_short_tad_dist_logFC_DT.Rdata"
    tad_dist_logFC_DT <- eval(parse(text = load(outFile)))
  }
  txt <- paste0("... total # of gene pairs: ", nrow(tad_dist_logFC_DT), "\n")
  printAndLog(txt, logFile)
  
  tad_dist_logFC_DT <- tad_dist_logFC_DT[!is.na(tad_dist_logFC_DT$sameFamily),]
  stopifnot(!is.na(tad_dist_logFC_DT))
  txt <- paste0("... total # of gene pairs with family annotation: ", nrow(tad_dist_logFC_DT), "\n")
  printAndLog(txt, logFile)
  
  binSize <- distLimit/nbrBins
  tad_dist_logFC_DT$distBin <- floor(tad_dist_logFC_DT$dist/binSize)
  
  tad_dist_logFC_DT$topTAD <- unlist(sapply(tad_dist_logFC_DT$TAD, function(x) as.numeric(x %in% topTADs)))
  
  aggBin_tad_dist_logFC_DT <- aggregate(rep(1, nrow(tad_dist_logFC_DT)), 
                                        by = list(sameFamily = tad_dist_logFC_DT$sameFamily, distBin= tad_dist_logFC_DT$distBin), sum)
  
  stopifnot(sum(tad_dist_logFC_DT$distBin== 0 & tad_dist_logFC_DT$sameFamily == 0) == 
              aggBin_tad_dist_logFC_DT$x[aggBin_tad_dist_logFC_DT$sameFamily == 0 & aggBin_tad_dist_logFC_DT$distBin == 0 ])
  
  stopifnot(sum(tad_dist_logFC_DT$distBin== 1 & tad_dist_logFC_DT$sameFamily == 1) == 
              aggBin_tad_dist_logFC_DT$x[aggBin_tad_dist_logFC_DT$sameFamily == 1 & aggBin_tad_dist_logFC_DT$distBin == 1 ])
  
  
  topTAD_aggBin_tad_dist_logFC_DT <- aggregate(rep(1, nrow(tad_dist_logFC_DT)), 
                                        by = list(sameFamily = tad_dist_logFC_DT$sameFamily, distBin= tad_dist_logFC_DT$distBin, topTAD = tad_dist_logFC_DT$topTAD),
                                        sum)
  
  stopifnot(sum(tad_dist_logFC_DT$distBin== 0 & tad_dist_logFC_DT$sameFamily == 0 & tad_dist_logFC_DT$topTAD == 0) == 
              topTAD_aggBin_tad_dist_logFC_DT$x[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 0 & topTAD_aggBin_tad_dist_logFC_DT$distBin == 0 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0 ])
  
  stopifnot(sum(tad_dist_logFC_DT$distBin== 1 & tad_dist_logFC_DT$sameFamily == 1 & tad_dist_logFC_DT$topTAD == 1) == 
              topTAD_aggBin_tad_dist_logFC_DT$x[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$distBin == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1 ])
  
  
  my_xlab <- paste0("Distance bin (", distLimit, "/", nbrBins, ") between the 2 genes")
  my_ylab <- paste0("Nbr of gene pairs")
  
  #============================================================= PLOT GENE PAIRS FROM SAME TAD, SAME FAMILY VS. DIFF. FAMILY [ALL TADs]
  
  
  outFile <- file.path(outFold, paste0(i_fam, "_same_TAD_same_family_vs_diff_family.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  plot(NULL,
       xlim = c(0, max(aggBin_tad_dist_logFC_DT$distBin)),
       ylim = range(aggBin_tad_dist_logFC_DT$x),
       xlab = my_xlab,
       ylab = my_ylab,
       main ="Gene pairs from same TAD")
  
  mtext(text = paste0(curr_dataset, " - ", i_fam, " - same family vs. diff. families - all TADs"), side=3)
  
  lines(x ~ distBin, aggBin_tad_dist_logFC_DT[aggBin_tad_dist_logFC_DT$sameFamily == 1,],
          col = col1)
  points(x ~ distBin, aggBin_tad_dist_logFC_DT[aggBin_tad_dist_logFC_DT$sameFamily == 1,],
        col = col1, pch=16, cex=0.7)
  lines(x ~ distBin, aggBin_tad_dist_logFC_DT[aggBin_tad_dist_logFC_DT$sameFamily == 0,],
          col = col2)
  points(x ~ distBin, aggBin_tad_dist_logFC_DT[aggBin_tad_dist_logFC_DT$sameFamily == 0,],
        col = col2, pch=16, cex=0.7)
  
  xidx <- aggBin_tad_dist_logFC_DT$sameFamily == 1
  yidx <- aggBin_tad_dist_logFC_DT$sameFamily == 0
  
  legTxt <- c(paste0("same TAD + same family (# pairs = ", sum(aggBin_tad_dist_logFC_DT$x[xidx]), ")"), 
              paste0("same TAD + diff. families (# pairs = ", sum(aggBin_tad_dist_logFC_DT$x[yidx]), ")"))
  
  # legTxt <- c(paste0("same TAD + same family (# pairs = ", sum(aggBin_tad_dist_logFC_DT$sameFamily == 1), ")"), 
  #             paste0("same TAD + diff. families (# pairs = ", sum(aggBin_tad_dist_logFC_DT$sameFamily == 0), ")"))
  # 
  legend("topright", legend = legTxt, bty="n", lty=2, col = c(col1, col2))
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  #============================================================= PLOT GENE PAIRS FROM SAME TAD, SAME FAMILY topTADs VS. SAME FAMILY otherTADs  VS. DIFF. FAMILY [ALL TADs]
  
  outFile <- file.path(outFold, paste0(i_fam, "_same_TAD_same_family_otherTADs_vs_same_family_diffTADs_vs_diff_family.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  plot(NULL,
       # xlim = c(0, max(topTAD_aggBin_tad_dist_logFC_DT$distBin)),
       # ylim = range(topTAD_aggBin_tad_dist_logFC_DT$x),
       # for the diff family -> needs to retain limits from the other DT
       xlim = c(0, max(aggBin_tad_dist_logFC_DT$distBin)),
       ylim = range(aggBin_tad_dist_logFC_DT$x),
       xlab = my_xlab,
       ylab = my_ylab,
       main ="Gene pairs from same TAD")
  
  mtext(text = paste0(curr_dataset, " - ", i_fam, " - same family topTADs vs. same family other TADs vs. diff. families"), side = 3)
  
  # same family + top TAD
  lines(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1,],
        col = col1)
  points(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1  & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1,],
         col = col1, pch=16, cex=0.7)
  
  # same family + other TAD
  lines(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0,],
        col = col2)
  points(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1  & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0,],
         col = col2, pch=16, cex=0.7)
  
  # diff family
  lines(x ~ distBin, aggBin_tad_dist_logFC_DT[aggBin_tad_dist_logFC_DT$sameFamily == 0,],
        col = col3)
  points(x ~ distBin, aggBin_tad_dist_logFC_DT[aggBin_tad_dist_logFC_DT$sameFamily == 0,],
         col = col3, pch=16, cex=0.7)
  
  # lines(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 0,],
  #       col = col3)
  # points(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 0,],
  #        col = col3, pch=16, cex=0.7)
  
  xidx <- topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1
  yidx <- topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0
  zidx <- topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 0
    
  legTxt <- c(paste0("same TAD + same family topTADs (# pairs = ", 
                     sum(topTAD_aggBin_tad_dist_logFC_DT$x[xidx]), ")"), 
              paste0("same TAD + same family otherTADs (# pairs = ",
                     sum(topTAD_aggBin_tad_dist_logFC_DT$x[yidx]) , ")"), 
              paste0("same TAD + diff. families (# pairs = ", 
                     sum(topTAD_aggBin_tad_dist_logFC_DT$x[zidx]) ,")"))
  
  # legTxt <- c(paste0("same TAD + same family topTADs (# pairs = ", 
  #                    sum(topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1), ")"), 
  #             paste0("same TAD + same family otherTADs (# pairs = ",
  #                    sum(topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0) , ")"), 
  #             paste0("same TAD + diff. families (# pairs = ", 
  #                    sum(topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 0) ,")"))
  
  legend("topright", legend = legTxt, bty="n", lty=2, col = c(col1, col2, col3))
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  #============================================================= PLOT GENE PAIRS FROM SAME TAD, SAME FAMILY topTADs VS. otherTADs 
  
  outFile <- file.path(outFold, paste0(i_fam, "_same_TAD_same_family_otherTADs_vs_diffTADs.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  plot(NULL,
       xlim = c(0, max(topTAD_aggBin_tad_dist_logFC_DT$distBin[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1])),
       ylim = range(topTAD_aggBin_tad_dist_logFC_DT$x[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1]),
       xlab = my_xlab,
       ylab = my_ylab,
       main ="Gene pairs from same TAD + same family")
  
  mtext(text = paste0(curr_dataset, " - ", i_fam, " - same family topTADs vs. otherTADs"),side = 3)
  
  # same family + top TAD
  lines(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1,],
        col = col1)
  points(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1  & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1,],
         col = col1, pch=16, cex=0.7)
  
  # same family + other TAD
  lines(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0,],
        col = col2)
  points(x ~ distBin, topTAD_aggBin_tad_dist_logFC_DT[topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1  & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0,],
         col = col2, pch=16, cex=0.7)
  
  xidx <- topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1
  yidx <- topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0
  
  legTxt <- c(paste0("same TAD + same family topTADs (# pairs = ", 
                     sum(topTAD_aggBin_tad_dist_logFC_DT$x[xidx]), ")"), 
              paste0("same TAD + same family otherTADs (# pairs = ",
                     sum(topTAD_aggBin_tad_dist_logFC_DT$x[yidx]) , ")"))
  
  # legTxt <- c(paste0("same TAD + same family topTADs (# pairs = ", 
  #                    sum(topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 1), ")"), 
  #             paste0("same TAD + same family otherTADs (# pairs = ",
  #                    sum(topTAD_aggBin_tad_dist_logFC_DT$sameFamily == 1 & topTAD_aggBin_tad_dist_logFC_DT$topTAD == 0) , ")"))
  
  legend("topright", legend = legTxt, bty="n", lty=2, col = c(col1, col2))
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  #============================================================= 
  
  my_ylab <- c("LogFC product")
  
  outFile <- file.path(outFold, paste0(i_fam, "_same_TAD_same_family_vs_diff_family_boxplot.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  
  # need factor to use add = TRUE ???
  tad_dist_logFC_DT$distBin <- factor(tad_dist_logFC_DT$distBin, levels = unique(sort(tad_dist_logFC_DT$distBin)))
  
  boxplot(prodLogFC ~ distBin, data = tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 1,],
          xlab = my_xlab,
          ylab = my_ylab,
          outline = FALSE)
  boxplot(prodLogFC ~ distBin, data = tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 0,], col ="red", add =T, outline=FALSE)
  
  # points(jitter(tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 1,]$distBin, factor=0.5)+1,
  #        tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 1,]$prodLogFC, pch=20, cex = 0.1)
  points(jitter(match( tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 1,]$distBin,levels(tad_dist_logFC_DT$distBin))-1, factor=0.5)+1,
         tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 1,]$prodLogFC, pch=20, cex = 0.1)
  # points(jitter(tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 0,]$distBin, factor=0.5)+1,
  #        tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 0,]$prodLogFC, pch=20, col ="red", cex = 0.1)
  points(jitter(match(tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 0,]$distBin,levels(tad_dist_logFC_DT$distBin))-1, factor=0.5)+1,
         tad_dist_logFC_DT[tad_dist_logFC_DT$sameFamily == 0,]$prodLogFC, pch=20, col ="red", cex = 0.1)
  
  title("logFC product between 2 genes from same TAD")
  mtext(paste0(curr_dataset, " - ", i_fam, " - same family vs. diff. families"), side = 3)
  

  legTxt <- c(paste0("same TAD + same family (# pairs = ", sum(tad_dist_logFC_DT$sameFamily == 1) , ")"),
              paste0("same TAD + diff. families (# pairs = ", sum(tad_dist_logFC_DT$sameFamily == 0), ")"))
  
  legend("bottomright", legend = legTxt, text.col = c("black", "red"), bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFold, paste0(i_fam, "_same_TAD_same_family_topTADs_vs_otherTADs_boxplot.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  boxplot(prodLogFC ~ distBin, data = tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 1,],
          xlab = my_xlab,
          ylab = my_ylab,
          outline = FALSE)
  boxplot(prodLogFC ~ distBin, data = tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 0,], col ="red", add =T, outline=FALSE)
  
  # points(jitter(tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 1,]$distBin, factor=0.5)+1,
  #        tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 1,]$prodLogFC, pch=20, cex = 0.1)
  points(jitter(match( tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 1,]$distBin,levels(tad_dist_logFC_DT$distBin))-1, factor=0.5)+1,
         tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 1,]$prodLogFC, pch=20, cex = 0.1)
  # points(jitter(tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 0,]$distBin, factor=0.5)+1,
  #        tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 0,]$prodLogFC, pch=20, col ="red", cex = 0.1)
  points(jitter(match( tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 0,]$distBin,levels(tad_dist_logFC_DT$distBin))-1, factor=0.5)+1,
         tad_dist_logFC_DT[tad_dist_logFC_DT$topTAD == 0,]$prodLogFC, pch=20, col ="red", cex = 0.1)
  
  title("logFC product between 2 genes from same TAD and same family")
  mtext(paste0(curr_dataset, " - ", i_fam, " - topTADs vs. other TADs"), side=3)
  
  legTxt <- c(paste0("same TAD + same family + topTADs (# pairs = ", sum(tad_dist_logFC_DT$topTAD == 1) , ")"), 
              paste0("same TAD + same families + other TADs (# pairs = ", sum(tad_dist_logFC_DT$topTAD == 0), ")"))
  
  legend("bottomright", legend = legTxt, text.col = c("black", "red"), bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
} # end iterating over familyData (family and family_short)

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

