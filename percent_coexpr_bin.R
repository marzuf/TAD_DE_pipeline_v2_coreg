startTime <- Sys.time()
cat(paste0("> Rscript percent_coexpr_bin.R\n"))

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggstatsplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(tools, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

axisLabSize <- 12
legendSize <- 10
plotTitSize <- 14

mytheme <- theme(
  # top, right, bottom and left
  plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
  plot.title = element_text(hjust = 0.5, face = "bold", size=plotTitSize, vjust=1),
  plot.subtitle = element_text(hjust = 0.5, face = "italic", size=plotTitSize-2, vjust=1),
  panel.background = element_rect(fill = "white", colour = NA), 
  panel.border = element_rect(fill = NA, colour = "grey20"), 
  panel.grid.major = element_line(colour = "grey92"), 
  panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
  strip.background = element_rect(fill = "grey85", colour = "grey20"), 
  #legend.key = element_rect(fill = "white", colour = NA), 
  axis.line.x = element_line(size = .3, color = "black"),
  axis.line.y = element_line(size = .3, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=axisLabSize),
  axis.text.x = element_text(color="black", hjust=0.5,vjust = 1, size=axisLabSize),
  axis.title.y = element_text(color="black", size=axisLabSize+1),
  axis.title.x = element_text(color="black", size=axisLabSize+1),
  legend.text = element_text(size=legendSize),
  legend.key.height = unit(1.5,"cm"),
  legend.key = element_blank()
)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/revision_utils.R"))

### HARD CODED
corMethod <- "pearson"
buildTable <- FALSE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 1000 * 10^3
# distLimitKb <- distLimit/1000
fitMeth <- "loess"
nbrLoessPoints <- 1000
nBinMax <- 5

# in bp
binSepSize <- 50000

coexprThresh <- 0.5

plotType <- "svg"
# this is ggsave -> only inches !
myHeight <- 7
myWidth <- 9

caller <- "TopDom"

cat(paste0("!!! WARNING HARD-CODED:\n"))
cat(paste0("... corMethod =\t", corMethod, "\n"))
cat(paste0("... buildTable =\t", buildTable, "\n"))
# cat(paste0("... distLimitKb =\t",distLimitKb, "\n"))
cat(paste0("... fitMeth =\t", fitMeth, "\n"))
cat(paste0("... nbrLoessPoints =\t", nbrLoessPoints, "\n"))
cat(paste0("... nBinMax =\t", nBinMax, "\n"))
cat(paste0("... coexprThresh =\t", coexprThresh, "\n"))

### RETRIEVE FROM COMMAND LINE
# Rscript percent_coexpr_bin.R <dataset> 
# Rscript percent_coexpr_bin.R TCGAcrc_msi_mss
curr_dataset <- "TCGAcrc_msi_mss"
args <- commandArgs(trailingOnly = T)
stopifnot(length(args) == 1)
curr_dataset <- args[1]

noBuildTableInFold <- file.path(setDir, 
                                paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/AUC_COEXPRDIST_SORTNODUP"),
                                curr_dataset)
stopifnot(file.exists(noBuildTableInFold))

################# UPDATE 30.08.2018 
# => USE SORTED NOT DUPLICATED TABLE
# => CHECK GENE1 < GENE2
# => DO NOT COMPUTE LENGTH(UNIQUE(C())) THAT TAKES FOREVER

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
stopifnot(file.exists(entrezDT_file))

outFold <- file.path(
  setDir, 
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/PERCENT_COEXPR_BIN",
  paste0("PERCENT_COEXPR_", distLimit/1000, "kb_", nBinMax, "BinsSep_", coexprThresh), curr_dataset)  
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("coexpr_dist_plot_logFile_buildTable", as.character(buildTable), ".txt"))  
system(paste0("rm -f ", logFile))

txt <- paste0("> ! Hard-coded parameters:\n")
printAndLog(txt, logFile)
txt <- paste0("... corMethod = ",  corMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("... distLimit = ",  distLimit, "\n")
printAndLog(txt, logFile)
txt <- paste0("... fitMeth = ",  fitMeth, "\n")
printAndLog(txt, logFile)
txt <- paste0("\n")
printAndLog(txt, logFile)

# mycols <- c("sameTAD" ="darkorange1" , "diffTAD"="darkslateblue",  "sameFam+sameTAD"="violetred1", "sameFam+diffTAD" = "lightskyblue")
mycols <- c("same TAD" ="darkorange1" , "diff. TAD"="darkslateblue",  "same Fam. + same TAD"="violetred1", "same Fam. + diff. TAD" = "lightskyblue")

sameTADcol_daniele <- "chocolate2"
diffTADcol_daniele <- "gray36"

################################################ DATA PREPARATION

if(buildTable) {
  
  cat(paste0("... load DIST data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_DIST_SORTNODUP/all_dist_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_dist_pairs)
  nrow(all_dist_pairs)
  all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
  all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
  ### UPDATE 30.06.2018
  stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
  
  cat(paste0("... load TAD data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_TAD_SORTNODUP/all_TAD_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_TAD_pairs)
  nrow(all_TAD_pairs)
  all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
  all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
  ### UPDATE 30.06.2018
  stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
  
  cat(paste0("... load COEXPR data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_COEXPR_SORTNODUP"),  paste0(curr_dataset, "_", corMethod), "coexprDT.Rdata"))
  cat(paste0(Sys.time(), "\n"))
  head(coexprDT)
  nrow(coexprDT)
  coexprDT$gene1 <- as.character(coexprDT$gene1)
  coexprDT$gene2 <- as.character(coexprDT$gene2)
  ### UPDATE 30.06.2018
  stopifnot(coexprDT$gene1 < coexprDT$gene2)
  
  #============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
  dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)
  
  script0_name <- "0_prepGeneData"
  pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
  pipeline_geneList <- as.character(pipeline_geneList)
  
  
  dataset_dist_pair <- all_dist_pairs[all_dist_pairs$gene1 %in% pipeline_geneList & 
                                        all_dist_pairs$gene2 %in% pipeline_geneList,]
  
  dataset_dist_pairs_limit <- dataset_dist_pair[dataset_dist_pair$dist <= distLimit,]
  head(dataset_dist_pairs_limit)
  nrow(dataset_dist_pairs_limit)
  
  dataset_TAD_pairs <- all_TAD_pairs[all_TAD_pairs$gene1 %in% pipeline_geneList & 
                                       all_TAD_pairs$gene2 %in% pipeline_geneList,]
  
  # START MERGING DATA 
  
  cat(paste0("... merge DIST - TAD data\t", Sys.time(), "\t"))
  dataset_dist_TAD_DT <- left_join(dataset_dist_pairs_limit, dataset_TAD_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  dataset_dist_TAD_DT$sameTAD <- ifelse(is.na(dataset_dist_TAD_DT$region), 0, 1)
  
  cat(paste0("... merge COEXPR data\t", Sys.time(), "\t"))
  
  dataset_dist_TAD_coexpr_DT <- left_join(dataset_dist_TAD_DT, coexprDT, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  allData_dt <- dataset_dist_TAD_coexpr_DT
  allData_dt$region <- NULL
  allData_dt <- na.omit(allData_dt)
  
  outFile <-file.path(outFold, paste0("allData_dt.Rdata"))
  save(allData_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))    
}  else{
  outFile <-file.path(noBuildTableInFold, paste0("allData_dt.Rdata"))
  allData_dt <- eval(parse(text = load(outFile)))
  # load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3/TCGAcrc_msi_mss_hgnc/hgnc_family_allData_dt.Rdata")
}

# load gene coordinate 
geneCoordDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F) 
head(geneCoordDT)
geneCoordDT$entrezID <- as.character(geneCoordDT$entrezID)
stopifnot(is.numeric(geneCoordDT$start))
stopifnot(is.numeric(geneCoordDT$end))
stopifnot(!any(duplicated(geneCoordDT$symbol)))
stopifnot(!any(duplicated(geneCoordDT$entrezID)))
geneCoordDT <- geneCoordDT[order(geneCoordDT$chromo, geneCoordDT$start, geneCoordDT$end),]
geneCoordDT$nGene <- 1:nrow(geneCoordDT)

# if binSepSize = 500
# 0 if 0-499; 1 if 500-999
allData_dt$nBinSep <- floor(allData_dt$dist/binSepSize)

allData_dt_withPos <- left_join(allData_dt, geneCoordDT[,c("entrezID", "nGene")], by=c("gene1" ="entrezID"))
colnames(allData_dt_withPos)[colnames(allData_dt_withPos) == "nGene"] <- "nGene1"

allData_dt_withPos <- left_join(allData_dt_withPos, geneCoordDT[,c("entrezID", "nGene")], by=c("gene2" ="entrezID"))
colnames(allData_dt_withPos)[colnames(allData_dt_withPos) == "nGene"] <- "nGene2"

allData_dt_withPos$nBinSep <- abs(allData_dt_withPos$nGene2 - allData_dt_withPos$nGene1) - 1 


# plot only the coexpressed -> positive coexpression value
plotDT <- allData_dt_withPos

plotDT_all <- allData_dt_withPos[allData_dt_withPos$nBinSep <= nBinMax,]
plotDT_all <- aggregate(sameTAD ~ nBinSep, FUN=length, data=plotDT_all)
colnames(plotDT_all)[2] <- "nGenes"
plotDT_all
write.table(plotDT_all, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_all_posCorr <- allData_dt_withPos[allData_dt_withPos$nBinSep <= nBinMax,]
plotDT_all_posCorr_agg <- aggregate(coexpr ~ nBinSep, FUN=function(x) sum(as.numeric(x > coexprThresh)), data=plotDT_all_posCorr)
colnames(plotDT_all_posCorr_agg)[2] <- "nPosCoexprAll"
plotDT_all_posCorr_agg
write.table(plotDT_all_posCorr_agg, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_diffTAD_all <- allData_dt_withPos[allData_dt_withPos$sameTAD == 0 & allData_dt_withPos$nBinSep <= nBinMax,]
plotDT_diffTAD_all <- aggregate(sameTAD ~ nBinSep, FUN=length, data=plotDT_diffTAD_all)
colnames(plotDT_diffTAD_all)[2] <- "nGenesDiffTAD"
plotDT_diffTAD_all
write.table(plotDT_diffTAD_all, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_diffTAD <- allData_dt_withPos[allData_dt_withPos$sameTAD == 0  & allData_dt_withPos$nBinSep <= nBinMax,]
plotDT_diffTAD_posCorr_agg <- aggregate(coexpr ~ nBinSep, data=plotDT_diffTAD, FUN=function(x) sum(as.numeric(x > coexprThresh)))
colnames(plotDT_diffTAD_posCorr_agg)[2] <- "nPosCoexprDiffTAD"
plotDT_diffTAD_posCorr_agg
write.table(plotDT_diffTAD_posCorr_agg, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_diffTAD_negCorr_agg <- aggregate(coexpr ~ nBinSep, data=plotDT_diffTAD, FUN=function(x) sum(as.numeric(x <= coexprThresh)))
colnames(plotDT_diffTAD_negCorr_agg)[2] <- "nNegCoexprDiffTAD"
plotDT_diffTAD_negCorr_agg
write.table(plotDT_diffTAD_negCorr_agg, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_sameTAD <- allData_dt_withPos[allData_dt_withPos$sameTAD == 1 & allData_dt_withPos$nBinSep <= nBinMax,]
plotDT_sameTAD_posCorr_agg <- aggregate(coexpr ~ nBinSep, data=plotDT_sameTAD, FUN=function(x) sum(as.numeric(x > coexprThresh)))
colnames(plotDT_sameTAD_posCorr_agg)[2] <- "nPosCoexprSameTAD"
plotDT_sameTAD_posCorr_agg
write.table(plotDT_sameTAD_posCorr_agg, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_sameTAD_negCorr_agg <- aggregate(coexpr ~ nBinSep, data=plotDT_sameTAD, FUN=function(x) sum(as.numeric(x <= coexprThresh)))
colnames(plotDT_sameTAD_negCorr_agg)[2] <- "nNegCoexprSameTAD"
plotDT_sameTAD_negCorr_agg
write.table(plotDT_sameTAD_negCorr_agg, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

plotDT_sameTAD_all_agg <- aggregate(coexpr ~ nBinSep, data=plotDT_sameTAD, FUN=length)
colnames(plotDT_sameTAD_all_agg)[2] <- "nGenesSameTAD"
plotDT_sameTAD_all_agg
write.table(plotDT_sameTAD_all_agg, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

stopifnot(plotDT_sameTAD_negCorr_agg$nNegCoexpr + plotDT_sameTAD_posCorr_agg$nPosCoexpr == plotDT_sameTAD_all_agg$nAll)


plotDT_final <- merge(merge(merge(merge(merge(merge(plotDT_all, plotDT_diffTAD_all, by="nBinSep"), plotDT_sameTAD_all_agg, by="nBinSep"),
                plotDT_sameTAD_negCorr_agg, by="nBinSep"), plotDT_sameTAD_posCorr_agg, by ="nBinSep"),
                plotDT_diffTAD_negCorr_agg, by ="nBinSep"), plotDT_diffTAD_posCorr_agg, by = "nBinSep")
write.table(plotDT_final, row.names=FALSE, col.names=TRUE, file = logFile, append=TRUE, sep="\t")
txt <- paste0("\n")
printAndLog(txt, logFile)

stopifnot( plotDT_final$nGenesSameTAD + plotDT_final$nGenesDiffTAD == plotDT_final$nGenes)
stopifnot( plotDT_final$nNegCoexprSameTAD + plotDT_final$nPosCoexprSameTAD == plotDT_final$nGenesSameTAD)
stopifnot( plotDT_final$nNegCoexprDiffTAD + plotDT_final$nPosCoexprDiffTAD == plotDT_final$nGenesDiffTAD)

stopifnot( plotDT_final$nPosCoexprDiffTAD + plotDT_final$nPosCoexprSameTAD == plotDT_all_posCorr_agg$nPosCoexprAll)


# at a gene sep of 1 
# red bar = among all the correlated genes, how many % are in same TAD ??
# grey bar = among all genes, how many % are in same TAD ?

plotDT_final$redBars <- 100 * plotDT_final$nPosCoexprSameTAD/(plotDT_final$nPosCoexprSameTAD + plotDT_final$nPosCoexprDiffTAD)

plotDT_final$greyBars <- 100 * plotDT_final$nGenesSameTAD/(plotDT_final$nGenesSameTAD + plotDT_final$nGenesDiffTAD)

plotDT_final$nBinSep <- factor(plotDT_final$nBinSep, levels=as.character(0:nBinMax))

plotDT_m <- melt(plotDT_final[,c("nBinSep", "redBars", "greyBars")], id="nBinSep")

corrBarPlot <- ggplot(plotDT_m, aes(x=nBinSep, y=value, fill=variable)) +
  ggtitle("Percentage of correlated gene pairs", 
          subtitle = paste0(curr_dataset, " - binSepSize = ", binSepSize/1000, "kb - coexpr. ", coexprThresh))+
  geom_bar(stat="identity", position="dodge")+
  scale_y_continuous(name = "% of genes in same TAD")+
  # scale_x_continuous(name = "# of genes separation", breaks = c(0:nBinMax), labels=as.character(c(0:nBinMax)))+
  scale_x_discrete(name = "# of bins separation", labels=c(0:nBinMax))+
  scale_fill_manual(name="",
                    labels = c( "greyBars"="all", "redBars" = "corr."),
                    values=c("redBars" = "brown2", "greyBars" = "darkgray"))+
  # labs(fill="")+
  mytheme

outFile <- file.path(outFold, paste0("barplot_percentGenesSameTAD.", plotType))
ggsave(filename = outFile, plot=corrBarPlot, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))


  
######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

                                          
