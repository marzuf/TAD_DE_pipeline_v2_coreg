startTime <- Sys.time()
cat(paste0("> Rscript percent_coexpr.R\n"))

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
# distLimit <- 1000 * 10^3
# distLimitKb <- distLimit/1000
# fitMeth <- "loess"
# nbrLoessPoints <- 1000
# nSepMax <- 5

coexprThresh <- 0.5 # filter: allData_dt$coexpr > coexprThresh

plotType <- "svg"
# this is ggsave -> only inches !
myHeight <- 7
myWidth <- 9

caller <- "TopDom"

cat(paste0("!!! WARNING HARD-CODED:\n"))
cat(paste0("... corMethod =\t", corMethod, "\n"))
cat(paste0("... buildTable =\t", buildTable, "\n"))
# cat(paste0("... distLimitKb =\t",distLimitKb, "\n"))
# cat(paste0("... fitMeth =\t", fitMeth, "\n"))
# cat(paste0("... nbrLoessPoints =\t", nbrLoessPoints, "\n"))
# cat(paste0("... nSepMax =\t", nSepMax, "\n"))

### RETRIEVE FROM COMMAND LINE
# Rscript fisher_sameTADcoexpr.R <dataset> 
# Rscript fisher_sameTADcoexpr.R TCGAcrc_msi_mss
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
  paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/FISHER_SAMETAD_COEXPR_", coexprThresh),
  curr_dataset)  
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("coexpr_dist_plot_logFile_buildTable", as.character(buildTable), ".txt"))  
system(paste0("rm -f ", logFile))

txt <- paste0("> ! Hard-coded parameters:\n")
printAndLog(txt, logFile)
txt <- paste0("... corMethod = ",  corMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)
# txt <- paste0("... distLimit = ",  distLimit, "\n")
# printAndLog(txt, logFile)
# txt <- paste0("... fitMeth = ",  fitMeth, "\n")
# printAndLog(txt, logFile)


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


allData_dt$posCoexpr <- as.numeric(allData_dt$coexpr > coexprThresh)
allData_dt$posCoexpr <- factor(as.character(allData_dt$posCoexpr), levels=c("1", "0"))
allData_dt$sameTAD <- factor(as.character(allData_dt$sameTAD), levels=c("1", "0"))


fisherMatrix <- table(allData_dt$sameTAD, allData_dt$posCoexpr, dnn=c("sameTAD", "posCoexpr"))

ft <- fisher.test(fisherMatrix, alternative="greater")

resultsFT <- c(odds_ratio = as.numeric(ft$estimate), p_value = as.numeric(ft$p.value))

outFile <- file.path(outFold, "resultsFT.Rdata")
save(resultsFT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

  
######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

                                          
