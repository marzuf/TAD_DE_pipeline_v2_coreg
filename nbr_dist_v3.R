startTime <- Sys.time()
cat(paste0("> Rscript nbr_dist_v3.R\n"))

options(scipen=100)

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggstatsplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


axisLabSize <- 12
legendSize <- 10
plotTitSize <- 14

mytheme <- theme(
  # top, right, bottom and left
  plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
  plot.title = element_text(hjust = 0.5, face = "bold", size=plotTitSize, vjust=1),
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

### HARD CODED
caller <- "TopDom"

buildTable <- FALSE
# for plotting:
# look at nbr ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"


### RETRIEVE FROM COMMAND LINE
# Rscript nbr_dist_v3.R
# Rscript nbr_dist_v3.R hgnc
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 1)
  
  familyData <- args[1]
} else{
  txt <- paste0("> Default parameters:\n")
  familyData <- "hgnc"
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "NBR_DIST_v3",  paste0(familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("nbr_dist_v3_logFile.txt"))  
system(paste0("rm -f ", logFile))

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
txt <- paste0("... fitMeth = ",  fitMeth, "\n")
printAndLog(txt, logFile)

# mycols <- c("sameTAD" ="darkorange1" , "diffTAD"="darkslateblue",  "sameFam+sameTAD"="violetred1", "sameFam+diffTAD" = "lightskyblue")
mycols <- c("same TAD" ="darkorange1" , "diff. TAD"="darkslateblue",  "same Fam. + same TAD"="violetred1", "same Fam. + diff. TAD" = "lightskyblue")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

toprankingScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking")
source(paste0(toprankingScriptDir, "/", "get_topTADs.R"))

utilsDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")
source(file.path(utilsDir, "coreg_utils_ggscatterhist.R"))

################################################ DATA PREPARATION

if(buildTable) {
  
  cat(paste0("... load DIST data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_DIST/all_dist_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_dist_pairs)
  nrow(all_dist_pairs)
  all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
  all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
  
  cat(paste0("... load TAD data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_TAD/all_TAD_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_TAD_pairs)
  nrow(all_TAD_pairs)
  all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
  all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
  

  
  #============================== THIS VERSION IS DATASET INDEPENDENT -> DO NOT FILTER PIPELINE GENE LIST 

  dist_pair <- all_dist_pairs
  dist_pairs_limit <- dist_pair[dist_pair$dist <= distLimit,]
  head(dist_pairs_limit)
  nrow(dist_pairs_limit)
  
  TAD_pairs <- all_TAD_pairs
  head(TAD_pairs)
  nrow(TAD_pairs)
  
  
  # START MERGING DATA 
  
  cat(paste0("... merge DIST - TAD data\t", Sys.time(), "\t"))
  dist_TAD_DT <- left_join(dist_pairs_limit, TAD_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  dist_TAD_DT$sameTAD <- ifelse(is.na(dist_TAD_DT$region), 0, 1)
  
}

all_familyData <- paste0(familyData, c("_family", "_family_short"))

for(i_fam in all_familyData) {
  
  if(buildTable){
    
    
    cat(paste0("... load FAMILY data\t", Sys.time(), "\t"))
    load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_FAMILY/", i_fam, "_all_family_pairs.Rdata")))
    cat(paste0(Sys.time(), "\n"))
    head(all_family_pairs)
    nrow(all_family_pairs)
    all_family_pairs$gene1 <- as.character(all_family_pairs$gene1)
    all_family_pairs$gene2 <- as.character(all_family_pairs$gene2)
    
    family_pairs <- all_family_pairs
    head(family_pairs)
    
    cat(paste0("... merge FAMILY data\t", Sys.time(), "\t"))
    
    dist_TAD_fam_DT <- left_join(dist_TAD_DT, family_pairs, by=c("gene1", "gene2"))
    cat(paste0(Sys.time(), "\n"))
    
    dist_TAD_fam_DT$sameFamily <- ifelse(is.na(dist_TAD_fam_DT$family), 0, 1)
    
    # DATASET INDEPENDENT ! NO COEXPR
    # cat(paste0("... merge COEXPR data\t", Sys.time(), "\t"))
    # dist_TAD_fam_coexpr_DT <- left_join(dist_TAD_fam_DT, coexprDT, by=c("gene1", "gene2"))
    # cat(paste0(Sys.time(), "\n"))
    
    allData_dt <- dist_TAD_fam_DT
    allData_dt$region <- NULL
    allData_dt$family <- NULL
    allData_dt <- na.omit(allData_dt)
    
    outFile <-file.path(outFold, paste0(i_fam, "_", "allData_dt.Rdata"))
    save(allData_dt, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))    
    
  } else{
    outFile <-file.path(outFold, paste0(i_fam, "_", "allData_dt.Rdata"))
    load(outFile)
    # load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3/TCGAcrc_msi_mss_hgnc/hgnc_family_allData_dt.Rdata")
    
  }
  
  nrow(allData_dt)
  allData_dt$dist_kb <- allData_dt$dist/1000
  
  allData_dt$curve1 <-  ifelse(allData_dt$sameTAD == "0", "diff. TAD", "same TAD")
  
  allData_dt$curve2 <-  ifelse(allData_dt$sameFamily == "0", NA,
                               ifelse(allData_dt$sameTAD == "0", "same Fam. + diff. TAD", "same Fam. + same TAD"))
  
  sameTAD_DT <- allData_dt[allData_dt$sameTAD == 1,c("gene1", "gene2", "dist_kb")]
  sameTAD_DT <- na.omit(sameTAD_DT)
  sameTAD_DT <- sameTAD_DT[order(sameTAD_DT$dist_kb),]
  # sameTAD_DT$cumdist <- cumsum(sameTAD_DT$dist_kb)
  sameTAD_DT$nPair <- 1:nrow(sameTAD_DT)
  sameTAD_DT$label <- "same TAD"
  
  diffTAD_DT <- allData_dt[allData_dt$sameTAD == 0,c("gene1", "gene2",  "dist_kb")]
  diffTAD_DT <- na.omit(diffTAD_DT)
  diffTAD_DT <- diffTAD_DT[order(diffTAD_DT$dist_kb),]
  # diffTAD_DT$cumdist <- cumsum(diffTAD_DT$dist_kb)
  diffTAD_DT$nPair <- 1:nrow(diffTAD_DT)
  diffTAD_DT$label <- "diff. TAD"
  
  sameFam_sameTAD_DT <- allData_dt[allData_dt$sameFamily == 1 & allData_dt$sameTAD == 1 ,c("gene1", "gene2", "dist_kb")]
  sameFam_sameTAD_DT <- na.omit(sameFam_sameTAD_DT)
  sameFam_sameTAD_DT <- sameFam_sameTAD_DT[order(sameFam_sameTAD_DT$dist_kb),]
  # sameFam_sameTAD_DT$cumdist <- cumsum(sameFam_sameTAD_DT$dist_kb)
  sameFam_sameTAD_DT$nPair <- 1:nrow(sameFam_sameTAD_DT)
  sameFam_sameTAD_DT$label <- "same Fam. + same TAD"
  
  sameFam_diffTAD_DT <- allData_dt[allData_dt$sameFamily == 1 & allData_dt$sameTAD == 0 ,c("gene1", "gene2","dist_kb")]
  sameFam_diffTAD_DT <- na.omit(sameFam_diffTAD_DT)
  sameFam_diffTAD_DT <- sameFam_diffTAD_DT[order(sameFam_diffTAD_DT$dist_kb),]
  # sameFam_diffTAD_DT$cumdist <- cumsum(sameFam_diffTAD_DT$dist_kb)
  sameFam_diffTAD_DT$nPair <- 1:nrow(sameFam_diffTAD_DT)
  sameFam_diffTAD_DT$label <- "same Fam. + diff. TAD"
  
  sum_4_curves_DT <- rbind(rbind(sameTAD_DT, diffTAD_DT),
                           rbind(sameFam_sameTAD_DT,sameFam_diffTAD_DT))
  
  stopifnot(sum_4_curves_DT$dist_kb <= distLimit/1000)
  
  sum_4_curves_DT$label <- factor(sum_4_curves_DT$label,
                                  levels=names(mycols))
  
  stopifnot(!any(is.na(sum_4_curves_DT$label)))
  
  my_xlab <- paste0("Distance between the 2 genes (kb)")
  my_ylab <- "Nbr of gene pairs"
  my_sub <- paste0("(", i_fam, " - all data)")
  
  p_nbr <- ggplot(sum_4_curves_DT, aes(x=dist_kb, y = nPair, color=label))+
    ggtitle("Number of pairs along distance betw. genes", subtitle = my_sub) + 
    xlab(my_xlab)+
    ylab(my_ylab)+
    geom_smooth(se=FALSE, method = fitMeth)+
    scale_color_manual(values=mycols, name="") + 
    mytheme
  
  outFile <- file.path(outFold, paste0(i_fam, "_sameTAD_diffTAD_sameTADfam_diffTADfam_4curves_nbr_dist.", plotType))
  ggsave(p_nbr, height = myHeight, width = myWidth, filename=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  p_nbr1 <- my_ggscatterhist(
    sum_4_curves_DT, 
    x = "dist_kb", 
    y = "nPair",
    xlab=my_xlab,
    ylab=my_ylab,
    title = paste0("Number of pairs along distance betw. genes"),
    subtitle=my_sub,
    legend.title="",
    point = FALSE,
    # rug=TRUE,
    add = fitMeth,
    color = "label", size = 3, alpha = 0.6,
    palette = mycols,
    x_margin.params = list(fill = "label", color = "black", size = 0.2),
    x_margin.ggtheme = theme_minimal(),
    x_margin.plot = "density",
    y_margin.params = list(fill = "label", color = "black", size = 0.2),
    y_margin.ggtheme = theme_minimal(),
    y_margin.plot = "boxplot",
    plot_xmargin = TRUE,
    plot_ymargin=TRUE,
    ymargin_as_xmargin = FALSE,
    # global_margins_cm = c(0, 0, 0.25, 0.25)
    global_margins_cm = c(0.25, 0.25, 0.25, 0.25)
  )
  outFile <- file.path(outFold, paste0(i_fam, "_sameTAD_diffTAD_sameTADfam_diffTADfam_4curves_nbr1_dist.", plotType))
  ggsave(p_nbr1, height = myHeight*1.2, width = myWidth*1.2, filename=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  

  
} # end iterating over familyData (family and family_short)

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

