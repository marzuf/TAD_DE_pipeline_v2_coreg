startTime <- Sys.time()
cat(paste0("> Rscript coexpr_dist_v2_4curves.R\n"))

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
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

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
corMethod <- "pearson"
buildTable <- FALSE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"


### RETRIEVE FROM COMMAND LINE
# Rscript coexpr_dist_v2.R
# Rscript coexpr_dist_v2.R TCGAcrc_msi_mss 50 hgnc
#  Rscript coexpr_dist_v2.R <dataset> <nTADs> <family>
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

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "COEXPR_DIST_v2_4curves",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("coexpr_dist_v2_4curves_logFile_plot.txt"))  
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
txt <- paste0("... corMethod = ",  corMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("... distLimit = ",  distLimit, "\n")
printAndLog(txt, logFile)
txt <- paste0("... fitMeth = ",  fitMeth, "\n")
printAndLog(txt, logFile)


# col1 <- "darkslateblue"
# col2 <-  "darkorange1"

mycols <- c("sameTAD" ="darkorange1" , "diffTAD"="darkslateblue",  "sameFam+sameTAD"="violetred1", "sameFam+diffTAD" = "lightskyblue")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT$midPos <- (gene2tadDT$start + gene2tadDT$end)/2


pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

toprankingScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking")
source(paste0(toprankingScriptDir, "/", "get_topTADs.R"))

utilsDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")
source(file.path(utilsDir, "coreg_utils_ggscatterhist.R"))




inFoldPairs <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_ALL_GENES"))


#============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
script0_name <- "0_prepGeneData"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline: ", length(pipeline_geneList), "\n")
printAndLog(txt, logFile)

stopifnot(pipeline_geneList %in% gene2tadDT$entrezID)

qqnormDT <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
stopifnot(names(pipeline_geneList) %in% rownames(qqnormDT))
qqnormDT <- qqnormDT[rownames(qqnormDT) %in% names(pipeline_geneList),]
stopifnot(nrow(qqnormDT) == length(pipeline_geneList))
rownames(qqnormDT) <- pipeline_geneList[rownames(qqnormDT)]
stopifnot(setequal(pipeline_geneList, rownames(qqnormDT)))

dim(qqnormDT)

cor_qqnormMat <- cor(t(qqnormDT), method = corMethod)
stopifnot(nrow(cor_qqnormMat) == length(pipeline_geneList))
stopifnot(ncol(cor_qqnormMat) == length(pipeline_geneList))


# cor_qqnormDT <- melt(cor_qqnormMat)

# TEST VERSION !!!
cor_qqnormDT <- melt(cor_qqnormMat[1:5,1:5])
colnames(cor_qqnormDT) <- c("gene1", "gene2", "coexpr")
cor_qqnormDT <- cor_qqnormDT[cor_qqnormDT$gene1 != cor_qqnormDT$gene2,]
cor_qqnormDT$gene_pair <- unlist(sapply(1:nrow(cor_qqnormDT), function(x) 
  paste0(sort(c(cor_qqnormDT$gene1[x], cor_qqnormDT$gene2[x])), collapse = "_")))


all_familyData <- paste0(familyData, c("_family", "_family_short"))






for(i_fam in all_familyData) {
  
  
  load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/PREP_ALL_GENES///hgnc_family_short_pairsDT.Rdata")
  pairsDT <- eval(parse(text = load(file.path(inFoldPairs, paste0(familyData, "_pairsDT.Rdata")))))
  pairsDT$gene1 <- as.character(pairsDT$gene1)
  pairsDT$gene2 <- as.character(pairsDT$gene2)
  # gene1     gene2 distGenes sameTAD sameFamily
  pairsDT$gene_pair <- unlist(sapply(1:nrow(pairsDT), function(x) 
    paste0(sort(c(pairsDT$gene1[x], pairsDT$gene2[x])), collapse = "_")))
  
  
  pairsDT<- pairsDT[pairsDT$gene1 %in% pipeline_geneList & pairsDT$gene2 %in% pipeline_geneList,]
  
  # TEST VERSION !!!
  pairsDT <- pairsDT[1:nrow(cor_qqnormDT),]
  pairsDT$gene1 <- cor_qqnormDT$gene1
  pairsDT$gene2 <- cor_qqnormDT$gene2
  pairsDT$gene_pair <- cor_qqnormDT$gene_pair
  
  stopifnot(pairsDT$gene_pair %in% cor_qqnormDT$gene_pair)
  stopifnot(cor_qqnormDT$gene_pair %in% pairsDT$gene_pair)
  
  fullDT <- left_join(cor_qqnormDT, pairsDT, by="gene_pair")
  
  stopifnot(fullDT$gene1.x == fullDT$gene2.x | fullDT$gene1.x == fullDT$gene2.y)
  stopifnot(fullDT$gene1.y == fullDT$gene2.x | fullDT$gene1.y == fullDT$gene2.y)
  

  fullDT$curve1 <-  ifelse(fullDT$sameTAD == "0", "diff. TAD", "same TAD")
  
  fullDT$curve2 <-  ifelse(fullDT$sameFamily == "0", NA,
                           ifelse(fullDT$sameTAD == "0", "same Fam. + diff. TAD", "same Fam. + same TAD"))
                           
  fullDT$dist_kb <- fullDT$distGenes/1000
  
  
  sameTAD_DT <- fullDT[fullDT$sameTAD == 1,c("gene_pair", "dist_kb")]
  sameTAD_DT <- na.omit(sameTAD_DT)
  sameTAD_DT <- sameTAD_DT[order(sameTAD_DT$dist_kb),]
  # sameTAD_DT$cumdist <- cumsum(sameTAD_DT$dist_kb)
  sameTAD_DT$nPair <- 1:nrow(sameTAD_DT)
  sameTAD_DT$label <- "sameTAD"
  
  diffTAD_DT <- fullDT[fullDT$sameTAD == 0,c("gene_pair", "dist_kb")]
  diffTAD_DT <- na.omit(diffTAD_DT)
  diffTAD_DT <- diffTAD_DT[order(diffTAD_DT$dist_kb),]
  # diffTAD_DT$cumdist <- cumsum(diffTAD_DT$dist_kb)
  diffTAD_DT$nPair <- 1:nrow(diffTAD_DT)
  diffTAD_DT$label <- "diffTAD"
  
  sameFam_sameTAD_DT <- fullDT[fullDT$sameFamily == 1 & fullDT$sameTAD == 1 ,c("gene_pair", "dist_kb")]
  sameFam_sameTAD_DT <- na.omit(sameFam_sameTAD_DT)
  sameFam_sameTAD_DT <- sameFam_sameTAD_DT[order(sameFam_sameTAD_DT$dist_kb),]
  # sameFam_sameTAD_DT$cumdist <- cumsum(sameFam_sameTAD_DT$dist_kb)
  sameFam_sameTAD_DT$nPair <- 1:nrow(sameFam_sameTAD_DT)
  sameFam_sameTAD_DT$label <- "sameFam+sameTAD"
  
  sameFam_diffTAD_DT <- fullDT[fullDT$sameFamily == 1 & fullDT$sameTAD == 0 ,c("gene_pair", "dist_kb")]
  sameFam_diffTAD_DT <- na.omit(sameFam_diffTAD_DT)
  sameFam_diffTAD_DT <- sameFam_diffTAD_DT[order(sameFam_diffTAD_DT$dist_kb),]
  # sameFam_diffTAD_DT$cumdist <- cumsum(sameFam_diffTAD_DT$dist_kb)
  sameFam_diffTAD_DT$nPair <- 1:nrow(sameFam_diffTAD_DT)
  sameFam_diffTAD_DT$label <- "sameFam+diffTAD"
  
  cumsum_4_curves_DT <- rbind(rbind(sameTAD_DT, diffTAD_DT),
                      rbind(sameFam_sameTAD_DT,sameFam_diffTAD_DT))
  
  cumsum_4_curves_DT <- cumsum_4_curves_DT[cumsum_4_curves_DT$dist_kb <= distLimit/1000,] # distLimit was given in bp !
  
  cumsum_4_curves_DT$label <- factor(cumsum_4_curves_DT$label, 
                                     levels=names(mycols))
  
  
  ggplot(cumsum_4_curves_DT, aes(x=dist_kb, y = nPair, color=label))+
    geom_smooth(se=FALSE, method = fitMeth)+
    scale_color_manual(values=mycols)
  
  
  ggplot() +
    geom_smooth(data=fullDT, aes(x=dist_kb, y=coexpr, color=curve1)) +
    geom_smooth(data=fullDT[fullDT$sameFamily==1,], aes(x=dist_kb, y=coexpr, color=curve2)) + 
    mytheme
  
  
  
  
 
 outFile <- file.path(outFold, paste0(i_fam, "_same_family_sameTADtopTADs_vs_sameTADotherTADs_vs_diffTADs.", plotType))
 ggsave(outFile, height = myHeight, width = myWidth)
 cat(paste0("... written: ", outFile, "\n"))
 
    
} # end iterating over familyData (family and family_short)

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

