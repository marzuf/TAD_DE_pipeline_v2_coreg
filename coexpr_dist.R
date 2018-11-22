startTime <- Sys.time()
cat(paste0("> Rscript coexpr_dist.R\n"))

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

### HARD CODED
caller <- "TopDom"
corMethod <- "pearson"
buildTable <- TRUE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 2 * 10^6

### RETRIEVE FROM COMMAND LINE
# Rscript coexpr_dist.R
# Rscript coexpr_dist.R TCGAcrc_msi_mss 50 hgnc
#  Rscript coexpr_dist.R <dataset> <nTADs> <family>
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

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "COEXPR_DIST",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("coexpr_dist_logFile.txt"))  
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


col1 <- "red"
col2 <- "green"
col3 <- "blue"

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

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))

familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)
# for all families
# -> for all gene pairs
#  --> distance
#  --> coexpression
#  --> same TAD ?

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

# take only genes from the pipeline
ds_familyDT <- familyDT[familyDT$entrezID %in% pipeline_geneList,]
# just ensure that I kept only the "TRUE" family data (not the single "unassigned" family)
stopifnot(!any(grepl("assigned", ds_familyDT$hgnc_family_short)))

txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline with family info.: ", nrow(ds_familyDT), "\n")
printAndLog(txt, logFile)

all_familyData <- paste0(familyData, c("_family", "_family_short"))

for(i_fam in all_familyData) {
  
  if(buildTable) {
    
    cat(paste0("> START: ", i_fam, "\n"))
    
    all_families <- unique(as.character(ds_familyDT[, i_fam]))
    cat(paste0("... # unique families: ", length(all_families), "\n"))
    
    family_dist_coexpr_DT <- foreach(curr_family = all_families, .combine='rbind') %dopar% {
      
      fam_genes <- ds_familyDT$entrezID[ds_familyDT[,i_fam] == curr_family]  
      stopifnot(fam_genes %in% pipeline_geneList)
      
      if(length(fam_genes) == 1) {
        return(data.frame(
          gene1 = fam_genes,
          gene2 = fam_genes,
          dist = NA,
          coexpr = NA,
          sameTAD = NA,
          family = curr_family,
          stringsAsFactors = FALSE
        ))
      }
      
      # create combn of gene pairs
      genePairs <- combn(x = fam_genes, m = 2)
      
      stopifnot(ncol(genePairs) == (length(fam_genes) * (length(fam_genes) -1) *0.5))
      
      genepairs_DT <- foreach(i_pair = 1:ncol(genePairs), .combine='rbind') %do% {
        gene1 <- genePairs[1,i_pair]
        gene2 <- genePairs[2,i_pair]
        stopifnot(gene1 %in% rownames(qqnormDT), gene1 %in% gene2tadDT$entrezID, 
                  gene2 %in% rownames(qqnormDT), gene2 %in% gene2tadDT$entrezID)
        chromo1 <- gene2tadDT$chromo[gene2tadDT$entrezID == gene1]
        chromo2 <- gene2tadDT$chromo[gene2tadDT$entrezID == gene2]
        stopifnot(length(chromo1) == length(chromo2))
        stopifnot(length(chromo1) == 1)
        pair_dist <- ifelse(chromo1 == chromo2, 
                            abs(gene2tadDT$midPos[gene2tadDT$entrezID == gene1] - gene2tadDT$midPos[gene2tadDT$entrezID == gene2]), NA)
        
        pair_coexpr <- cor(qqnormDT[gene1, ], qqnormDT[gene2,], method = corMethod)
        
        pair_sameTAD <- as.numeric(gene2tadDT$region[gene2tadDT$entrezID == gene1] == gene2tadDT$region[gene2tadDT$entrezID == gene2])
        
        data.frame(
          gene1 = gene1,
          gene2 = gene2,
          dist = pair_dist,
          coexpr = pair_coexpr,
          sameTAD = pair_sameTAD,
          family = curr_family,
          stringsAsFactors = FALSE
        )
      } # end iterating over gene pairs
      genepairs_DT
    } # end iterating over the families
    outFile <- file.path(outFold, paste0(i_fam, "_family_dist_coexpr_DT.Rdata"))
    save(family_dist_coexpr_DT, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))

  } else {
    # data already built in previous run
    outFile <- file.path(outFold, paste0(i_fam, "_family_dist_coexpr_DT.Rdata"))
    # outFile <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST/TCGAcrc_msi_mss/hgnc_family_family_dist_coexpr_DT.Rdata"
    family_dist_coexpr_DT <- eval(parse(text = load(outFile)))
  }
  txt <- paste0("... # of gene pairs: ", nrow(family_dist_coexpr_DT), "\n")
  printAndLog(txt, logFile)
  txt <- paste0("... # of families: ", length(unique(family_dist_coexpr_DT$family)), "\n")
  printAndLog(txt, logFile)
  txt <- paste0("... family with 1 gene: ", sum(family_dist_coexpr_DT$gene1 == family_dist_coexpr_DT$gene2), "\n")
  printAndLog(txt, logFile)
  
  family_dist_coexpr_DT <- family_dist_coexpr_DT[family_dist_coexpr_DT$gene1 != family_dist_coexpr_DT$gene2,]
  txt <- paste0("... # of gene pairs (diff. genes only): ", nrow(family_dist_coexpr_DT), "\n")
  printAndLog(txt, logFile)
  
  family_dist_coexpr_DT <- na.omit( family_dist_coexpr_DT)
  txt <- paste0("... # of gene pairs (same chromosome only): ", nrow(family_dist_coexpr_DT), "\n")
  printAndLog(txt, logFile)
  
  family_dist_coexpr_DT <- family_dist_coexpr_DT[family_dist_coexpr_DT$dist <= distLimit,]
  txt <- paste0("... # of gene pairs within ", distLimit, " bp: ", nrow(family_dist_coexpr_DT), "\n")
  printAndLog(txt, logFile)
  
  family_dist_coexpr_DT$tad1 <- unlist(sapply(family_dist_coexpr_DT$gene1, function(x)
    gene2tadDT$region[gene2tadDT$entrezID == x]))
  
  family_dist_coexpr_DT$tad2 <- unlist(sapply(family_dist_coexpr_DT$gene2, function(x)
    gene2tadDT$region[gene2tadDT$entrezID == x]))
  
  family_dist_coexpr_DT$same_tad <- as.numeric(family_dist_coexpr_DT$tad1 == family_dist_coexpr_DT$tad2)
  stopifnot(family_dist_coexpr_DT$same_tad == family_dist_coexpr_DT$sameTAD)
  
  topTADs <- get_topTADs(mydataset = curr_dataset, nTop = nTopTADs, TADonly = TRUE)
  family_dist_coexpr_DT$topTAD1 <- unlist(sapply(family_dist_coexpr_DT$tad1, function(x) as.numeric(x %in% topTADs )))
  family_dist_coexpr_DT$topTAD2 <- unlist(sapply(family_dist_coexpr_DT$tad2, function(x) as.numeric(x %in% topTADs )))
  
  
  ##################################### PLOT1: BOXPLOT COEXPRESSION 
  ##################################### SAME FAMILY: SAME TADs VS DIFF TADs
  
  my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
  
  outFile <- file.path(outFold, paste0(i_fam, "_boxplot_coexpr_sameTAD_diffTAD.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  boxplot(coexpr ~ sameTAD, data = family_dist_coexpr_DT, 
          ylab = my_ylab,
          xlab="sameTAD")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  dataframeDT <- family_dist_coexpr_DT
  dataframeDT$sameTAD <- factor(as.character(dataframeDT$sameTAD), levels = c("0", "1"))
  
  my_comparisons <- list( c("0", "1") )
  
  violinP <- ggviolin(dataframeDT, x = "sameTAD", y = "coexpr", fill = "sameTAD",
                      legend.title = "",
                      # yscale = "log10",
                      # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                      palette = c("#00AFBB", "#FC4E07"),
                      xlab ="",
                      ylab =paste0(corMethod, "'s corr. coexpr."),
                      add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
    # stat_compare_means(label.y.npc="top", label.x.npc = "left")    # Add global the p-value 
    stat_compare_means(label.y = max(dataframeDT$coexpr), label.x = 0.5)
  
  violinP <- violinP +
    geom_hline(yintercept = mean(dataframeDT$coexpr[dataframeDT$type=="0"]), linetype=2, color = "#00AFBB", size=1)+
    geom_hline(yintercept = mean(dataframeDT$coexpr[dataframeDT$type=="1"]), linetype=2, color = "#FC4E07", size=1)
  
  if(SSHFS) violinP
  outFile <- file.path(outFold, paste0(i_fam, "_cmp_violin_sameTAD_diffTAD.svg"))
  ggsave(filename = outFile, plot = violinP, height = 7, width=8)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ##################################### PLOT2: COEXPRESSION ~ DISTANCE
  ##################################### SAME FAMILY: SAME TADs VS DIFF TADs
  
  my_xlab <- paste0("Distance between the 2 genes (bp)")
                    
  outFile <- file.path(outFold, paste0(i_fam, "_same_family_sameTADtopTADs_vs_sameTADotherTADs_vs_diffTADs.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  plot(x = family_dist_coexpr_DT$dist,
        # x = log10(family_dist_coexpr_DT$dist),
        y = family_dist_coexpr_DT$coexpr,
       ylab = my_ylab,
       xlab = my_xlab,
       cex = 0.5,
       pch = 16,
       main = "Coexpr~dist for gene pairs from same family")
  mtext(paste0(curr_dataset, " - ", i_fam, " - same TAD topTAD / same TAD otherTADs / diff. TADs"), side = 3)
  lines(smooth.spline(x = family_dist_coexpr_DT$dist[family_dist_coexpr_DT$same_tad == 1 & family_dist_coexpr_DT$topTAD1 == 1 ],
                     y = family_dist_coexpr_DT$coexpr[family_dist_coexpr_DT$same_tad == 1  & family_dist_coexpr_DT$topTAD1 == 1]), 
       col =col1)
  
  lines(smooth.spline(x = family_dist_coexpr_DT$dist[family_dist_coexpr_DT$same_tad == 1 & family_dist_coexpr_DT$topTAD1 == 0 ],
                      y = family_dist_coexpr_DT$coexpr[family_dist_coexpr_DT$same_tad == 1  & family_dist_coexpr_DT$topTAD1 == 0]), 
        col =col2)
  
  
  lines(smooth.spline(x = family_dist_coexpr_DT$dist[family_dist_coexpr_DT$same_tad == 0],
                      y = family_dist_coexpr_DT$coexpr[family_dist_coexpr_DT$same_tad == 0]), 
        col =col3)
  
  legTxt <- c(paste0("same TAD topTADs (# pairs = ", sum(family_dist_coexpr_DT$same_tad == 1  & family_dist_coexpr_DT$topTAD1 == 1 ), ")"),
              paste0("same TAD otherTADs (# pairs = ", sum(family_dist_coexpr_DT$same_tad == 1  & family_dist_coexpr_DT$topTAD1 == 0 ), ")"), 
              paste0("diff. TADs (# pairs = ", sum(family_dist_coexpr_DT$same_tad == 0), ")"))
  legend("bottomright", legend = legTxt, lty=1, col = c(col1, col2, col3), bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ##################################### PLOT3: COEXPRESSION ~ DISTANCE
  ##################################### SAME FAMILY: SAME TADs VS DIFF TADs
  
  my_xlab <- paste0("Distance between the 2 genes (bp)")
  my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
  
  outFile <- file.path(outFold, paste0(i_fam, "_same_family_sameTAD_vs_diffTADs.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  plot(x = family_dist_coexpr_DT$dist,
       # x = log10(family_dist_coexpr_DT$dist),
       y = family_dist_coexpr_DT$coexpr,
       ylab = my_ylab,
       xlab = my_xlab,
       cex = 0.5,
       pch = 16,
       main = "Coexpr~dist for gene pairs from same family")
  mtext(paste0(curr_dataset, " - ", i_fam, " - same TADs / diff. TADs"), side = 3)
  lines(smooth.spline(x = family_dist_coexpr_DT$dist[family_dist_coexpr_DT$same_tad == 1],
                      y = family_dist_coexpr_DT$coexpr[family_dist_coexpr_DT$same_tad == 1]), 
        col =col1)
  lines(smooth.spline(x = family_dist_coexpr_DT$dist[family_dist_coexpr_DT$same_tad == 0],
                      y = family_dist_coexpr_DT$coexpr[family_dist_coexpr_DT$same_tad == 0]), 
        col =col2)
  legTxt <- c(paste0("same TAD (# pairs = ", sum(family_dist_coexpr_DT$same_tad == 1), ")"), paste0("diff. TADs (# pairs = ", sum(family_dist_coexpr_DT$same_tad == 0), ")"))
  legend("bottomright", legend = legTxt, lty=1, col = c(col1, col2), bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ##################################### PLOT4: COEXPRESSION ~ DISTANCE
  ##################################### SAME FAMILY + SAME TAD: ALL TADs VS TOP TADs
  
 sameTAD_family_dist_coexpr_DT <- family_dist_coexpr_DT[family_dist_coexpr_DT$sameTAD == 1,]
 txt <- paste0("... # of gene pairs in same TAD: ",  nrow(sameTAD_family_dist_coexpr_DT), "\n")
 printAndLog(txt, logFile)
 
 stopifnot(sameTAD_family_dist_coexpr_DT$tad1 == sameTAD_family_dist_coexpr_DT$tad2)
 
 sameTAD_family_dist_coexpr_DT$topTAD <- unlist(sapply(sameTAD_family_dist_coexpr_DT$tad1, function(x) as.numeric(x %in% topTADs )))
 
 outFile <- file.path(outFold, paste0(i_fam, "_same_family_sameTAD_topTADs_vs_otherTADs.", plotType))
 do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x =  sameTAD_family_dist_coexpr_DT$dist,
       # x = log10( sameTAD_family_dist_coexpr_DT$dist),
       y =  sameTAD_family_dist_coexpr_DT$coexpr,
       ylab = my_ylab,
       xlab = my_xlab,
       cex = 0.5,
       pch = 16,
       main = "Coexpr~dist for gene pairs from same family + same TAD")
  mtext(paste0(curr_dataset, " - ", i_fam, " - topTADs / other TADs"), side = 3)

  lines(smooth.spline(x = sameTAD_family_dist_coexpr_DT$dist[sameTAD_family_dist_coexpr_DT$topTAD == 1],
                      y = sameTAD_family_dist_coexpr_DT$coexpr[sameTAD_family_dist_coexpr_DT$topTAD == 1]), 
        col =col1)
  lines(smooth.spline(x = sameTAD_family_dist_coexpr_DT$dist[sameTAD_family_dist_coexpr_DT$topTAD == 0],
                      y = sameTAD_family_dist_coexpr_DT$coexpr[sameTAD_family_dist_coexpr_DT$topTAD == 0]), 
        col =col2)
  legTxt <- c(paste0("topTADs (# pairs = ", sum(sameTAD_family_dist_coexpr_DT$topTAD == 1), ")"), paste0("diff. TADs (# pairs = ", sum(sameTAD_family_dist_coexpr_DT$topTAD == 0), ")"))
  legend("bottomright", legend = legTxt, lty=1, col = c(col1, col2), bty="n")
  
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

