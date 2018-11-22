startTime <- Sys.time()

# Rscript enhancers_TADs_Genes_GeneHancers.R

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

curr_dataset <- "TCGAcrc_msi_mss"
plotType <- "png"
if(plotType == "png") {
  myHeight <- 400
  myWidth <- 400
  myWidthDensity <- 600
} else {
  myHeight <- 7
  myWidth <- 7
  myWidthDensity <- 10
}

source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/main_settings.R"))

nTopTADs <- 50

addLeg <- function(xvect, yvect, mypos="topright", corMeth="Pearson") {
  nona_idx = which(!is.na(xvect) & !is.na(yvect) & !is.infinite(xvect) & !is.infinite(yvect)  )
  corTest <- cor.test(xvect[nona_idx], yvect[nona_idx])
  # legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%.2f", corTest$p.value), ")")
  legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%1.2e", corTest$p.value), ")")
  legend(mypos, legTxt, bty="n")
}


printAndLog <- function(txt, logFile=NULL){
  cat(txt)
  if(!is.null(logFile)) cat(txt, file = logFile, append=T)
}

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  meanDens <- unlist(lapply(size_list, function(x) mean(na.omit(x))))
  
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")\n(mean=", sprintf("%.2f", meanDens), ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset), "ENHANCERS_TADs_GENES_GENEHANCERS")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "enhancers_TADs_genes_logFile.txt")
system(paste0("rm -f ", logFile))

#####################################################
##################################################### load data
#####################################################

pvalCombFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/11_runEmpPvalCombined"), "emp_pval_combined.Rdata")
pval_comb <- eval(parse(text = load(pvalCombFile)))
pval_comb <- p.adjust(pval_comb, method="BH")

topTADs <- names(sort(pval_comb))[1:nTopTADs]

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TAD_DT <- read.table(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
TAD_DT <- TAD_DT[grepl("_TAD", TAD_DT$region),]

meanLogfcFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/3_runMeanTADLogFC"), "all_meanLogFC_TAD.Rdata")
meanLogFC <- eval(parse(text = load(meanLogfcFile)))

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_GENES_MAPPED/prep_enhancer_DT.Rdata"))
load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_TADS_MAPPED/enhancer_tad_DT.Rdata"))

txt <- paste0("... total number of enhancers (all): ", length(unique(prep_enhancer_DT$enhancer)), "\n")
printAndLog(txt, logFile)

txt <- paste0("... total number of enhancers within TADs: ", length(unique(enhancer_tad_DT[grepl("_TAD",enhancer_tad_DT$region),]$enhancer)), "\n")
printAndLog(txt, logFile)

txt <- paste0("... total number of enhancers with mapped genes: ", length(unique(na.omit(prep_enhancer_DT)$enhancer)), "\n")
printAndLog(txt, logFile)

#####################################################
##################################################### prepare data
#####################################################

head(prep_enhancer_DT)
colnames(prep_enhancer_DT)[colnames(prep_enhancer_DT) == "chromo"] <- "enhancer_chromo"
colnames(prep_enhancer_DT)[colnames(prep_enhancer_DT) == "start"] <- "enhancer_start"
colnames(prep_enhancer_DT)[colnames(prep_enhancer_DT) == "end"] <- "enhancer_end"

# they might be duplicated based on gene_remapped !!!
prep_enhancer_DT$connected_gene <- NULL
prep_enhancer_DT <- prep_enhancer_DT[!duplicated(prep_enhancer_DT),]

head(enhancer_tad_DT)
colnames(enhancer_tad_DT)[colnames(enhancer_tad_DT) == "chromo"] <- "enhancer_chromo"
colnames(enhancer_tad_DT)[colnames(enhancer_tad_DT) == "start"] <- "enhancer_start"
colnames(enhancer_tad_DT)[colnames(enhancer_tad_DT) == "end"] <- "enhancer_end"

enhancer_tad_genes_DT <- left_join(prep_enhancer_DT, enhancer_tad_DT, by=c("enhancer_chromo", "enhancer_start", "enhancer_end", "enhancer"))
colnames(enhancer_tad_genes_DT)[colnames(enhancer_tad_genes_DT) == "region"] <- "enhancer_region"

prep_enhancer_DT$connected_gene_renamed <- as.character(prep_enhancer_DT$connected_gene_renamed)
enhancer_tad_genes_tad_DT <- left_join(enhancer_tad_genes_DT, gene2tadDT[,c("entrezID", "region")], by=c("connected_gene_renamed" = "entrezID"))
colnames(enhancer_tad_genes_tad_DT)[colnames(enhancer_tad_genes_tad_DT) == "region"] <- "gene_region"

enhancer_tad_genes_tad_DT$enhancer_midpos <- (enhancer_tad_genes_tad_DT$enhancer_start + enhancer_tad_genes_tad_DT$enhancer_end)*0.5
enhancer_tad_genes_tad_DT$gene_midpos <- (enhancer_tad_genes_tad_DT$connected_gene_start + enhancer_tad_genes_tad_DT$connected_gene_end)*0.5

enhancer_tad_genes_tad_DT$enhancer_gene_dist <- abs(enhancer_tad_genes_tad_DT$enhancer_midpos - enhancer_tad_genes_tad_DT$gene_midpos)

enhancer_tad_genes_tad_DT$enhancer_gene_dist[enhancer_tad_genes_tad_DT$enhancer_chromo != enhancer_tad_genes_tad_DT$connected_gene_chromo] <- NA

enhancer_tad_genes_tad_DT$sameTAD <- as.numeric(enhancer_tad_genes_tad_DT$enhancer_region == enhancer_tad_genes_tad_DT$gene_region)
nrow(enhancer_tad_genes_tad_DT)

#####################################################
##################################################### remove the data for which I could not map any genes
#####################################################

noNA_DT <- na.omit(enhancer_tad_genes_tad_DT)
nrow(noNA_DT)
sum(noNA_DT$sameTAD == 0)/nrow(noNA_DT)
sum(noNA_DT$sameTAD == 1)/nrow(noNA_DT)

txt <- paste0("... number of enhancer-gene pairs (all): ", nrow(noNA_DT), "\n")
printAndLog(txt, logFile)

txt <- paste0("... number of enhancer-gene pairs (diff. TAD): ", sum(noNA_DT$sameTAD == 0), 
              " (", round(sum(noNA_DT$sameTAD == 0)/nrow(noNA_DT) * 100, 2), " %)\n")
printAndLog(txt, logFile)

txt <- paste0("... number of enhancer-gene pairs (same TAD): ", sum(noNA_DT$sameTAD == 1), 
              " (", round(sum(noNA_DT$sameTAD == 1)/nrow(noNA_DT) * 100, 2), " %)\n")
printAndLog(txt, logFile)

#####################################################
##################################################### number of pairs in the same TAD
#####################################################

# pairs in same TAD
sameTAD_noNA_DT <- noNA_DT[noNA_DT$sameTAD == 1,]

tmp <- sameTAD_noNA_DT[,c("enhancer", "connected_gene_renamed")]
stopifnot(!any(duplicated(tmp)))

stopifnot(sameTAD_noNA_DT$enhancer_region == sameTAD_noNA_DT$gene_region)
agg_sameTAD_noNA_DT <- aggregate(sameTAD ~ enhancer_region, data = sameTAD_noNA_DT, FUN=sum)
colnames(agg_sameTAD_noNA_DT)[colnames(agg_sameTAD_noNA_DT) == "sameTAD"] <- "nPairs"
colnames(agg_sameTAD_noNA_DT)[colnames(agg_sameTAD_noNA_DT) == "enhancer_region"] <- "region"

# pairs in diff TAD
diffTAD_noNA_DT <- noNA_DT[noNA_DT$sameTAD == 0,]

tmp <- diffTAD_noNA_DT[,c("enhancer", "connected_gene_renamed")]
stopifnot(!any(duplicated(tmp)))

stopifnot(diffTAD_noNA_DT$enhancer_region != diffTAD_noNA_DT$gene_region)
agg_diffTAD_noNA_DT <- aggregate(sameTAD ~ enhancer_region, data = diffTAD_noNA_DT, FUN=sum)
colnames(agg_diffTAD_noNA_DT)[colnames(agg_diffTAD_noNA_DT) == "sameTAD"] <- "nPairs"
colnames(agg_diffTAD_noNA_DT)[colnames(agg_diffTAD_noNA_DT) == "enhancer_region"] <- "region"

#####################################################
##################################################### multiDensity topTADs
#####################################################
subTit <- paste0(curr_dataset)

outFile <- file.path(outFold, paste0("pvalComb_presence_absence_geneEnhancerPairs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
plot_multiDens(list(
  withEnhancerGenePair = -log10(pval_comb[names(pval_comb) %in% sameTAD_noNA_DT$enhancer_region]),
    noEnhancerGenePair = -log10(pval_comb[!names(pval_comb) %in% sameTAD_noNA_DT$enhancer_region])
),
my_xlab = "-log10 pval combined",
plotTit = paste0("Pval combined and presence/absence enhancer-gene pair(s) in TAD")
)
mtext(text = paste0(subTit), side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("meanLogFC_presence_absence_geneEnhancerPairs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
plot_multiDens(list(
  withEnhancerGenePair = meanLogFC[names(meanLogFC) %in% sameTAD_noNA_DT$enhancer_region],
  noEnhancerGenePair = meanLogFC[!names(meanLogFC) %in% sameTAD_noNA_DT$enhancer_region]
),
my_xlab = "mean log FC",
plotTit = paste0("mean log FC and presence/absence enhancer-gene pair(s) in TAD")
)
mtext(text = paste0(subTit), side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("nbrGeneEnhancerPairs_topTADs_otherTADs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
plot_multiDens(list(
  topTADs = log10(agg_sameTAD_noNA_DT$nPairs[agg_sameTAD_noNA_DT$region %in% topTADs]),
  otherTADs = log10(agg_sameTAD_noNA_DT$nPairs[!agg_sameTAD_noNA_DT$region %in% topTADs])
),
my_xlab = "# gene-enhancer pairs within TAD",
plotTit = paste0("Nbr enhancer-gene pair(s) by TAD")
)
mtext(text = paste0(subTit), side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#####################################################
##################################################### count number of elements by TAD - all data
#####################################################

nEnhancer_TAD_DT_allData <- aggregate(enhancer~enhancer_region, data=noNA_DT, FUN=function(x) length(unique(x)))

nMappedGenes_TAD_DT_allData <- aggregate(connected_gene_renamed~gene_region, data=noNA_DT, FUN=function(x) length(unique(x)))

nEnhancer_TAD_DT_sameTAD <- aggregate(enhancer~enhancer_region, 
                                      data=noNA_DT[noNA_DT$enhancer_region == noNA_DT$gene_region,],
                                      FUN=function(x) length(unique(x)))

nMappedGenes_TAD_DT_sameTAD <- aggregate(connected_gene_renamed~gene_region, 
                                         data=noNA_DT[noNA_DT$enhancer_region == noNA_DT$gene_region,],
                                         FUN=function(x) length(unique(x)))


nEnhancer_TAD_DT_diffTAD <- aggregate(enhancer~enhancer_region, 
                                      data=noNA_DT[noNA_DT$enhancer_region != noNA_DT$gene_region,],
                                      FUN=function(x) length(unique(x)))

nMappedGenes_TAD_DT_diffTAD <- aggregate(connected_gene_renamed~gene_region, 
                                         data=noNA_DT[noNA_DT$enhancer_region != noNA_DT$gene_region,],
                                         FUN=function(x) length(unique(x)))


enhancer_tad_genes_tad_DT_allData <- enhancer_tad_genes_tad_DT
enhancer_tad_genes_tad_DT_sameTAD <- enhancer_tad_genes_tad_DT[enhancer_tad_genes_tad_DT$gene_region == enhancer_tad_genes_tad_DT$enhancer_region,]
enhancer_tad_genes_tad_DT_diffTAD <- enhancer_tad_genes_tad_DT[enhancer_tad_genes_tad_DT$gene_region != enhancer_tad_genes_tad_DT$enhancer_region,]

all_datasets <- c("allData", "sameTAD", "diffTAD")

for(data_type in all_datasets) {
  
  nEnhancer_TAD_DT <- eval(parse(text = paste0("nEnhancer_TAD_DT", "_", data_type)))
  nMappedGenes_TAD_DT <- eval(parse(text = paste0("nMappedGenes_TAD_DT", "_", data_type)))
    
  colnames(nEnhancer_TAD_DT)[colnames(nEnhancer_TAD_DT) == "enhancer"] <- "nbrUniqueEnhancers"
  colnames(nEnhancer_TAD_DT)[colnames(nEnhancer_TAD_DT) == "enhancer_region"] <- "region"
  nEnhancer_TAD_DT <- nEnhancer_TAD_DT[grepl("_TAD", nEnhancer_TAD_DT$region),]
  
  
  colnames(nMappedGenes_TAD_DT)[colnames(nMappedGenes_TAD_DT) == "connected_gene_renamed"] <- "nbrUniqueGenes"
  colnames(nMappedGenes_TAD_DT)[colnames(nMappedGenes_TAD_DT) == "gene_region"] <- "region"
  nMappedGenes_TAD_DT <- nMappedGenes_TAD_DT[grepl("_TAD", nMappedGenes_TAD_DT$region),]
  
  nElements_TAD_DT <- inner_join(nMappedGenes_TAD_DT, nEnhancer_TAD_DT, by="region")
  nElements_TAD_DT$pvalComb <- unlist(sapply(nElements_TAD_DT$region, function(x) pval_comb[x]))
  nElements_TAD_DT$meanLogFC <- unlist(sapply(nElements_TAD_DT$region, function(x) meanLogFC[x]))
  nElements_TAD_DT$nbrUniqueEnhancers_log10 <- log10(nElements_TAD_DT$nbrUniqueEnhancers)
  nElements_TAD_DT$nbrUniqueGenes_log10 <- log10(nElements_TAD_DT$nbrUniqueGenes)
  nElements_TAD_DT$ratio_enhancer_over_genes <- nElements_TAD_DT$nbrUniqueEnhancers/nElements_TAD_DT$nbrUniqueGenes
  nElements_TAD_DT$ratio_log10enhancer_over_log10genes <- nElements_TAD_DT$nbrUniqueEnhancers_log10/nElements_TAD_DT$nbrUniqueGenes_log10
  
  curr_enhancer_tad_genes_tad_DT <- eval(parse(text = paste0("enhancer_tad_genes_tad_DT", "_", data_type)))
    
  #####################################################
  ##################################################### start plotting
  #####################################################
  
  subTit <- paste0(curr_dataset, " - ", data_type)
  
  ### GENERAL PLOTS _ all pairs
  
  outFile <- file.path(outFold, paste0("dist_genes_enhancers_density_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
  plot(density(na.omit(curr_enhancer_tad_genes_tad_DT$enhancer_gene_dist)), 
       xlab="Distance enhancers/genes", main="Distance betw. genes and enhancers")
  legend("topright", legend = paste0("n = ", length(na.omit(curr_enhancer_tad_genes_tad_DT$enhancer_gene_dist))))
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFold, paste0("dist_genes_enhancers_log10_density_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
  plot(density(log10(na.omit(curr_enhancer_tad_genes_tad_DT$enhancer_gene_dist))), 
       xlab="Distance (bp)  [log10]", main="Distance betw. genes and enhancers")
  mtext(text = paste0(subTit), side = 3)
  legend("topright", legend = paste0("n = ", length(log10(na.omit(curr_enhancer_tad_genes_tad_DT$enhancer_gene_dist)))) )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFold, paste0("nUniqueEnhancers_nUniqueGenes_perTAD_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$nbrUniqueEnhancers,
       y = nElements_TAD_DT$nbrUniqueGenes,
       pch=16, cex=0.7,
       xlab="nbr unique enhancers",
       ylab="nbr unique genes",
       main ="Nbr elements per TAD")
  addLeg(xvect = nElements_TAD_DT$nbrUniqueEnhancers, yvect = nElements_TAD_DT$nbrUniqueGenes,
         mypos="topright", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("nUniqueEnhancers_nUniqueGenes_perTAD_log10_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$nbrUniqueEnhancers_log10,
       y = nElements_TAD_DT$nbrUniqueGenes_log10,
       pch=16, cex=0.7,
       xlab="nbr unique enhancers (log10)",
       ylab="nbr unique genes (log10)",
       main ="Nbr elements per TAD (log10)")
  addLeg(xvect = nElements_TAD_DT$nbrUniqueEnhancers_log10, yvect = nElements_TAD_DT$nbrUniqueGenes_log10,
         mypos="topleft", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nEnhancers_over_nGenes_perTAD_density_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
  plot(density(na.omit(nElements_TAD_DT$ratio_enhancer_over_genes)), 
       xlab="Enhancers/genes ratio", main="Ratio enhancers/genes by TAD")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nEnhancersLog10_over_nGenesLog10_perTAD_density_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidthDensity))
  plot(density(na.omit(nElements_TAD_DT$ratio_log10enhancer_over_log10genes)), 
       xlab="Enhancers/genes ratio (log10)", main="Ratio enhancers (log10)/genes (log10) by TAD")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nUniqueEnhancersLog10_vs_pvalComb_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$nbrUniqueEnhancers_log10,
       y = -log10(nElements_TAD_DT$pvalComb),
       pch=16, cex=0.7,
       xlab="nbr unique enhancers (log10)",
       ylab="-log10 pval comb.",
       main ="Nbr enhancers and combined pval")
  addLeg(xvect = nElements_TAD_DT$nbrUniqueEnhancers_log10, yvect = -log10(nElements_TAD_DT$pvalComb),
         mypos="topleft", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nUniqueGenesLog10_vs_pvalComb_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$nbrUniqueGenes_log10,
       y = -log10(nElements_TAD_DT$pvalComb),
       pch=16, cex=0.7,
       xlab="nbr unique genes (log10)",
       ylab="-log10 pval comb.",
       main ="Nbr genes and combined pval")
  addLeg(xvect = nElements_TAD_DT$nbrUniqueGenes_log10, yvect = -log10(nElements_TAD_DT$pvalComb),
         mypos="topleft", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))

  outFile <- file.path(outFold, paste0("ratio_nEnhancers_over_nGenes_vs_pvalComb_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$ratio_enhancer_over_genes,
       y = -log10(nElements_TAD_DT$pvalComb),
       pch=16, cex=0.7,
       xlab="ratio enhancers/genes",
       ylab="-log10 pval comb.",
       main ="Ratio enhancers/genes and combined pval")
  addLeg(xvect = nElements_TAD_DT$ratio_enhancer_over_genes, yvect = -log10(nElements_TAD_DT$pvalComb),
         mypos="topright", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))

  
  outFile <- file.path(outFold, paste0("ratio_nEnhancersLog10_over_nGenesLog10_vs_pvalComb_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$ratio_log10enhancer_over_log10genes,
       y = -log10(nElements_TAD_DT$pvalComb),
       pch=16, cex=0.7,
       xlab="ratio enhancers(log10)/genes(log10)",
       ylab="-log10 pval comb.",
       main ="Ratio enhancers/genes and combined pval")
  addLeg(xvect = nElements_TAD_DT$ratio_log10enhancer_over_log10genes, yvect = -log10(nElements_TAD_DT$pvalComb),
         mypos="topright", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nUniqueEnhancersLog10_vs_meanLogFC_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$nbrUniqueEnhancers_log10,
       y = nElements_TAD_DT$meanLogFC,
       pch=16, cex=0.7,
       xlab="nbr unique enhancers (log10)",
       ylab= "mean log FC",
       main ="Nbr enhancers and mean log FC")
  addLeg(xvect = nElements_TAD_DT$nbrUniqueEnhancers_log10, yvect = nElements_TAD_DT$meanLogFC,
        mypos="topleft", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nUniqueGenesLog10_vs_meanLogFC_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$nbrUniqueGenes_log10,
       y = nElements_TAD_DT$meanLogFC,
       pch=16, cex=0.7,
       xlab="nbr unique genes (log10)",
       ylab="mean log FC",
       main ="Nbr genes and mean log FC")
  addLeg(xvect = nElements_TAD_DT$nbrUniqueGenes_log10, yvect = nElements_TAD_DT$meanLogFC,
         mypos="topleft", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0("ratio_nEnhancersLog10_over_nGenesLog10_vs_meanLogFC_scatter_", data_type, ".", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = nElements_TAD_DT$ratio_log10enhancer_over_log10genes,
       y = nElements_TAD_DT$meanLogFC,
       pch=16, cex=0.7,
       xlab="ratio enhancers(log10)/genes(log10)",
       ylab= "mean log FC",
       main ="Ratio enhancers/genes and mean log FC")
  addLeg(xvect = nElements_TAD_DT$ratio_log10enhancer_over_log10genes, yvect = nElements_TAD_DT$meanLogFC,
         mypos="topright", corMeth="Pearson")
  mtext(text = paste0(subTit), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))

}
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
