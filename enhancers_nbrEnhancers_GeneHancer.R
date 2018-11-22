startTime <- Sys.time()

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript enhancers_nbrEnhancers_GeneHancers.R

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/main_settings.R"))

curr_dataset <- "TCGAcrc_msi_mss"
plotType <- "svg"
myHeight <- 7
myWidth <- 7
nTopTADs <- 50

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


outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset), "NBR_ENHANCERS_GENEHANCERS")
system(paste0("mkdir -p ", outFold))

# genehancer_DT <- read.delim("genehancer.csv", header=T, stringsAsFactors = F, sep=",")
# nrow(genehancer_DT)
# genehancer_DT <- genehancer_DT[grepl("enhancer", tolower(genehancer_DT$feature.name)),]
# nrow(genehancer_DT)

pvalCombFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/11_runEmpPvalCombined"), "emp_pval_combined.Rdata")
pval_comb <- eval(parse(text = load(pvalCombFile)))
pval_comb <- p.adjust(pval_comb, method="BH")
topTADs <- names(sort(pval_comb))[1:nTopTADs]


TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TAD_DT <- read.table(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
TAD_DT <- TAD_DT[grepl("_TAD", TAD_DT$region),]

meanLogfcFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/3_runMeanTADLogFC"), "all_meanLogFC_TAD.Rdata")
meanLogFC <- eval(parse(text = load(meanLogfcFile)))

load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_TADS_MAPPED/enhancer_tad_DT.Rdata"))
stopifnot(!any(duplicated(enhancer_tad_DT$enhancer)))
agg_enhancer_TAD <- aggregate(enhancer ~ region, FUN=length, data=enhancer_tad_DT)
colnames(agg_enhancer_TAD)[colnames(agg_enhancer_TAD) == "enhancer"] <- "nEnhancers"
agg_enhancer_TAD$region <- as.character(agg_enhancer_TAD$region)

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

curr_tads <- intersect(TAD_DT$region, names(pval_comb))

TAD_pval_enhancers_DT <- TAD_DT[TAD_DT$region %in% curr_tads,]
TAD_pval_enhancers_DT$adj_pval_comb <- unlist(sapply(TAD_pval_enhancers_DT$region, function(x) pval_comb[x]))

TAD_pval_enhancers_DT$mean_logFC <- unlist(sapply(TAD_pval_enhancers_DT$region, function(x) meanLogFC[x]))

TAD_pval_enhancers_DT$nEnhancers <- unlist(sapply(TAD_pval_enhancers_DT$region, function(x) ifelse(x %in% agg_enhancer_TAD$region, 
                                                                                                   agg_enhancer_TAD$nEnhancers[agg_enhancer_TAD$region == x], 0.1)))

subTit <- paste0(curr_dataset, " - ", length(curr_tads), " TADs")      

addLeg <- function(xvect, yvect, mypos="topright", corMeth="Pearson") {
  corTest <- cor.test(xvect, yvect)
  # legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%.2f", corTest$p.value), ")")
  legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%1.2e", corTest$p.value), ")")
  legend(mypos, legTxt, bty="n")
}


outFile <- file.path(outFold, paste0(curr_dataset, "_", "nbrEnhancers", "_", "density.svg"))
svg(outFile, height = myHeight, width = myWidth)
plot( 
  density(log10(TAD_pval_enhancers_DT$nEnhancers)),
  xlab="# enhancers/TAD (log10)",
  main=""
)
mtext(text = subTit, side=3, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(curr_dataset, "_", "nbrEnhancers", "_topTADs_", "density.svg"))
svg(outFile, height = myHeight, width = myWidth)
plot_multiDens( list(
  topTADs = log10(TAD_pval_enhancers_DT$nEnhancers[TAD_pval_enhancers_DT$region %in% topTADs]) ,
  otherTADs =log10(TAD_pval_enhancers_DT$nEnhancers[!TAD_pval_enhancers_DT$region %in% topTADs])),
  my_xlab = "# enhancers/TAD (log10)",
  plotTit = paste0(""),
  legPos="topleft"
)
mtext(text = subTit, side=3, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




curr_x <- log10(TAD_pval_enhancers_DT$nEnhancers)
curr_y <- -log10(TAD_pval_enhancers_DT$adj_pval_comb)

outFile <- file.path(outFold, paste0(curr_dataset, "_", "nbrEnhancers", "_", "adjPvalComb.svg"))
svg(outFile, height = myHeight, width = myWidth)
plot( 
  x = curr_x,
  y = curr_y,
  xlab="nbr enhancers (log10)",
  ylab="TAD adj. pval combined (-log10)",
  main="# enhancers and TAD pval combined",
  pch=16, cex=0.7
)
mtext(text = subTit, side=3, font=3)
addLeg(curr_x, curr_y, mypos="topleft")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

curr_x <- log10(TAD_pval_enhancers_DT$nEnhancers)
curr_y <- TAD_pval_enhancers_DT$mean_logFC

outFile <- file.path(outFold, paste0(curr_dataset, "_", "nbrEnhancers", "_", "meanLogFC.svg"))
svg(outFile, height = myHeight, width = myWidth)
plot( 
  x = curr_x,
  y = curr_y,
  xlab="nbr enhancers (log10)",
  ylab="TAD mean log FC",
  main="# enhancers and TAD mean logFC",
  pch=16, cex=0.7
)
mtext(text = subTit, side=3, font=3)
addLeg(curr_x, curr_y, mypos="topleft")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))






######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))