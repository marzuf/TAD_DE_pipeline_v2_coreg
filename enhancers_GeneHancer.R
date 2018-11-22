startTime <- Sys.time()

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

curr_dataset <- "TCGAcrc_msi_mss"
plotType <- "svg"
myHeight <- 7
myWidth <- 7


genehancer_DT <- read.delim("genehancer.csv", header=T, stringsAsFactors = F, sep=",")
nrow(genehancer_DT)
genehancer_DT <- genehancer_DT[grepl("enhancer", tolower(genehancer_DT$feature.name)),]
nrow(genehancer_DT)

pvalCombFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/11_runEmpPvalCombined"), "emp_pval_combined.Rdata")
pval_comb <- eval(parse(text = load(pvalCombFile)))
pval_comb <- p.adjust(pval_comb, method="BH")

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TAD_DT <- read.table(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
TAD_DT <- TAD_DT[grepl("_TAD", TAD_DT$region),]

meanLogfcFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/", curr_dataset, "/3_runMeanTADLogFC"), "all_meanLogFC_TAD.Rdata")
meanLogFC <- eval(parse(text = load(meanLogfcFile)))


source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/main_settings.R"))
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)


load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_GENES_MAPPED/prep_enhancer_DT.Rdata"))
load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_TADS_MAPPED/enhancer_tad_DT.Rdata"))


nEnhancers_vect <- foreach(i = 1:nrow(TAD_DT), .combine='c') %dopar% {
  tmpDT <- genehancer_DT[genehancer_DT$chrom == TAD_DT$chromo[i] & 
                          genehancer_DT$start >= TAD_DT$start[i] &
                           genehancer_DT$end <= TAD_DT$end[i],]
  nrow(tmpDT)
}

TAD_DT$nEnhancers <- nEnhancers_vect

plot(density(TAD_DT$nEnhancers), xlab="# of enhancers/TAD", main="Number of enhancers per TAD")

curr_tads <- intersect(TAD_DT$region, names(pval_comb))

TAD_pval_enhancers_DT <- TAD_DT[TAD_DT$region %in% curr_tads,]
TAD_pval_enhancers_DT$adj_pval_comb <- unlist(sapply(TAD_pval_enhancers_DT$region, function(x) pval_comb[x]))

TAD_pval_enhancers_DT$mean_logFC <- unlist(sapply(TAD_pval_enhancers_DT$region, function(x) meanLogFC[x]))

subTit <- paste0(curr_dataset, " - ", length(curr_tads), " TADs")      

addLeg <- function(xvect, yvect, mypos="topright", corMeth="Pearson") {
  corTest <- cor.test(curr_x, curr_y)
  # legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%.2f", corTest$p.value), ")")
  legTxt <- paste0(corMeth, "'s CC = ", sprintf("%.2f", corTest$estimate), "\n(pval=", sprintf("%1.2e", corTest$p.value), ")")
  legend(mypos, legTxt, bty="n")
}

curr_x <- log10(TAD_pval_enhancers_DT$nEnhancers)
curr_y <- -log10(TAD_pval_enhancers_DT$adj_pval_comb)

outFile <- file.path(outFold, paste0(curr_dataset, "_", "nbrEnhancers", "_", "adjPvalComb.svg"))
svg(outFile, height = myHeight, width = myWidth)
plot( 
      x = curr_x,
      y = curr_y,
      xlab="",
      ylab="",
      main=""
      )
mtext(text = subTit, side=3, font=3)
addLeg(curr_x, curr_y)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


curr_x <- log10(TAD_pval_enhancers_DT$nEnhancers)
curr_y <- TAD_pval_enhancers_DT$mean_logFC

outFile <- file.path(outFold, paste0(curr_dataset, "_", "nbrEnhancers", "_", "meanLogFC.svg"))
svg(outFile, height = myHeight, width = myWidth)
plot( 
      x = curr_x,
      y = curr_y,
      xlab="",
      ylab="",
      main=""
)
mtext(text = subTit, side=3, font=3)
addLeg(curr_x, curr_y)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))






######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))