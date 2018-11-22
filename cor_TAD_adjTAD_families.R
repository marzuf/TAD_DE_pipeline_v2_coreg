startTime <- Sys.time()
cat(paste0("> Rscript corr_TAD_adjTAD_families.R\n"))

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
buildTable <- TRUE

### RETRIEVE FROM COMMAND LINE
# Rscript corTAD_adjTAD_families.R
# Rscript corTAD_adjTAD_families.R TCGAcrc_msi_mss 50 hgnc
# Rscript corTAD_adjTAD_families.R <dataset> <nTADs> <family>
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

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CORR_TAD_adjTAD_FAM",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("corr_TAD_adjTAD_families.txt"))  
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


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TADpos_DT <- read.delim(TADpos_file, header=F, stringsAsFactors = F, col.names=c("chromo", "region", "start", "end"))
TADpos_DT <- TADpos_DT[grepl("_TAD", TADpos_DT$region),]
TADpos_DT <- TADpos_DT[order(TADpos_DT$chromo, TADpos_DT$start),]

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)

#============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
script0_name <- "0_prepGeneData"
geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
txt <- paste0("... number of genes used in the ", curr_dataset, " pipeline: ", length(geneList), "\n")
printAndLog(txt, logFile)
stopifnot(geneList %in% gene2tadDT$entrezID)

fpkmDT <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "rna_fpkmDT.Rdata"))))
stopifnot(names(geneList) %in% rownames(fpkmDT))
fpkmDT <- fpkmDT[rownames(fpkmDT) %in% names(geneList),]
stopifnot(nrow(fpkmDT) == length(geneList))
rownames(fpkmDT) <- geneList[rownames(fpkmDT)]
stopifnot(setequal(geneList, rownames(fpkmDT)))

regionList_file <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/TCGAcrc_msi_mss/0_prepGeneData/pipeline_regionList.Rdata")
stopifnot(file.exists(regionList_file))
regionList <- eval(parse(text = load(regionList_file)))
#==============================


hasAdjacent <- function(mytad, alltads, tadposDT) {
  idx <- which(tadposDT$region == mytad)
  allidx <- as.numeric(sapply(alltads, function(x) which(tadposDT$region == x)))
  chrom <- tadposDT$chromo[idx]
  allchrom <- tadposDT$chromo[allidx]
  return(any(abs(idx-allidx) == 1 & chrom == allchrom))
}

getAdjacent <- function(mytad, alltads, tadposDT) {
  idx <- which(tadposDT$region == mytad)
  allidx <- as.numeric(sapply(alltads, function(x) which(tadposDT$region == x)))
  chrom <- tadposDT$chromo[idx]
  allchrom <- tadposDT$chromo[allidx]
  adjIdx <- which(abs(idx-allidx) == 1 & chrom == allchrom)
  return(alltads[adjIdx])
}

# getAdjacent(mytad="chr10_TAD1", c("chr5_TAD3", "chr4_TAD2", "chr9_TAD1"), tadposDT=TADpos_DT)
# getAdjacent(mytad="chr10_TAD1", c("chr10_TAD2", "chr4_TAD2", "chr9_TAD1"), tadposDT=TADpos_DT)
# getAdjacent(mytad="chr1_TAD364", c("chr10_TAD1", "chr4_TAD2", "chr9_TAD1"), tadposDT=TADpos_DT)

# hasAdjacent(mytad="chr10_TAD1", c("chr5_TAD3", "chr4_TAD2", "chr9_TAD1"), tadposDT=TADpos_DT)
# hasAdjacent(mytad="chr10_TAD1", c("chr10_TAD2", "chr4_TAD2", "chr9_TAD1"), tadposDT=TADpos_DT)
# hasAdjacent(mytad="chr1_TAD364", c("chr10_TAD1", "chr4_TAD2", "chr9_TAD1"), tadposDT=TADpos_DT)

curr_gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
rm(gene2tadDT)

curr_familyDT <- familyDT[familyDT$entrezID %in% geneList,]
rm(familyDT)

all_fams <- paste0(familyData, c("_family", "_family_short"))

# iterate over _family and _family_short
for(i_fam in all_fams) {
  
  cat(paste0("> START ", i_fam, "\n"))
  
  all_fams <- unique(curr_familyDT[,i_fam])
  
  # iterate over all families
  famCorrList <- foreach(curr_fam = all_fams) %dopar% {
    
    cat(paste0("... Start family: ", curr_fam, "\n"))
    
    curr_fam_tads <- unique(curr_familyDT$region[curr_familyDT[,i_fam] == curr_fam])
    
    tads_with_adj <- curr_fam_tads[sapply(curr_fam_tads, function(x) hasAdjacent(x, curr_fam_tads, TADpos_DT))]
    
    if(length(tads_with_adj) == 0) {
      cat("...... no adjacent TAD --> SKIP \n")
      return(NULL)
    } 
    # iterate over TADs that have adjacent TADs
    tad_corr_list <- foreach(tad = tads_with_adj) %do% {
      
      fam_tad_genes <- intersect(curr_familyDT$entrezID[curr_familyDT[,i_fam] == curr_fam], curr_gene2tadDT$entrezID[curr_gene2tadDT$region == tad])
      
      if(length(fam_tad_genes) == 1) {
        cat("...... no adjacent TAD --> ONLY ONE GENE IN TAD \n")
        return(NULL)
      }
      
      adjTADs <- getAdjacent(tad, alltads=curr_fam_tads, tadposDT=TADpos_DT)
      adjTADs <- c(tad, adjTADs)
      
      fam_adjTAD_genes <- intersect(curr_familyDT$entrezID[curr_familyDT[,i_fam] == curr_fam], curr_gene2tadDT$entrezID[curr_gene2tadDT$region %in% adjTADs])
      
      stopifnot(fam_tad_genes %in% rownames(fpkmDT))
      stopifnot(fam_adjTAD_genes %in% rownames(fpkmDT))
      
      log2fpkm_tad <- log2(fpkmDT[fam_tad_genes,]+0.0001)
      log2fpkm_adjTad <- log2(fpkmDT[fam_adjTAD_genes,]+0.0001)
      
      corr_tad <- cor(t(log2fpkm_tad))
      corr_adjTad <- cor(t(log2fpkm_adjTad))
      
      corr_TAD_tri <-  corr_tad[lower.tri(corr_tad)]
      corr_adjTAD_tri <- corr_adjTad[lower.tri(corr_adjTad)]
      
      list(
        corr_TAD_tri = corr_TAD_tri,
        corr_adjTAD_tri = corr_adjTAD_tri
      )
      
    } # end iterating over TADs for the current families
    
    names(tad_corr_list) <- tads_with_adj
    
    tad_corr_list
    
  } # end iterating over families
  
  names(famCorrList) <- all_fams
  
  outFile <- file.path(outFold, paste0(i_fam, "_famCorrList.Rdata"))
  save(famCorrList, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  # load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CORR_TAD_adjTAD_FAM/TCGAcrc_msi_mss_nTop50_hgnc/hgnc_family_famCorrList.Rdata")
  
  unlist_fam <- unlist(famCorrList, recursive = F)
  unlist_fam <- Filter(f =function(x)!is.null(x), unlist_fam)
  unlist_corr_TAD <- unlist(lapply(unlist_fam, function(x) x[["corr_TAD_tri"]]))
  unlist_corr_adjTAD <- unlist(lapply(unlist_fam, function(x) x[["corr_adjTAD_tri"]]))
  
  outFile <- file.path(outFold, paste0(i_fam, "cmp_density_TAD_adjTAD_log10.svg"))
  svg(outFile, height = 7, width = 10)
  plot_multiDens(list(
    tad = log10(unlist_corr_TAD),
    adj_tad = log10(unlist_corr_adjTAD)
  ),
  my_xlab = "lower tri. corr log2FPKM [log10]",
  plotTit = "gene pairwise corr. log2FPKM"
  )
  mtext(text = paste0(i_fam, " - ", curr_dataset), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(i_fam, "cmp_density_TAD_adjTAD.svg"))
  svg(outFile, height = 7, width = 10)
  plot_multiDens(list(
    tad = unlist_corr_TAD,
    adj_tad = unlist_corr_adjTAD
  ),
  my_xlab = "lower tri. corr log2FPKM",
  plotTit = "gene pairwise corr. log2FPKM"
  )
  mtext(text = paste0(i_fam, " - ", curr_dataset), side = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  dataframeDT <- data.frame(
    values = c(unlist_corr_TAD, unlist_corr_adjTAD),
    type = c(rep("TAD", length(unlist_corr_TAD)), rep("adjTAD", length(unlist_corr_adjTAD)))
  )
  dataframeDT$type <- factor(dataframeDT$type, levels = c("adjTAD", "TAD"))
  
  my_comparisons <- list( c("TAD", "adjTAD") )
  
  violinP <- ggviolin(dataframeDT, x = "type", y = "values", fill = "type",
           title = paste0(i_fam, " - ", curr_dataset),            
           legend.title = "",
           # yscale = "log10",
           # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
           palette = c("#00AFBB", "#FC4E07"),
           xlab ="",
           ylab ="gene log2FPKM correlation",
           add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
    # stat_compare_means(label.y.npc="top", label.x.npc = "left")    # Add global the p-value 
    stat_compare_means(label.y = max(dataframeDT$values), label.x = 0.5)
  
  violinP <- violinP +
    geom_hline(yintercept = mean(dataframeDT$values[dataframeDT$type=="adjTAD"]), linetype=2, color = "#00AFBB", size=1)+
    geom_hline(yintercept = mean(dataframeDT$values[dataframeDT$type=="TAD"]), linetype=2, color = "#FC4E07", size=1)
  
  if(SSHFS) violinP
    
  outFile <- file.path(outFold, paste0(i_fam, "cmp_violin_TAD_adjTAD.svg"))
  ggsave(filename = outFile, plot = violinP, height = 7, width=8)
  cat(paste0("... written: ", outFile, "\n"))
}




######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



