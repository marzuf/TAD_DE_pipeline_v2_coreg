startTime <- Sys.time()
cat(paste0("> Rscript corr_TAD_adjTAD.R\n"))

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
# Rscript cor_TAD_adjTAD.R
# Rscript cor_TAD_adjTAD.R TCGAcrc_msi_mss 50 hgnc
#  Rscript cor_TAD_adjTAD.R <dataset> <nTADs> <family>
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 1)
  curr_dataset <- args[1]
} else{
  txt <- paste0("> Default parameters:\n")
  curr_dataset <- "TCGAcrc_msi_mss"
}

#outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CORR_TAD_adjTAD",  paste0(curr_dataset, "_nTop", nTopTADs, "_", familyData))
outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CORR_TAD_adjTAD",  paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("corr_TAD_adjTAD.txt"))  
system(paste0("rm -f ", logFile))

printAndLog(txt, logFile)
txt <- paste0("... curr_dataset = ",  curr_dataset, "\n")
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

curr_gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
rm(gene2tadDT)

# curr_TADpos_DT <- TADpos_DT[TADpos_DT$region %in% regionList,]
# all_z_scores <- foreach(curr_tad = regionList) %dopar% {
#   # should I keep only the TADs from regionList ???
#   i_tad <- which(TADpos_DT$region == curr_tad)

all_corr_scores <- foreach(i_tad = seq_len(nrow(TADpos_DT))) %dopar% {  
  
  cat(paste0("... start: ", i_tad, "/", nrow(TADpos_DT), "\n"))
  
  curr_region <- TADpos_DT$region[i_tad]
  curr_chromo <- TADpos_DT$chromo[i_tad]
  
  if(i_tad == 1) {
    region_before <- NULL
  } else {
    if(TADpos_DT$chromo[i_tad-1] == curr_chromo) {
      region_before <- TADpos_DT$region[i_tad-1]
    } else{
      region_before <- NULL
    }
  }
  
  if(i_tad == nrow(TADpos_DT)) {
    region_after <- NULL
  } else {
    if(TADpos_DT$chromo[i_tad+1] == curr_chromo) {
      region_after <- TADpos_DT$region[i_tad+1]
    } else{
      region_after <- NULL
    }
  }
  
  # might be true if I use only regions from regionList and filter adjacent?
  if(is.null(region_after) & is.null(region_before)) stop("error")
  
  all_reg <- c(region_before, curr_region, region_after)  
  
  tad_genes <- curr_gene2tadDT$entrezID[curr_gene2tadDT$region == curr_region]
  adj_tad_genes <- curr_gene2tadDT$entrezID[curr_gene2tadDT$region %in% all_reg]
  
  stopifnot(tad_genes %in% rownames(fpkmDT))
  stopifnot(adj_tad_genes %in% rownames(fpkmDT))
  
  log2fpkm_tad <- log2(fpkmDT[tad_genes,]+0.0001)
  log2fpkm_adjTad <- log2(fpkmDT[adj_tad_genes,]+0.0001)
  
  corr_tad <- cor(t(log2fpkm_tad))
  corr_adjTad <- cor(t(log2fpkm_adjTad))
  
  
  corr_TAD <-  corr_tad[lower.tri(corr_tad)]
  corr_adjTAD <- corr_adjTad[lower.tri(corr_adjTad)]

  list(corr_TAD = corr_TAD,
       corr_adjTAD = corr_adjTAD
       )  
}

names(all_corr_scores) <- TADpos_DT$region

outFile <- file.path(outFold, "all_corr_scores.Rdata")
save(all_corr_scores, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CORR_TAD_adjTAD/TCGAcrc_msi_mss_nTop50_hgnc/all_corr_scores.Rdata")

unlist_corr_TAD <- unlist(lapply(all_corr_scores, function(x) x[["corr_TAD"]]))
# unlist_corr_TAD <- unname(abs(unlist_corr_TAD))
unlist_corr_TAD <- unname(unlist_corr_TAD)

unlist_corr_adjTAD <- unlist(lapply(all_corr_scores, function(x) x[["corr_adjTAD"]]))
# unlist_corr_adjTAD <- unname(abs(unlist_corr_adjTAD))
unlist_corr_adjTAD <- unname(unlist_corr_adjTAD)

outFile <- file.path(outFold, "cmp_density_TAD_adjTAD_log10.svg")
svg(outFile, height = 7, width = 10)
plot_multiDens(list(
  tad = log10(unlist_corr_TAD),
  adj_tad = log10(unlist_corr_adjTAD)
),
my_xlab = "lower tri. corr log2FPKM",
plotTit = "gene pairwise corr. log2FPKM"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "cmp_density_TAD_adjTAD.svg")
svg(outFile, height = 7, width = 10)
plot_multiDens(list(
  tad = unlist_corr_TAD,
  adj_tad = unlist_corr_adjTAD
),
my_xlab = "lower tri. corr log2FPKM",
plotTit = "gene pairwise corr. log2FPKM"
)



dataframeDT <- data.frame(
  values = c(unlist_corr_TAD, unlist_corr_adjTAD),
  type = c(rep("TAD", length(unlist_corr_TAD)), rep("adjTAD", length(unlist_corr_adjTAD)))
)
dataframeDT$type <- factor(dataframeDT$type, levels = c("adjTAD", "TAD"))

my_comparisons <- list( c("TAD", "adjTAD") )

violinP <- ggviolin(dataframeDT, x = "type", y = "values", fill = "type",
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

outFile <- file.path(outFold, paste0("cmp_violin_TAD_adjTAD.svg"))
ggsave(filename = outFile, plot = violinP, height = 7, width=8)
cat(paste0("... written: ", outFile, "\n"))


#foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



