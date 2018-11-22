startTime <- Sys.time()
cat(paste0("> Rscript coexpr_dist_v2.R\n"))

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
buildTable <- TRUE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"


### RETRIEVE FROM COMMAND LINE
# Rscript coexpr_dist_v2.R
# Rscript coexpr_dist_v2.R TCGAcrc_msi_mss hgnc
#  Rscript coexpr_dist_v2.R <dataset> <family>
# top-ranking:
# Rscript coexpr_dist_v2.R TCGAcrc_msi_mss hgnc
# Rscript coexpr_dist_v2.R GSE74927_neg_pos hgnc
# Rscript coexpr_dist_v2.R GSE102073_stic_nostic hgnc   
# worst-ranking
# Rscript coexpr_dist_v2.R GSE65540_before_after hgnc   # ok
# Rscript coexpr_dist_v2.R GSE84231_lhb_rhb hgnc    # running
# Rscript coexpr_dist_v2.R GSE86356_tibMD1_tibNorm hgnc    # waiting !s
#  Rscript coexpr_dist_v2.R <dataset> <family>


args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 2)
  curr_dataset <- args[1]

  familyData <- args[2]
} else{
  txt <- paste0("> Default parameters:\n")

  curr_dataset <- "TCGAcrc_msi_mss"
  familyData <- "hgnc"
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "COEXPR_DIST_v2",  paste0(curr_dataset,  "_", familyData))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("coexpr_dist_v2_logFile_plot.txt"))  
system(paste0("rm -f ", logFile))

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


col1 <- "darkslateblue"
col2 <-  "darkorange1"

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
    # outFold <- "/media/electron//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST/TCGAcrc_msi_mss_nTop50_hgnc"
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
  
  # consider only pairs up to distLimit
  family_dist_coexpr_DT <- family_dist_coexpr_DT[family_dist_coexpr_DT$dist <= distLimit,]
  txt <- paste0("... # of gene pairs within ", distLimit, " bp: ", nrow(family_dist_coexpr_DT), "\n")
  printAndLog(txt, logFile)
  
  family_dist_coexpr_DT$tad1 <- unlist(sapply(family_dist_coexpr_DT$gene1, function(x)
    gene2tadDT$region[gene2tadDT$entrezID == x]))
  
  family_dist_coexpr_DT$tad2 <- unlist(sapply(family_dist_coexpr_DT$gene2, function(x)
    gene2tadDT$region[gene2tadDT$entrezID == x]))
  
  family_dist_coexpr_DT$same_tad <- as.numeric(family_dist_coexpr_DT$tad1 == family_dist_coexpr_DT$tad2)
  stopifnot(family_dist_coexpr_DT$same_tad == family_dist_coexpr_DT$sameTAD)
  
  

  
  ##################################### PLOT2: COEXPRESSION ~ DISTANCE
  ##################################### SAME FAMILY: SAME TADs VS DIFF TADs
  my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
  
  my_xlab <- paste0("Distance between the 2 genes (kb)")
  
                    
  family_dist_coexpr_DT$sameTAD <- as.character(family_dist_coexpr_DT$sameTAD)
  family_dist_coexpr_DT$sameTAD_lab <- ifelse(family_dist_coexpr_DT$sameTAD == "0", "diff. TAD", "same TAD")
  family_dist_coexpr_DT$dist_kb <- family_dist_coexpr_DT$dist/1000
  
 
 coexpr_dist_plot <- my_ggscatterhist(
   family_dist_coexpr_DT, 
   x = "dist_kb", 
   xlab=my_xlab,
   ylab=my_ylab,
   title = paste0("Gene pair expr. corr. vs. dist. - ", curr_dataset),
   subtitle=paste0("(", i_fam,")"),
   legend.title="",
   y = "coexpr",
   point = FALSE,
   # rug=TRUE,
   add = fitMeth,
   color = "sameTAD_lab", size = 3, alpha = 0.6,
   palette = c(col1, col2),
   x_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
   x_margin.ggtheme = theme_minimal(),
   x_margin.plot = "density",
   y_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
   y_margin.ggtheme = theme_minimal(),
   y_margin.plot = "boxplot",
   plot_xmargin = TRUE,
   plot_ymargin=TRUE,
   ymargin_as_xmargin = FALSE,
   # global_margins_cm = c(0, 0, 0.25, 0.25)
   global_margins_cm = c(0.25, 0.25, 0.25, 0.25)
 )
 
 # to add title at bottom:
 # coexpr_dist_plot+labs(caption="Bottom Title") + 
 #   theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))
 
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

