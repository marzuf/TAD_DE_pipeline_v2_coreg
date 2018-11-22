startTime <- Sys.time()

## TOP 3
# Rscript gene_enhancer_TAD_rank.R TCGAcrc_msi_mss
# Rscript gene_enhancer_TAD_rank.R GSE74927_neg_pos
# Rscript gene_enhancer_TAD_rank.R GSE102073_stic_nostic
# # LAST 3
# Rscript gene_enhancer_TAD_rank.R GSE65540_before_after
# Rscript gene_enhancer_TAD_rank.R GSE84231_lhb_rhb
# Rscript gene_enhancer_TAD_rank.R GSE86356_tibMD1_tibNorm

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- "TCGAcrc_msi_mss"
stopifnot(length(args) == 1)
curr_dataset <- args[1]

caller <- "TopDom"
nTop <- 50
topRatio <- 5
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)
nTopTADs <- 50


cat(paste0("... for dataset =\t", curr_dataset, "\n"))
cat("!!! HARD-CODED:")
cat(paste0("... plot top pval TADs nTop =\t", nTop, "\n"))    
cat(paste0("... print top ratio TADs topRatio =\t", topRatio, "\n"))    

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "GENE_ENHANCER_TAD_RANK",  paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))

topTableFile <- file.path(outFold, paste0("top", topRatio, "TADs_dt.txt"))
system(paste0("rm -f ", topTableFile))

# !! source(settingFile) here INSTEAD
# settingFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput/run_settings_", curr_dataset, ".R"))
# stopifnot(file.exists(settingFile))
mainSettingFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/main_settings.R"))
stopifnot(file.exists(mainSettingFile))
source(mainSettingFile)
# gene2tadDTFile <- "../../gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt"
stopifnot(file.exists(gene2tadDT_file))
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)

enhancerTAD_dt <- eval(parse(text = load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_TADS_MAPPED/enhancer_tad_DT.Rdata"))))
stopifnot(!any(duplicated(enhancerTAD_dt$enhancer)))
enhancerTAD_dt$region <- as.character(enhancerTAD_dt$region)

pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"))

geneListFile <- file.path(pipDir, curr_dataset, "0_prepGeneData", "pipeline_geneList.Rdata")
stopifnot(file.exists(geneListFile))
geneList <- eval(parse(text = load(geneListFile)))


regionListFile <- file.path(pipDir, curr_dataset, "0_prepGeneData", "pipeline_regionList.Rdata")
stopifnot(file.exists(regionListFile))
regionList <- eval(parse(text = load(regionListFile)))

pvalFile <- file.path(pipDir, curr_dataset, "11_runEmpPvalCombined", "emp_pval_combined.Rdata")
stopifnot(file.exists(pvalFile))
pvalComb <- eval(parse(text = load(pvalFile)))
adj_pvalComb <- p.adjust(pvalComb, method="BH")
adj_pvalComb_sorted <- sort(adj_pvalComb)
adj_pvalComb_rank <- rank(adj_pvalComb, ties="min")

curr_enhancerTAD_DT <- enhancerTAD_dt[enhancerTAD_dt$region %in% regionList,]
nbrEnhancersByTAD_DT <- data.frame(
  TAD = names(table(curr_enhancerTAD_DT$region)),
  nbrEnhancers = as.numeric(table(curr_enhancerTAD_DT$region)),
  stringsAsFactors = FALSE
)


curr_g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
nbrGenesByTAD_DT <- data.frame(
  TAD = names(table(curr_g2t_DT$region)),
  nbrGenes = as.numeric(table(curr_g2t_DT$region)),
  stringsAsFactors = FALSE
)

genesEnhancersTAD_DT <- merge(nbrEnhancersByTAD_DT, nbrGenesByTAD_DT, by="TAD")
genesEnhancersTAD_DT$geneEnhancerRatio <- genesEnhancersTAD_DT$nbrGenes/genesEnhancersTAD_DT$nbrEnhancers
genesEnhancersTAD_DT <- genesEnhancersTAD_DT[as.character(genesEnhancersTAD_DT$TAD) %in% names(adj_pvalComb),]

# order according to pvalComb
genesEnhancersTAD_DT <- genesEnhancersTAD_DT[match(names(adj_pvalComb_sorted)[names(adj_pvalComb_sorted) %in% genesEnhancersTAD_DT$TAD], genesEnhancersTAD_DT$TAD),]

outFile <- file.path(outFold, paste0("genes_enhancers_ratio_ranked_TADs_all.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(genesEnhancersTAD_DT$geneEnhancerRatio, 
     ylab = "TAD genes/enhancers ratio",
     xlab = "TAD adj. combined pval rank",
     main = paste0(curr_dataset, "TADS: genes/enhancers vs. pval rank"),
     type="h")
mtext(text = paste0("all TADs (n = ", nrow(genesEnhancersTAD_DT), ")"), side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


stopifnot(nTop>=1)
outFile <- file.path(outFold, paste0("genes_enhancers_ratio_ranked_TADs_top", nTop, ".", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(genesEnhancersTAD_DT$geneEnhancerRatio[1:nTop], 
     ylab = "TAD genes/enhancers ratio",
     xlab = "TAD adj. combined pval rank",
     main = paste0(curr_dataset, "TADS: genes/enhancers vs. pval rank"),
     type="h")
mtext(text = paste0("top TADs (n = ", nTop, ")"), side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("> Top ", topRatio, " TADs with highest genes/enhancers ratio:\n"))
tmpSort <- genesEnhancersTAD_DT[order(genesEnhancersTAD_DT$geneEnhancerRatio, decreasing = T),]
tmpSort$geneEnhancerRatio <- round(tmpSort$geneEnhancerRatio, 2)
head(tmpSort$TAD)
stopifnot(tmpSort$TAD %in% names(adj_pvalComb_rank))
tmpSort$tadRank <- adj_pvalComb_rank[tmpSort$TAD]
write.table(tmpSort[seq_len(topRatio),], quote=F, sep="\t", col.names=T, row.names=F)

write.table(tmpSort[seq_len(topRatio),], quote=F, sep="\t", col.names=T, row.names=F, file = topTableFile)
cat(paste0("... written: ", topTableFile, "\n"))
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




