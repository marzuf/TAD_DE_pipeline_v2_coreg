startTime <- Sys.time()
cat(paste0("> Rscript build_final_TAD_table.R\n"))

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)


# Rscript build_final_TAD_table.R TCGAcrc_msi_mss 

familyData <- "hgnc"
family_type <- "family_short"

caller = "TopDom"
curr_dataset <- "TCGAcrc_msi_mss"
args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- args[1]
stopifnot(length(args) == 1)

if(curr_dataset != "TCGAcrc_msi_mss") {
  stop(paste0("! not implemented for ", curr_dataset, " (missing enhancer data) !\n"))
}

# plotType <- "png"
# myHeight <- ifelse(plotType=="png", 400, 7)
# myWidth <- ifelse(plotType=="png", 400, 7)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "BUILD_FINAL_TAD_TABLE", paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))

#*********************************************************************
##### LOAD DATA FROM THE PIPELINE -> meanFC, meanCorr, FCC
#*********************************************************************
dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8_name <- "8c_runAllDown"
script11_name <- "11_runEmpPvalCombined"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
pipeline_regionList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_regionList.Rdata"))))
meanLogFC_TAD <- eval(parse(text = load(file.path(dataset_pipDir, script3_name, "all_meanLogFC_TAD.Rdata"))))
meanCorr_TAD <- eval(parse(text = load(file.path(dataset_pipDir, script4_name, "all_meanCorr_TAD.Rdata"))))
fcc_TAD  <- eval(parse(text = load(file.path(dataset_pipDir, script8_name, "all_obs_prodSignedRatio.Rdata"))))
pvalComb_TAD <- eval(parse(text = load(file.path(dataset_pipDir, script11_name, "emp_pval_combined.Rdata"))))
pvalComb_TAD_adj <- p.adjust(pvalComb_TAD, method="BH")
pvalComb_TAD_adj_sort <- sort(pvalComb_TAD_adj)

stopifnot(setequal(pipeline_regionList, names(meanLogFC_TAD)))
stopifnot(setequal(pipeline_regionList, names(meanCorr_TAD)))
stopifnot(setequal(pipeline_regionList, names(fcc_TAD)))
stopifnot(setequal(pipeline_regionList, names(pvalComb_TAD)))

curr_pip_DT <-  foreach(curr_tad = pipeline_regionList, .combine='rbind') %dopar% {
  pvalCombRank <- which(names(pvalComb_TAD_adj_sort) == curr_tad)
    data.frame(
    region = curr_tad,
    meanLogFC_TAD = meanLogFC_TAD[curr_tad],
    meanCorr_TAD = meanCorr_TAD[curr_tad],
    fcc_TAD = fcc_TAD[curr_tad],
    pvalComb_TAD_adj = pvalComb_TAD_adj[curr_tad],
    pvalCombRank = pvalCombRank,
    stringsAsFactors = FALSE
  )
}
rownames(curr_pip_DT) <- NULL
# region meanLogFC_TAD meanCorr_TAD     fcc_TAD pvalComb_TAD_adj pvalCombRank
# 1  chr10_TAD1   -0.19256795   0.12654739  0.39984532      0.044543338         1610
# 2  chr10_TAD2   -0.02026378   0.26506321  0.08854313      0.006870273          881

stopifnot(setequal(pipeline_regionList, curr_pip_DT$region))

#*********************************************************************
##### LOAD DATA FROM PIP SETTINGS -> NBR GENES
#*********************************************************************

# !! source(settingFile) here INSTEAD
# settingFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput/run_settings_", curr_dataset, ".R"))
# stopifnot(file.exists(settingFile))
mainSettingFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/main_settings.R"))
stopifnot(file.exists(mainSettingFile))
source(mainSettingFile)
# gene2tadDTFile <- "../../gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt"
stopifnot(file.exists(gene2tadDT_file))
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)

TADposDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region", "start", "end"), stringsAsFactors = F)

entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
stopifnot(gene2tadDT$entrezID %in% entrezDT$entrezID)
stopifnot(!duplicated(entrezDT$entrezID))

stopifnot(pipeline_geneList %in% gene2tadDT$entrezID)
stopifnot(pipeline_regionList %in% TADposDT$region)

curr_g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList,]
curr_g2t_DT$gene <- unlist(sapply(curr_g2t_DT$entrezID, function(x) entrezDT$symbol[entrezDT$entrezID == x]))

curr_tad_DT <- TADposDT[TADposDT$region %in% pipeline_regionList,]
curr_tad_DT$datasetGenes <- unlist(sapply(curr_tad_DT$region, function(x) {
  xgenes <- paste0(sort(as.character(curr_g2t_DT$gene[as.character(curr_g2t_DT$region) == as.character(x)])), collapse=",")
}))
curr_tad_DT$nGenes <- unlist(sapply(curr_tad_DT$region, function(x) {
  xgenes <- as.character(curr_g2t_DT$gene[as.character(curr_g2t_DT$region) == as.character(x)])
  length(xgenes)
}))
curr_tad_DT$combPvalRank <- unlist(sapply(curr_tad_DT$region, function(x) which(names(pvalComb_TAD_adj_sort) == x)))
curr_tad_DT <- curr_tad_DT[order(curr_tad_DT$combPvalRank),]

# chromo      region    start      end                                                                                                                                            datasetGenes nGenes combPvalRank
# 2045   chr1 chr1_TAD150 89440001 89920000                                                                                                  CCBL2,GBP1,GBP1P1,GBP2,GBP3,GBP4,GBP5,GBP6,GBP7,RBMXL1     10            1
# 591   chr12 chr12_TAD81 54160001 54600000                                                  FLJ12825,HOTAIR,HOXC10,HOXC11,HOXC12,HOXC13,HOXC4,HOXC5,HOXC6,HOXC8,HOXC9,LOC100240735,LOC400043,SMUG1     14            2

stopifnot(setequal(pipeline_regionList, curr_tad_DT$region))

#*********************************************************************
##### PREPARE FAMILY DATA -> nbr families
#*********************************************************************
inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)
# entrezID chromo   start     end     region                                                                              hgnc_family                      hgnc_family_short
# 1    347688  chr10   92828   95178 chr10_TAD1                                                                                 Tubulins                               Tubulins
# 3     10771  chr10  180405  300577 chr10_TAD1 Zinc fingers MYND-type|PHD finger proteins|PWWP domain containing|Bromodomain containing                 Zinc fingers MYND-type
# 8     23185  chr10  852854  977645 chr10_TAD2                                                   La ribonucleoprotein domain containing La ribonucleoprotein domain containing
familyDT$family <- familyDT[,paste0(familyData, "_", family_type)]
stopifnot(!any(is.na(familyDT$family)))
stopifnot(familyDT$family != "")

curr_fam_DT <- foreach(curr_tad = pipeline_regionList, .combine='rbind') %dopar% {
  tad_genes <- curr_g2t_DT$entrezID[curr_g2t_DT$region == curr_tad]
  avFamInfo <- sum(tad_genes %in% familyDT$entrezID)
  naFamInfo <- sum(!tad_genes %in% familyDT$entrezID)
  if(avFamInfo > 0) {
    nUniqFam <- length(unique(familyDT$family[familyDT$entrezID %in% tad_genes]))
  } else{
    nUniqFam <- NA
  }
  data.frame(
    region = curr_tad,
    avFamInfo=avFamInfo,
    naFamInfo=naFamInfo,
    nUniqFam=nUniqFam,
    stringsAsFactors = FALSE
  )
}
# region avFamInfo naFamInfo nUniqFam
# 1  chr10_TAD1         2         2        2
# 2  chr10_TAD2         2         4        2

stopifnot(setequal(pipeline_regionList, curr_fam_DT$region))

#*********************************************************************
##### PREPARE ENHANCER DATA
#*********************************************************************

enhancerTAD_dt <- eval(parse(text = load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_TADS_MAPPED/enhancer_tad_DT.Rdata"))))
rm(enhancer_tad_DT)
enhancerTAD_dt <- na.omit(enhancerTAD_dt)
stopifnot(!any(duplicated(enhancerTAD_dt$enhancer)))
enhancerTAD_dt$region <- as.character(enhancerTAD_dt$region)
# chromo     start       end    enhancer       region
# 1   chr6 133518045 133518161 GH06I133518  chr6_TAD180
# 2  chr13  51630569  51633727 GH13I051630  chr13_TAD58
# 3  chr17  63664756  63666583 GH17I063664 chr17_TAD104
# 4   chr2  23052345  23052647 GH02I023052   chr2_TAD33

curr_enh_DT <- foreach(curr_tad = pipeline_regionList, .combine='rbind') %dopar% {
  nEnh <- nrow(enhancerTAD_dt[enhancerTAD_dt$region == curr_tad,])
  data.frame(
    region = curr_tad,
    nbrGeneHancer = nEnh,
    stringsAsFactors = FALSE
  )
}

# region nbrGeneHancer
# 1  chr10_TAD1           141
# 2  chr10_TAD2            26

stopifnot(setequal(pipeline_regionList, curr_enh_DT$region))

#*********************************************************************
##### ASSEMBLE
#*********************************************************************

stopifnot(setequal(pipeline_regionList, curr_pip_DT$region))

stopifnot(setequal(pipeline_regionList, curr_tad_DT$region))

stopifnot(setequal(pipeline_regionList, curr_fam_DT$region))

stopifnot(setequal(pipeline_regionList, curr_enh_DT$region))


out_dt <- inner_join(curr_pip_DT, curr_tad_DT, by="region")
out_dt <- inner_join(out_dt, curr_fam_DT, by="region")
out_dt <- inner_join(out_dt, curr_enh_DT, by="region")

# stopifnot(!is.na(out_dt))
stopifnot(out_dt$combPvalRank == out_dt$pvalCombRank)

out_dt <- out_dt[,c("region","chromo", "start", "end", "meanLogFC_TAD","meanCorr_TAD","fcc_TAD","pvalComb_TAD_adj", "pvalCombRank","datasetGenes","nGenes","avFamInfo","naFamInfo",        
                    "nUniqFam","nbrGeneHancer")]
out_dt <- out_dt[order(out_dt$pvalCombRank),]
outFile <- file.path(outFold, paste0(curr_dataset,  "_all_TADs_recapTable.txt"))
write.table(out_dt, file = outFile, sep="\t", quote=F, col.names=T, row.names=F)
cat(paste0("... written: ", outFile, "\n"))

#*********************************************************************
#*********************************************************************
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


