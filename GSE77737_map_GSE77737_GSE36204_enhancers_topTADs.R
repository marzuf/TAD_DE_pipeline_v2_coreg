startTime <- Sys.time()
cat(paste0("> Rscript map_GSE77737_GSE36204_enhancers_topTADs.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

# Rscript map_GSE77737_GSE36204_enhancers_topTADs.R TCGAcrc_msi_mss 10

caller = "TopDom"
curr_dataset = "TCGAcrc_msi_mss"
nTopTADs = 10
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
curr_dataset <- args[1]
nTopTADs <- as.numeric(args[2])
stopifnot(!is.na(nTopTADs))

# just to ensure after the source() settingFile
caller_cmdl = caller
curr_dataset_cmdl = curr_dataset 
nTopTADs_cmdl = nTopTADs


if(curr_dataset != "TCGAcrc_msi_mss") {
  stop(paste0("! enhancer data not available for: ", curr_dataset, "!\n"))
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "MAP_GSE77737_GSE36204_ENHANCERS_topTADS", paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))


clAnnotFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg", "h3k27ac_msi_mss_data",
                         "GSE_data", "GSE77737_GSE36401_samples_info.csv")
stopifnot(file.exists(clAnnotFile))

GSE77737_VELsFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg", "h3k27ac_msi_mss_data",
                                 "GSE_data", "GSE77737_VELs")
stopifnot(file.exists(GSE77737_VELsFolder))
GSE77737_VELsFiles <- list.files(GSE77737_VELsFolder, full.names = T, pattern = ".csv$")

GSE36204_VELsFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg", "h3k27ac_msi_mss_data",
                                 "GSE_data", "GSE36204_VELs")
stopifnot(file.exists(GSE36204_VELsFolder))
GSE36204_VELsFiles <- list.files(GSE36204_VELsFolder, full.names = T, pattern = ".txt$")


#*********************************************************************
##### LOAD AND PREP DATA FROM THE PIPELINE -> retrieve topTADs
#*********************************************************************
dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

script0_name <- "0_prepGeneData"
script11_name <- "11_runEmpPvalCombined"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
pipeline_regionList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_regionList.Rdata"))))
pvalComb_TAD <- eval(parse(text = load(file.path(dataset_pipDir, script11_name, "emp_pval_combined.Rdata"))))

pvalComb_TAD_adj <- p.adjust(pvalComb_TAD, method="BH")
pvalComb_TAD_adj_sort <- sort(pvalComb_TAD_adj)

topTADs <- names(pvalComb_TAD_adj_sort)[1:nTopTADs]

#*********************************************************************
##### LOAD AND PREP DATA FROM PIP SETTINGS -> gene and TAD positions
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
stopifnot(topTADs %in% TADposDT$region)
stopifnot(topTADs %in% pipeline_regionList)

# curr_g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList,]
curr_g2t_DT <- gene2tadDT[gene2tadDT$region %in% topTADs,]
curr_g2t_DT$gene <- unlist(sapply(curr_g2t_DT$entrezID, function(x) entrezDT$symbol[entrezDT$entrezID == x]))

curr_tad_DT <- TADposDT[TADposDT$region %in% topTADs,]
curr_tad_DT$datasetGenes <- unlist(sapply(curr_tad_DT$region, function(x) {
  xgenes <- paste0(sort(as.character(curr_g2t_DT$gene[as.character(curr_g2t_DT$region) == as.character(x)])), collapse=",")
}))
curr_tad_DT$combPvalRank <- unlist(sapply(curr_tad_DT$region, function(x) which(names(pvalComb_TAD_adj_sort) == x)))
curr_tad_DT <- curr_tad_DT[order(curr_tad_DT$combPvalRank),]

caller <- caller_cmdl 
curr_dataset <- curr_dataset_cmdl 
nTopTADs <- nTopTADs_cmdl 

#*********************************************************************
##### LOAD CELL LINE DATA
#*********************************************************************
clAnnotDT <- read.delim(clAnnotFile, header=T, stringsAsFactors = FALSE, sep=",")
my_status_cl <- setNames(clAnnotDT$MS_status, clAnnotDT$cell_line)

#*********************************************************************
##### ITERATE OVER ENHANCER DATA
#*********************************************************************
# i_file = GSE77737_VELsFiles[1]
# i_file = GSE77737_VELsFiles[3]
# i_file = GSE77737_VELsFiles[5]

for(i_file in GSE77737_VELsFiles) {
  curr_file <- basename(i_file)
  cat(paste0("... Start GSE77737 file:\t", basename(i_file), "\n"))
  
  curr_dt <- read.delim(i_file, header=T, stringsAsFactors = F, sep=",")
  stopifnot("chr" %in% colnames(curr_dt), "start" %in% colnames(curr_dt), "end" %in% colnames(curr_dt))
  stopifnot(is.numeric(curr_dt$start))
  stopifnot(is.numeric(curr_dt$end))
  
  # select the ones in topTADs
  match_tad <- unlist(sapply(seq_len(nrow(curr_dt)), function(x) {
    enh_chr <- as.character(curr_dt$chr[x])
    enh_start <- curr_dt$start[x]
    enh_end <- curr_dt$end[x]
    
    idx <- which(curr_tad_DT$chromo == enh_chr &
                   curr_tad_DT$start <= enh_start &
                   curr_tad_DT$end >= enh_end)
    if(length(idx) == 0) {
      return(NA)
    } else{
      stopifnot(length(idx) == 1)
      return(curr_tad_DT$region[idx])
    }
  }))
  cat(paste0("...... found enhancers in topTADs:\t", sum(!is.na(match_tad)), "/", nrow(curr_dt), "\n"))
  
  if(all(is.na(match_tad))){
    next
  } else{
    curr_dt$matchTAD <- match_tad
    curr_dt <- curr_dt[!is.na(curr_dt$matchTAD),]
    stopifnot(!any(is.na(curr_dt)))
  }
  
  if(curr_file == "supp2_gainedVELs.csv" | curr_file == "supp3_lostVELs.csv") {
    # binary matrix with cell lines as colnames [chr, start, end, <...>, <cl>, <...>, recurrence]
    stopifnot(any(names(my_status_cl) %in% colnames(curr_dt)))
    av_cl <- colnames(curr_dt)[colnames(curr_dt) %in% names(my_status_cl)]
    av_msi <- av_cl[my_status_cl[av_cl] == "MSI"]
    av_mss <- av_cl[my_status_cl[av_cl] == "MSS"]
    msi_tmpDT <- curr_dt[,av_msi]
    mss_tmpDT <- curr_dt[,av_mss]
    stopifnot(nrow(msi_tmpDT) == nrow(mss_tmpDT))
    msi_mss_DT <- foreach(i_r = seq_len(nrow(msi_tmpDT)), .combine="rbind") %dopar% {
      nMSI_match <- sum(msi_tmpDT[i_r,])
      nMSS_match <- sum(mss_tmpDT[i_r,])
      clMSI_match <- colnames(msi_tmpDT)[which( as.numeric(msi_tmpDT[i_r,]) == 1 )]
      clMSS_match <- colnames(mss_tmpDT)[which( as.numeric(mss_tmpDT[i_r,]) == 1 )]
      stopifnot(length(clMSI_match) == nMSI_match)
      data.frame(nMSI_match = nMSI_match,
                 nMSS_match = nMSS_match,
                 sampleMSI_match = paste0(clMSI_match, collapse=","),
                 sampleMSS_match = paste0(clMSS_match, collapse=","),
                 stringsAsFactors = FALSE
                 )
    }
    out_dt <- curr_dt[,c("chr", "start", "end")]
    colnames(out_dt) <- c("enh_chr", "enh_start", "enh_end")
    out_dt <- cbind(out_dt, msi_mss_DT)
    out_dt$matchTAD <- curr_dt$matchTAD
  } else if(curr_file == "supp6_CRCadenomaVELs.csv" | curr_file == "supp6_CRCspecificVELs.csv"){
    # VEL_chr,VEL_start,VEL_end,genes
    # -> modif to chr, start, end,genes directly in file
    out_dt <- curr_dt[,c("chr", "start", "end", "genes", "matchTAD")]
    colnames(out_dt) <- c("enh_chr", "enh_start", "enh_end", "enh_genes", "matchTAD")
  } else if(curr_file == "supp7_recurrentGainedVELs.csv" | curr_file == "supp8_recurrentLostVELs.csv") {
    # chr,start,end,gene,recurrence
    out_dt <- curr_dt[,c("chr", "start", "end", "gene", "recurrence", "matchTAD")]
    colnames(out_dt) <- c("enh_chr", "enh_start", "enh_end", "enh_gene", "recurrence", "matchTAD")
  } else {
    stop("unknown file\n")
  }
  out_dt <- left_join(out_dt, curr_tad_DT, by =c("matchTAD"="region"))
  stopifnot(out_dt$combPvalRank <= nTopTADs)
  out_dt <- out_dt[order(out_dt$combPvalRank),]
  outFile <- file.path(outFold, paste0(curr_dataset, "_topTADs", nTopTADs, "_", "GSE77737", "_", gsub(".csv", ".txt", curr_file)))
  write.table(out_dt, file = outFile, sep="\t", quote=F, row.names=F, col.names=T)
  cat(paste0("... written: ", outFile, "\n"))
}

# i_file = GSE36204_VELsFiles[1]
# i_file = GSE36204_VELsFiles[2]

for(i_file in GSE36204_VELsFiles) {
  
  curr_file <- basename(i_file)
  cat(paste0("... Start GSE36204 file:\t", basename(i_file), "\n"))
  
  curr_dt <- read.delim(i_file, header=T, stringsAsFactors = F, sep="\t")
  stopifnot("chr" %in% colnames(curr_dt), "start" %in% colnames(curr_dt), "end" %in% colnames(curr_dt))
  # there are some rows that look like
  # sample enhancer_type chr start end
  # 103823   V429       control CHR START END
  curr_dt$start <- as.numeric(curr_dt$start)
  curr_dt$end <- as.numeric(curr_dt$end)
  curr_dt <- na.omit(curr_dt)
  stopifnot(is.numeric(curr_dt$start))
  stopifnot(is.numeric(curr_dt$end))
  
  # select the ones in topTADs
  match_tad <- unlist(sapply(seq_len(nrow(curr_dt)), function(x) {
    enh_chr <- as.character(curr_dt$chr[x])
    enh_start <- curr_dt$start[x]
    enh_end <- curr_dt$end[x]
    
    idx <- which(curr_tad_DT$chromo == enh_chr &
                   curr_tad_DT$start <= enh_start &
                   curr_tad_DT$end >= enh_end)
    if(length(idx) == 0) {
      return(NA)
    } else{
      stopifnot(length(idx) == 1)
      return(curr_tad_DT$region[idx])
    }
  }))
  cat(paste0("...... found enhancers in topTADs:\t", sum(!is.na(match_tad)), "/", nrow(curr_dt), "\n"))
  
  if(all(is.na(match_tad))){
    next
  } else{
    curr_dt$matchTAD <- match_tad
    curr_dt <- curr_dt[!is.na(curr_dt$matchTAD),]
    stopifnot(!any(is.na(curr_dt)))
  }
  if(curr_file == "GSE36204_common_VELs.txt") {
    # chr	start	stop	type_of_VEL	num_samples_shared_by 
    # -> modif stop to end directly in file
    out_dt <- curr_dt
  } else if(curr_file == "GSE36204_all_samples_control_enhancers_and_VELs.txt"){
    # sample	enhancer_type	enh_chr	enh_start	enh_stop # sample == cell line
    # -> modif to sample, enhancer, chr, start, end in file
  
    # select sample
    curr_dt <- curr_dt[curr_dt$sample %in% names(my_status_cl),]
    curr_dt$MS_status <- unlist(sapply(curr_dt$sample, function(x){
      as.character(my_status_cl[as.character(x)])
    }))
    
    msi_names <- names(my_status_cl[my_status_cl == "MSI"])
    mss_names <- names(my_status_cl[my_status_cl == "MSS"])
    
    tmp_msi_DT <- aggregate(sample ~  chr + start + end + enhancer_type + matchTAD, data = curr_dt, function(x){
      msi_s <- x[x %in% msi_names]
      paste0(msi_s, collapse=",")
    }) 
    colnames(tmp_msi_DT)[colnames(tmp_msi_DT) == "sample"] <- "sampleMSI_match"
    
    tmp_mss_DT <- aggregate(sample ~  chr + start + end + enhancer_type + matchTAD, data = curr_dt, function(x){
      mss_s <- x[x %in% mss_names]
      paste0(mss_s, collapse=",")
    }) 
    colnames(tmp_mss_DT)[colnames(tmp_mss_DT) == "sample"] <- "sampleMSS_match"
    
    
    tmp_nMSI_DT <- aggregate(MS_status ~  chr + start + end + enhancer_type + matchTAD, data = curr_dt, function(x) sum(x == "MSI"))
    colnames(tmp_nMSI_DT)[colnames(tmp_nMSI_DT) == "MS_status"] <- "nMSI_match"
    
    tmp_nMSS_DT <- aggregate(MS_status ~  chr + start + end + enhancer_type + matchTAD, data = curr_dt, function(x) sum(x == "MSS"))
    colnames(tmp_nMSS_DT)[colnames(tmp_nMSS_DT) == "MS_status"] <- "nMSS_match"
    
    out_dt <- inner_join(tmp_nMSI_DT, tmp_nMSS_DT, by = c("chr", "start", "end", "enhancer_type", "matchTAD"))
    out_dt <- inner_join(out_dt, tmp_msi_DT, by = c("chr", "start", "end", "enhancer_type", "matchTAD"))
    out_dt <- inner_join(out_dt, tmp_mss_DT, by = c("chr", "start", "end", "enhancer_type", "matchTAD"))
    tmpMatchTAD <- out_dt$matchTAD
    out_dt$matchTAD <- NULL
    out_dt$matchTAD <- tmpMatchTAD
  } else {
    stop("unknown file\n")
  }
  colnames(out_dt)[colnames(out_dt) == "chr"] <- "enh_chr"
  colnames(out_dt)[colnames(out_dt) == "start"] <- "enh_start"
  colnames(out_dt)[colnames(out_dt) == "end"] <- "enh_end"
  out_dt <- left_join(out_dt, curr_tad_DT, by =c("matchTAD"="region"))
  stopifnot(out_dt$combPvalRank <= nTopTADs)
  out_dt <- out_dt[order(out_dt$combPvalRank),]
  outFile <- file.path(outFold, paste0(curr_dataset, "_topTADs", nTopTADs, "_", gsub(".csv", ".txt", curr_file)))
  # outFile <- file.path(outFold, paste0(curr_dataset, "_topTADs", nTopTADs,"_", "GSE36204", "_", gsub(".csv", ".txt", curr_file)))
  write.table(out_dt, file = outFile, sep="\t", quote=F, row.names=F, col.names=T)
  cat(paste0("... written: ", outFile, "\n"))
}


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

