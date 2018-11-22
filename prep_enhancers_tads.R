startTime <- Sys.time()

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

genehancer_DT <- read.delim("genehancer.csv", header=T, stringsAsFactors = F, sep=",")
nrow(genehancer_DT)
genehancer_DT <- genehancer_DT[grepl("enhancer", tolower(genehancer_DT$feature.name)),]
nrow(genehancer_DT)

head(genehancer_DT)

source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/main_settings.R"))
source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R"))

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_TADS_MAPPED")
system(paste0("mkdir -p ", outFold))

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)


TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
TAD_DT <- read.table(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
TAD_DT <- TAD_DT[grepl("_TAD", TAD_DT$region),]


# enhancer_tad_DT <- foreach(i = 1:100, .combine="rbind") %do% {
  enhancer_tad_DT <- foreach(i = 1:nrow(genehancer_DT), .combine="rbind") %dopar% {
  cat(paste0("... start enhancer: ", i, "/", nrow(genehancer_DT), "\n"))
  
  curr_att <- genehancer_DT$attributes[i]
  all_att <- strsplit(x=curr_att,split=";")[[1]]
  curr_enhancer <- all_att[grep("genehancer_id", all_att)]
  stopifnot(length(curr_enhancer) == 1)
  curr_enhancer <- gsub("genehancer_id=", "", curr_enhancer)
  
  mapTAD <- TAD_DT$region[ TAD_DT$chromo == genehancer_DT$chrom[i] & 
                            TAD_DT$start <= genehancer_DT$start[i] &
                           TAD_DT$end >= genehancer_DT$end[i]]
  
  
  if(length(mapTAD) == 0) mapTAD <- NA
  
    data.frame(
      chromo = genehancer_DT$chrom[i],
      start = genehancer_DT$start[i],
      end = genehancer_DT$end[i],
      enhancer = curr_enhancer,
      region = mapTAD,
      stringsAsFactors = F
    )
  
}
  

outFile <- file.path(outFold, "enhancer_tad_DT.Rdata")
save(enhancer_tad_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
