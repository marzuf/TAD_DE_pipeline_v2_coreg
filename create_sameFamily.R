startTime <- Sys.time()
cat(paste0("> Rscript create_sameFamily.R\n"))


stop("!!! use sortNoDup !!!\n")

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

caller <- "TopDom"
familyData <- "hgnc"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_SAME_FAMILY")
system(paste0("mkdir -p ", outFold))

inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)


all_familyData <- paste0(familyData, c("_family", "_family_short"))

for(i_fam in all_familyData) {
  
  all_fams <- as.character(unique(familyDT[, i_fam]))
  
  all_family_pairs <- foreach(fam = all_fams, .combine='rbind') %dopar% {
    
    fam_dt <- familyDT[familyDT[,i_fam] == fam,]
    if(nrow(fam_dt) == 1) return(NULL)
    famDT <- as.data.frame(t(combn(fam_dt$entrezID, m=2)))
    colnames(famDT) <- c("gene1", "gene2")
    famDT$family <- fam
    famDT
  }
  
  outFile <- file.path(outFold, paste0(i_fam, "_all_family_pairs.Rdata"))
  save(all_family_pairs, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))


}

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
