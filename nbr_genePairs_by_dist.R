startTime <- Sys.time()
cat(paste0("> Rscript nbr_genePairs_by_dist.R\n"))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

source(file.path("coreg_utils.R"))

mainDir <- file.path(setDir, "mnt/ed4/marie/scripts", "TAD_DE_pipeline_v2_TopDom")

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- myHeight
plotCex <- 1.2

# Rscript nbr_genePairs_by_dist.R GSE102073_stic_nostic
exprds <- "TCGAkich_norm_kich"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
hicds <- ""
exprds <- args[1]

outFold <- file.path("NBR_GENEPAIRS_BY_DIST", exprds)
dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, paste0(exprds, "_nbrGenePairsByDist_logFile.txt"))
cat("logFile=", logFile, "\n")
if(!SSHFS) file.remove(logFile)
if(SSHFS) logFile <- ""

foldSuffix <- paste0("_40kb")
family <- "hgnc"
familydata <- "family_short"
distLim <- 500

txt <- paste0("!!! HARD-CODED SETTINGS !!!\n")
printAndLog(txt, logFile)
txt <- paste0("... foldSuffix = ",  foldSuffix, "\n")
printAndLog(txt, logFile)
txt <- paste0("... family = ",  family, "\n")
printAndLog(txt, logFile)
txt <- paste0("... familydata = ",  familydata, "\n")
printAndLog(txt, logFile)
txt <- paste0("... distLim = ",  distLim, "\n")
printAndLog(txt, logFile)


sameTADcol <- "darkorange1"
diffTADcol <- "darkslateblue"

sameFamSameTADcol <- "violetred1"
sameFamDiffTADcol <- "lightskyblue"

myxlab <- paste0("Distance between genes (Kb)")
myylab <- paste0("Cumul. # gene pairs")
myTit <- paste0("# gene pairs along dist. sep.")

# retrieve the information from the loess model, so no need to recompute everything !
# (loess model computed in:  Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R)

sameTAD_loessModFile <- file.path(mainDir, "AUC_COEXPRDIST_WITHFAM_SORTNODUP",
                                  # hicds,
                                  paste0(exprds, "_", family), 
                                  paste0(family, "_", familydata), 
                                  "sameTAD_mod.Rdata")
stopifnot(file.exists(sameTAD_loessModFile))
cat("... load sameTAD_mod\n")
sameTAD_mod <- eval(parse(text = load(sameTAD_loessModFile)))

sameFam_sameTAD_loessModFile <- file.path(mainDir, "AUC_COEXPRDIST_WITHFAM_SORTNODUP",
                                  hicds, 
                                  paste0(exprds, "_", family), 
                                  paste0(family, "_", familydata), 
                                  "sameFamSameTAD_mod.Rdata")
stopifnot(file.exists(sameFam_sameTAD_loessModFile))
cat("... load sameFamSameTAD_mod\n")
sameTAD_mod <- eval(parse(text = load(sameFam_sameTAD_loessModFile)))


diffTAD_loessModFile <- file.path(mainDir, "AUC_COEXPRDIST_WITHFAM_SORTNODUP",
                                  hicds, 
                                  paste0(exprds, "_", family), 
                                  paste0(family, "_", familydata), 
                                  "diffTAD_mod.Rdata")
stopifnot(file.exists(diffTAD_loessModFile))
cat("... load diffTAD_mod\n")
diffTAD_mod <- eval(parse(text = load(diffTAD_loessModFile)))

sameFamDiffTAD_loessModFile <- file.path(mainDir, "AUC_COEXPRDIST_WITHFAM_SORTNODUP",
                                  hicds, 
                                  paste0(exprds, "_", family), 
                                  paste0(family, "_", familydata), 
                                  "sameFamDiffTAD_mod.Rdata")
stopifnot(file.exists(sameFamDiffTAD_loessModFile))
cat("... load sameFamDiffTAD_mod\n")
sameFamDiffTAD_mod <- eval(parse(text = load(sameFamDiffTAD_loessModFile)))



all_genePairs_dist <- sort(c(sameTAD_mod$x, diffTAD_mod$x))/1000
all_genePairs_nbr <- seq_along(all_genePairs_dist)

outFile <- file.path(outFold, paste0(exprds, "_sameTAD_diffTAD_nbrGenePairsByDist.", plotType))
sameTAD_genePairs_dist <- sort(sameTAD_mod$x)/1000
sameTAD_genePairs_nbr <- seq_along(sameTAD_genePairs_dist)
diffTAD_genePairs_dist <- sort(diffTAD_mod$x)/1000
diffTAD_genePairs_nbr <- seq_along(diffTAD_genePairs_dist)
mySub <- paste0(exprds)
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     xlim = c(0, distLim),
     ylim = c(0, max(c(sameTAD_genePairs_nbr, diffTAD_genePairs_nbr), na.rm=T)),
     xlab = myxlab,
     ylab = myylab,
     main = myTit,
     cex.axis = plotCex, cex.lab=plotCex,
     type="l")
lines(x = sameTAD_genePairs_dist,
     y = sameTAD_genePairs_nbr,
     col=sameTADcol)
lines(x = diffTAD_genePairs_dist,
      y = diffTAD_genePairs_nbr,
      col=diffTADcol)
mtext(text = mySub, side=3)
legend("topleft",
       c("sameTAD", "diffTAD"),
       lty=1, bty="n",
       col = c(sameTADcol,  diffTADcol)
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, paste0(exprds, "_sameFamSameTAD_sameFamDiffTAD_nbrGenePairsByDist.", plotType))
sameFamSameTAD_genePairs_dist <- sort(sameFamSameTAD_mod$x)/1000
sameFameSameTAD_genePairs_nbr <- seq_along(sameFamSameTAD_genePairs_dist)
sameFamDiffTAD_genePairs_dist <- sort(sameFamDiffTAD_mod$x)/1000
sameFamDiffTAD_genePairs_nbr <- seq_along(sameFamDiffTAD_genePairs_dist)
mySub <- paste0( exprds, " (sameFam)")
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     xlim = c(0, distLim),
     ylim = c(0, max(c(sameFameSameTAD_genePairs_nbr, sameFamDiffTAD_genePairs_nbr), na.rm=T)),
     xlab = myxlab,
     ylab = myylab,
     main = myTit,
     cex.axis = plotCex, cex.lab=plotCex,
     type="l")
lines(x = sameFamSameTAD_genePairs_dist,
      y = sameFameSameTAD_genePairs_nbr,
      col=sameFamSameTADcol)
lines(x = sameFamDiffTAD_genePairs_dist,
      y = sameFamDiffTAD_genePairs_nbr,
      col=sameFamDiffTADcol)
mtext(text = mySub, side=3)
legend("topleft",
       c("sameFamSameTAD","sameFamDiffTAD"),
       lty=1, bty="n",
       col = c(sameFamSameTADcol, sameFamDiffTADcol)
       )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

