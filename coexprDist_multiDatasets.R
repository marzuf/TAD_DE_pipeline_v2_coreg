SSHFS <- F

options(scipen=100)

library(foreach)

# Rscript coexprDist_multiDatasets.R

setDir <- ifelse(SSHFS, "/media/electron", "")

# setDir <- "~/Desktop/17.09/media/electron"

scoreDT_file <- file.path(setDir,
          "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
          "BUILDDT_OUTPUT_FOLDER/",
          "17ipreRatios_build_score_table_limited_pval_withShuffle",
          "all_scores_ratios_samples_DT.Rdata")


stopifnot(file.exists(scoreDT_file))


my_xlab <- "Distance between gene pairs (bp)"
my_ylab <- "Gene pairs expression (qqnorm) PCC"
# densityPolygon <- 90
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 350, 8)


scoreDT <- eval(parse(text = load(scoreDT_file)))
scoreDT <- scoreDT[order(scoreDT$prodSignedRatio_auc_permGenes, decreasing = T),]

myPatt <- NULL
# myPatt <- "cancer"

topOnly <- TRUE

if(!is.null(myPatt)) {
  scoreDT <- scoreDT[grepl(myPatt, scoreDT$process),]
}

nTop <- 5
# nTop <- nrow(scoreDT)


topRankingDs <- scoreDT$dataset[1:nTop]
botRankingDs <- scoreDT$dataset[(nrow(scoreDT) - nTop+1):nrow(scoreDT)]

colSameTAD_bot <- "darkorange4"
colSameTAD_top <- "darkorange1"

colDiffTAD_bot <- "deepskyblue4"
colDiffTAD_top <- "deepskyblue1"


mainDir <- file.path(setDir,
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                     "AUC_COEXPRDIST_SORTNODUP/")

allDs <- c(topRankingDs, botRankingDs)
allDs <- unique(allDs)

caller <- "TopDom"
outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "COEXPRDIST_MULTIDS_PLOT")
system(paste0("mkdir -p ", outFold))

cat(paste0("... will draw lines for ", length(allDs) , " datasets\n"))


cat("... start building DT for the datasets\n")
all_coexpr_vect <- foreach(ds = allDs, .combine='cbind') %dopar% {
  
  dist_vect_file <- file.path(mainDir,ds,
                         "distVect.Rdata")  
  stopifnot(file.exists(dist_vect_file))
  dist_vect <- eval(parse(text = load(dist_vect_file)))
  
  
  sameTADvect_file <- file.path(mainDir,ds,
                              "smooth_vals_sameTAD_distVect.Rdata")  
  stopifnot(file.exists(sameTADvect_file))
  sameTAD_vect <- eval(parse(text = load(sameTADvect_file)))
  
  
  diffTADvect_file <- file.path(mainDir,ds,
                                "smooth_vals_diffTAD_distVect.Rdata")  
  stopifnot(file.exists(diffTADvect_file))
  diffTAD_vect <- eval(parse(text = load(diffTADvect_file)))
  
  
  data.frame(dist_vect = dist_vect,
              sameTAD_vect = sameTAD_vect,
              diffTAD_vect = diffTAD_vect,
              stringsAsFactors = FALSE)
 
}

colnames(all_coexpr_vect) <- paste0(
                                    rep(allDs, each=3),
                                    "_",
                                    colnames(all_coexpr_vect))


str(all_coexpr_vect)

check_dist <- all_coexpr_vect[,grepl("dist_vect", colnames(all_coexpr_vect))]

stopifnot(nrow(unique(t(check_dist))) == 1)


cat("... start preparing values to plot\n")

distVect <- unique(t(check_dist))[1,]
stopifnot(length(distVect) == nrow(check_dist))

all_coexpr_vect <- na.omit(all_coexpr_vect)

diffTAD_dt <- all_coexpr_vect[,grepl("_diffTAD_", colnames(all_coexpr_vect))]
colnames(diffTAD_dt) <- gsub("_diffTAD_vect", "", colnames(diffTAD_dt))
diffTAD_dt <- cbind(diffTAD_dt, all_coexpr_vect[,paste0(allDs[1], "_dist_vect")])
diffTAD_dt <- na.omit(diffTAD_dt)
diffTAD_distVect <- diffTAD_dt[,grepl("_dist_vect", colnames(diffTAD_dt))]

diffTAD_dt_top <- diffTAD_dt[, colnames(diffTAD_dt) %in% topRankingDs]
diffTAD_minLine_topDs <- apply(diffTAD_dt_top, 1, min, na.rm=T)
diffTAD_maxLine_topDs <- apply(diffTAD_dt_top, 1, max, na.rm=T)

sameTAD_dt <- all_coexpr_vect[,grepl("_sameTAD_", colnames(all_coexpr_vect))]
colnames(sameTAD_dt) <- gsub("_sameTAD_vect", "", colnames(sameTAD_dt))
sameTAD_dt <- cbind(sameTAD_dt, all_coexpr_vect[,paste0(allDs[1], "_dist_vect")])
sameTAD_dt <- na.omit(sameTAD_dt)
sameTAD_distVect <- sameTAD_dt[,grepl("_dist_vect", colnames(sameTAD_dt))]

sameTAD_dt_top <- sameTAD_dt[, colnames(sameTAD_dt) %in% topRankingDs]
sameTAD_minLine_topDs <- apply(sameTAD_dt_top, 1, min, na.rm=T)
sameTAD_maxLine_topDs <- apply(sameTAD_dt_top, 1, max, na.rm=T)

if(length(topRankingDs) < nrow(scoreDT) & ! topOnly) {
  
  diffTAD_dt_bot <- diffTAD_dt[, colnames(diffTAD_dt) %in% botRankingDs]
  diffTAD_minLine_botDs <- apply(diffTAD_dt_bot, 1, min, na.rm=T)
  diffTAD_maxLine_botDs <- apply(diffTAD_dt_bot, 1, max, na.rm=T)
  
  sameTAD_dt_bot <- sameTAD_dt[, colnames(sameTAD_dt) %in% botRankingDs]
  sameTAD_minLine_botDs <- apply(sameTAD_dt_bot, 1, min, na.rm=T)
  sameTAD_maxLine_botDs <- apply(sameTAD_dt_bot, 1, max, na.rm=T)
  
  mytit <-paste0("Top and bottom ", nTop, 
                          " datasets: coexpr ~ dist loess fit")
  
  yrange <- range(c(diffTAD_minLine_topDs, diffTAD_maxLine_topDs, sameTAD_minLine_botDs, sameTAD_maxLine_botDs))
  
} else if (topOnly) {
  mytit <-paste0("Top ", nTop, 
                 " datasets: coexpr ~ dist loess fit")
  yrange <- range(c(diffTAD_minLine_topDs, diffTAD_maxLine_topDs))
} else {
  mytit <-paste0("All datasets: coexpr ~ dist loess fit")
  yrange <- range(c(diffTAD_minLine_topDs, diffTAD_maxLine_topDs))
}

if(!is.null(myPatt)) {
  outFile <- file.path(outFold, paste0("sameTAD_diffTAD_top", nTop, myPatt,"datasets_vectDist.", plotType))
  
} else{
  outFile <- file.path(outFold, paste0("sameTAD_diffTAD_top", nTop, "datasets_vectDist.", plotType))  
}

do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(NULL,
     xlim = range(c(diffTAD_distVect, sameTAD_distVect)), 
     ylim = yrange,
     # xlab="", 
     # ylab="",
     xlab=my_xlab,
     ylab=my_ylab,
     # main=paste0("Top and bottom ", nTop, 
     #             " datasets: coexpr ~ dist loess fit")
     main=paste0(mytit)
     )



nDatasets <- length(allDs)

if(length(topRankingDs) < nrow(scoreDT)) {
  if(topOnly){
    nDatasets <- nDatasets/2
  }
}

if(nDatasets > nrow(scoreDT))
  nDatasets <- nrow(scoreDT)

  
if(!is.null(myPatt)) {
  mtext(text = paste0("Only ", myPatt, " datasets (n = ", nDatasets, ")") , side = 3)
} else {
  mtext(text = paste0("n = ", nDatasets, "") , side = 3)  
}

for(i in 1:nTop) {
  
  lines( x = sameTAD_distVect, y = sameTAD_dt_top[,i], col = colSameTAD_top)  
  lines( x = diffTAD_distVect, y = diffTAD_dt_top[,i], col = colDiffTAD_top)  
  
  if(length(topRankingDs) < nrow(scoreDT) & !topOnly) {
    lines( x = sameTAD_distVect, y = sameTAD_dt_bot[,i], col = colSameTAD_bot)
    lines( x = diffTAD_distVect, y = diffTAD_dt_bot[,i], col = colDiffTAD_bot)
    
  }
}


if(length(topRankingDs) < nrow(scoreDT) & !topOnly ) {
  legend("topright", 
         legend=c("sameTAD_top", "sameTAD_bot", "diffTAD_top", "diffTAD_bot"),
         col = c(colSameTAD_top, colSameTAD_bot, colDiffTAD_top, colDiffTAD_bot),
         lty=1,
         bty = "n")
} else {
  legend("topright", 
         legend=c("sameTAD", "diffTAD"),
         col = c(colSameTAD_top, colDiffTAD_top),
         lty=1,
         bty = "n")
}

fooo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




# #### TOP DATASETS - SAME TADs
# polygon(x = c(sameTAD_distVect, rev(sameTAD_distVect)),
#         y = c(sameTAD_minLine_topDs,
#               rev(sameTAD_maxLine_topDs)),
#         col = colSameTAD_top,
#         density=densityPolygon)
# #### TOP DATASETS - DIFF TADs
# polygon(x = c(diffTAD_distVect, rev(diffTAD_distVect)),
#         y = c(diffTAD_minLine_topDs,
#               rev(diffTAD_maxLine_topDs)),
#         col = colDiffTAD_top,
#         density=densityPolygon)
# 
# #### BOT DATASETS - SAME TADs
# polygon(x = c(sameTAD_distVect, rev(sameTAD_distVect)),
#         y = c(sameTAD_minLine_botDs,
#               rev(sameTAD_maxLine_botDs)),
#         col = colSameTAD_bot,
#         density=densityPolygon)
# #### BOT DATASETS - DIFF TADs
# polygon(x = c(diffTAD_distVect, rev(diffTAD_distVect)),
#         y = c(diffTAD_minLine_botDs,
#               rev(diffTAD_maxLine_botDs)),
#         col = colDiffTAD_bot,
#         density=densityPolygon)


