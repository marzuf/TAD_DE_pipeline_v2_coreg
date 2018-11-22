SSHFS <- F

options(scipen=100)

library(foreach)
library(tools)
library(Hmisc)

# Rscript coexprDist_multiDatasets_rescaled.R all
# Rscript coexprDist_multiDatasets_rescaled.R top 5
# Rscript coexprDist_multiDatasets_rescaled.R bottom 5

setDir <- ifelse(SSHFS, "/media/electron", "")

source("set_dataset_colors.R")
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)

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
myPatt <- "cancer"

topOnly <- TRUE

if(!is.null(myPatt)) {
  scoreDT <- scoreDT[grepl(myPatt, scoreDT$process),]
}

toPlot <- "top"
nTop <- 5

args <- commandArgs(trailingOnly = TRUE)
toPlot <- args[1]  # all, top, bot

stopifnot(toPlot %in% c("top", "bot", "all"))

if(toPlot %in% c("top", "bot")){
  nTop <- as.numeric(args[2])
}
stopifnot(!is.na(nTop))

cat(paste0("... toPlot =\t", toPlot, "\n"))
cat(paste0("... nTop =\t", nTop, "\n"))
cat(paste0("... myPatt =\t", myPatt, "\n"))

sameTADcol <- "darkorange4"
diffTADcol <- "deepskyblue4"

if(toPlot == "all") {
  toplotRankingDs <- scoreDT$dataset
  nTop <- length(toplotRankingDs)
} else if(toPlot == "bot") {
  toplotRankingDs <- scoreDT$dataset[(nrow(scoreDT) - nTop+1):nrow(scoreDT)]
} else if(toPlot == "top") {
  toplotRankingDs <- scoreDT$dataset[1:nTop]
}else{
  stop("--error\n")
}

if(!is.null(myPatt)) {
  patt_ds <- scoreDT$dataset[grepl(myPatt, scoreDT$process)]  
  toplotRankingDs <- toplotRankingDs[toplotRankingDs %in% patt_ds]
}

stopifnot(length(toplotRankingDs) > 0)

mainDir <- file.path(setDir,
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                     "AUC_COEXPRDIST_SORTNODUP/")

allDs <- toplotRankingDs
allDs <- unique(allDs)

caller <- "TopDom"
outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "COEXPRDIST_MULTIDS_PLOT_RESCALED")
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


# str(all_coexpr_vect)

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

diffTAD_dt_toplot <- diffTAD_dt[, colnames(diffTAD_dt) %in% toplotRankingDs]
diffTAD_minLine_toplotDs <- apply(diffTAD_dt_toplot, 1, min, na.rm=T)
diffTAD_maxLine_toplotDs <- apply(diffTAD_dt_toplot, 1, max, na.rm=T)

sameTAD_dt <- all_coexpr_vect[,grepl("_sameTAD_", colnames(all_coexpr_vect))]
colnames(sameTAD_dt) <- gsub("_sameTAD_vect", "", colnames(sameTAD_dt))
sameTAD_dt <- cbind(sameTAD_dt, all_coexpr_vect[,paste0(allDs[1], "_dist_vect")])
sameTAD_dt <- na.omit(sameTAD_dt)
sameTAD_distVect <- sameTAD_dt[,grepl("_dist_vect", colnames(sameTAD_dt))]

sameTAD_dt_toplot <- sameTAD_dt[, colnames(sameTAD_dt) %in% toplotRankingDs]
sameTAD_minLine_toplotDs <- apply(sameTAD_dt_toplot, 1, min, na.rm=T)
sameTAD_maxLine_toplotDs <- apply(sameTAD_dt_toplot, 1, max, na.rm=T)

stopifnot(nrow(sameTAD_dt) == nrow(diffTAD_dt))
stopifnot(nrow(sameTAD_dt) == nrow(all_coexpr_vect))
stopifnot( (ncol(sameTAD_dt) - 1) == nTop )
stopifnot( ncol(sameTAD_dt_toplot) == nTop )
stopifnot( ncol(diffTAD_dt_toplot) == nTop )


stopifnot( length(diffTAD_minLine_toplotDs) == length(diffTAD_maxLine_toplotDs) )
stopifnot( length(sameTAD_minLine_toplotDs) == length(sameTAD_maxLine_toplotDs) )

stopifnot(length(toplotRankingDs) == nTop)  

mytit <- paste0("Coexpr. ~ dist.")

if(is.null(myPatt)){
  mysubtit <- paste0(toTitleCase(toPlot), " datasets (n=", length(toplotRankingDs), ")")
} else {
  mysubtit <- paste0(toTitleCase(toPlot), " datasets (n=", length(toplotRankingDs), "; patt=", myPatt, ")")
}




#==================================================================== START PLOTTING

##############################
### plot same TAD and diff TAD
##############################

all_lines_toplot = c("same", "diff", "same_diff")

sameLegTxt <- "sameTAD"
diffLegTxt <- "diffTAD"
diffLegTxt <-  expression(italic("diffTAD"))

for(line_toplot in all_lines_toplot){
  
  if(line_toplot == "same") {
    xrange <- range(c(sameTAD_distVect))
    yrange <- range(c(sameTAD_minLine_toplotDs, sameTAD_maxLine_toplotDs))
    legTxt <- c(sameLegTxt)
    legCol <- c(sameTADcol)
  } else if(line_toplot == "diff") {
    xrange <- range(c(diffTAD_distVect))
    yrange <- range(c(diffTAD_minLine_toplotDs, diffTAD_maxLine_toplotDs))
    legTxt <- c(diffLegTxt)
    legCol <- c(diffTADcol)
  } else if(line_toplot == "same_diff") {
    xrange <- range(c(diffTAD_distVect, sameTAD_distVect))
    yrange <- range(c(sameTAD_minLine_toplotDs, sameTAD_maxLine_toplotDs,diffTAD_minLine_toplotDs, diffTAD_maxLine_toplotDs))
    legTxt <- c(sameLegTxt, diffLegTxt)
    legCol <- c(sameTADcol, diffTADcol)
  } else {
    stop("--error\n")
  }
  
  if(!is.null(myPatt)){
    outFile <- file.path(outFold, paste0(line_toplot, "TAD_", toPlot, "_n", nTop, "_datasets_", myPatt, "_vectDist.", plotType))
  } else {
    outFile <- file.path(outFold, paste0(line_toplot, "TAD_", toPlot, "_n", nTop, "_datasets_vectDist.", plotType))  
  }
  
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(NULL,
       xlim = xrange,
       ylim = yrange,
       # xlab="", 
       # ylab="",
       xlab=my_xlab,
       ylab=my_ylab,
       # main=paste0("Top and bottom ", nTop, 
       #             " datasets: coexpr ~ dist loess fit")
       main=paste0(mytit)
  )
  mtext(text = mysubtit, side=3)
  
  same_curves <- list()
  diff_curves <- list()
  
  for(i in 1:nTop) {
    same_curves[[i]] <- list(x=sameTAD_distVect, y=sameTAD_dt_toplot[,i])    
    diff_curves[[i]] <- list(x=diffTAD_distVect, y=diffTAD_dt_toplot[,i])
    if(grepl("same", line_toplot)) {
      lines( x = sameTAD_distVect, y = sameTAD_dt_toplot[,i], col = sameTADcol)  
    }
    if(grepl("diff", line_toplot)) {
      lines( x = diffTAD_distVect, y = diffTAD_dt_toplot[,i], col = diffTADcol)  
    }
  }
  
  stopifnot(length(same_curves) == length(diff_curves))
  stopifnot( colnames(sameTAD_dt_toplot) == colnames(diffTAD_dt_toplot) )
  
  curve_labels_same <- colnames(sameTAD_dt_toplot)
  curve_labels_diff <- as.expression(sapply(colnames(diffTAD_dt_toplot), function(y) bquote(italic(.(y)))))
  
  color_labels <- dataset_proc_colors[curve_labels_same]
  stopifnot(!is.na(color_labels))
  
  if(grepl("same", line_toplot)) {
    labcurve(same_curves, curve_labels_same, tilt=TRUE, type="l", col=color_labels)
  } 
  
  if(grepl("diff", line_toplot)) {  
    labcurve(diff_curves, curve_labels_diff, tilt=TRUE, type="s", col=color_labels)
  }
  
  
  legend("topright", 
         legend=legTxt, 
         col = legCol,
         lty=1,
         bty = "n")
  
  fooo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}





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


