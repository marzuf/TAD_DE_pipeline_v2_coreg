library(foreach)
library(doMC)
library(tools)
library(DESeq2)

startTime <- Sys.time()

pointPch <- 16
pointCex <- 1
cexAxis <- 1.1
cexLab <- 1.1

# Rscript gene_variance.R <exprType>
# Rscript gene_variance.R fpkm
# Rscript gene_variance.R log2fpkm
# Rscript gene_variance.R NBvar
# Rscript gene_variance.R voom
cat("> START: gene_variance.R\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
exprType <- "log2fpkm"
nTopLast <- 1000
exprType <- args[1]
if(!is.na(args[2]))
  nTopLast <- as.numeric(args[2])
stopifnot(is.numeric(nTopLast))
stopifnot(exprType %in% c("fpkm", "log2fpkm", "NBvar", "voom"))

exprTypeName <- ifelse(exprType == "fpkm", "FPKM", 
				  ifelse(exprType == "log2fpkm", "log2 FPKM",
					ifelse(exprType == "NBvar", "negative binomial var.",
						ifelse(exprType == "voom", "voom transf.", NA))))
stopifnot(!is.na(exprTypeName))

cat("... START with exprType =\t", exprType, "\n")
cat("... START with nTopLast =\t", nTopLast, "\n")

buildTable <- F

# registerDoMC(ifelse(SSHFS, 2, 20))

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))

settingFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

pipMainFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

aucFCCfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUCfcc_vs_AUCcoexprdist/all_auc_FCC.Rdata")
aucCoexprDistFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUCfcc_vs_AUCcoexprdist/all_auc_CoexprDist.Rdata")
stopifnot(file.exists(aucFCCfile))
stopifnot(file.exists(aucCoexprDistFile))

all_setting_files <- list.files(settingFolder, full.names=T)
all_setting_files <- all_setting_files[grep(".R$", all_setting_files)]
stopifnot(length(all_setting_files) > 0)

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", paste0("GENE_VARIANCE_", toupper(exprType)))
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
# myHeight <- ifelse(plotType == "png", 600, 10)
# myWidth <- ifelse(plotType == "png", 600, 10)
myHeight <- ifelse(plotType == "png", 300, 5)
myWidth <- ifelse(plotType == "png", 300, 5)

if(buildTable) {
  all_ds_geneVarDT <- foreach(ds_file = all_setting_files, .combine="rbind") %dopar% {
  
    stopifnot(file.exists(ds_file))
    cat("... source settingFile", basename(ds_file), "\n")
    source(ds_file)
  
    curr_ds <- basename(pipOutFold)
    cat("... START", curr_ds, "\n")
  
    ds_pipFolder <- file.path(pipMainFolder, pipOutFold)
    cat(ds_pipFolder,"\n")
    stopifnot(file.exists(ds_pipFolder))
    
    cat("... load samp1\n")
    samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
    cat("... load samp2\n")
    samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))
    
    cat("... load geneList\n")
    geneList <- eval(parse(text = load(
      file.path(ds_pipFolder, "0_prepGeneData", "pipeline_geneList.Rdata")
    )))
    
    cat("... load exprDT\n")
    
    if(exprType == "fpkm" | exprType == "log2fpkm" | exprType == "NBvar") {
      exprDT <- eval(parse(text = load(
        file.path(ds_pipFolder, "0_prepGeneData", paste0("rna_", "fpkm", "DT.Rdata"))
      )))
    } else if(exprType == "voom") {
      exprDT <- eval(parse(text = load(
        file.path(ds_pipFolder, "1_runGeneDE", paste0(exprType, "_lmFitInputDT.Rdata"))
      )))
    } else{
      stop(paste0("!!! unimplemented exprType", exprType, "\n"))
    }
    
    stopifnot(samp1 %in% colnames(exprDT))
    stopifnot(samp2 %in% colnames(exprDT))
    stopifnot(names(geneList) %in% rownames(exprDT))
    
    curr_exprDT <- exprDT[rownames(exprDT) %in% names(geneList), c(samp1,samp2)]  
    stopifnot(is.numeric(curr_exprDT[1,1]))
    
    if(exprType == "log2fpkm") {
      curr_exprDT <- log2(curr_exprDT + 1)
      stopifnot(is.numeric(curr_exprDT[1,1]))
    } else if(exprType == "NBvar") {
        cat("!!! WARNING to use \"DESeq2:::varianceStabilizingTransformation\" the expression data are rounded !!!\n")
        curr_exprDT <- try(DESeq2:::varianceStabilizingTransformation(round(as.matrix(curr_exprDT))))
        if(class(curr_exprDT) == "try-error"){
          return(data.frame(
            nTopLast = nTopLast,
            data_type = exprType,
            dataset = curr_ds,
            meanMostVar = NA,
            meanLeastVar = NA,
            stringsAsFactors = FALSE
          ))
        }
          
        stopifnot(is.numeric(curr_exprDT[1,1]))
      }
    
    geneVar <- apply(curr_exprDT, 1,  var, na.rm=T)
    geneVar <- sort(geneVar, decreasing = TRUE)
      
    stopifnot(length(geneVar) >= nTopLast)
  
    mostVariant <- geneVar[1:nTopLast]
    leastVariant <- geneVar[(length(geneVar)-nTopLast+1):length(geneVar)]
    stopifnot(length(mostVariant) == nTopLast)
    stopifnot(length(leastVariant) == nTopLast)
    
    meanMostVar <- mean(mostVariant)
    meanLeastVar <- mean(leastVariant)
    
    data.frame(
      nTopLast = nTopLast,
      data_type = exprType,
      dataset = curr_ds,
      meanMostVar = meanMostVar,
      meanLeastVar = meanLeastVar,
      stringsAsFactors = FALSE
    )
    
  }
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  save(all_ds_geneVarDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  all_ds_geneVarDT <- eval(parse(text = load(outFile)))
}

myTit <- paste0("mean variance top and least ", nTopLast, " genes (", exprType, " data)")

mySubTit <- paste0("(n datasets = ", nrow(all_ds_geneVarDT), ")")

if(exprType == "NBvar") {
  all_ds_geneVarDT <- all_ds_geneVarDT[!is.na(all_ds_geneVarDT$meanMostVar) & !is.na(all_ds_geneVarDT$meanMostVar),]
}
######### plot density

outFile <- file.path(outFold, paste0("density_mostVar_log10_vs_leastVar_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))

plot_multiDens(list(
  meanMostVar = log10(all_ds_geneVarDT$meanMostVar),
  meanLeastVar = log10(all_ds_geneVarDT$meanLeastVar)
),
plotTit = myTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######### scatter plot most var ~ least var

outFile <- file.path(outFold, paste0("mostVar_log10_vs_leastVar_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x =  log10(all_ds_geneVarDT$meanMostVar),
     y =  log10(all_ds_geneVarDT$meanLeastVar),
     pch = 16, cex = 0.7,
     xlab = "mean most var (log10)",
     ylab = "mean least var (log10)",
     main = myTit
       )
text(x = log10(all_ds_geneVarDT$meanMostVar), 
     log10(all_ds_geneVarDT$meanLeastVar),
     labels = all_ds_geneVarDT$dataset,
     pos=3, cex = 0.7)

add_legend_narm(x = log10(all_ds_geneVarDT$meanMostVar), y = log10(all_ds_geneVarDT$meanLeastVar), mypos="topleft", mymet="Pearson")  

mtext(text = mySubTit, side = 3)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

####### scatter plot most var, least var, and ratio

all_auc <- c("FCC", "CoexprDist")

aucFCC <- eval(parse(text = load(aucFCCfile)))
aucCoexprDist <- eval(parse(text = load(aucCoexprDistFile)))

for(auc_type in all_auc) {
  
  cat("... start", auc_type, "\n")
  
  curr_auc <- eval(parse(text = paste0("auc", auc_type)))
  stopifnot(all_ds_geneVarDT$dataset %in% names(curr_auc))
  all_ds_geneVarDT[, auc_type] <- curr_auc[all_ds_geneVarDT$dataset]

  myTit <- paste0(auc_type, ": mean variance least ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_leastVar_log10.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(y =  all_ds_geneVarDT[, auc_type],
       x =  log10(all_ds_geneVarDT$meanLeastVar),
       pch = 16, cex = 0.7,
       ylab = paste0("AUC ", auc_type),
       xlab = "mean least var (log10)",
       main = myTit
  )
  text(x = log10(all_ds_geneVarDT$meanLeastVar),
       y = all_ds_geneVarDT[, auc_type],
       labels = all_ds_geneVarDT$dataset,
       pos=3, cex = 0.7)
  add_legend_narm(x = log10(all_ds_geneVarDT$meanLeastVar), y = all_ds_geneVarDT[, auc_type], mypos="topleft", mymet="Pearson")  
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  myTit <- paste0(auc_type, ": mean variance most ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVar_log10.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(y =  all_ds_geneVarDT[, auc_type],
       x =  log10(all_ds_geneVarDT$meanMostVar),
       pch = 16, cex = 0.7,
       ylab = paste0("AUC ", auc_type),
       xlab = "mean most var (log10)",
       main = myTit
  )
  text(x = log10(all_ds_geneVarDT$meanMostVar),
       y = all_ds_geneVarDT[, auc_type],
       labels = all_ds_geneVarDT$dataset,
       pos=3, cex = 0.7)
  add_legend_narm(x = log10(all_ds_geneVarDT$meanMostVar), y = all_ds_geneVarDT[, auc_type], mypos="topleft", mymet="Pearson")  
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  myTit <- paste0(auc_type, ": mean variance most/least ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVarleastVar_log10_ratio.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(y =  all_ds_geneVarDT[, auc_type],
       x =  log10(all_ds_geneVarDT$meanMostVar)/log10(all_ds_geneVarDT$meanLeastVar),
       pch = 16, cex = 0.7,
       ylab = paste0("AUC ", auc_type),
       xlab = "mean most/mean least var (log10)",
       main = myTit
  )
  text(x =  log10(all_ds_geneVarDT$meanMostVar)/log10(all_ds_geneVarDT$meanLeastVar),
       y = all_ds_geneVarDT[, auc_type],
       labels = all_ds_geneVarDT$dataset,
       pos=3, cex = 0.7)
  add_legend_narm(x =  log10(all_ds_geneVarDT$meanMostVar)/log10(all_ds_geneVarDT$meanLeastVar), y = all_ds_geneVarDT[, auc_type], mypos="topleft", mymet="Pearson")  
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}

for(auc_type in all_auc) {
  
  cat("... start", auc_type, "\n")
  
  curr_auc <- eval(parse(text = paste0("auc", auc_type)))
  stopifnot(all_ds_geneVarDT$dataset %in% names(curr_auc))
  all_ds_geneVarDT[, auc_type] <- curr_auc[all_ds_geneVarDT$dataset]
  
  myTit <- paste0(auc_type, ": mean variance least ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_leastVar_log10_withFit_noLab.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  xvar <- log10(all_ds_geneVarDT$meanLeastVar)
  yvar <- all_ds_geneVarDT[, auc_type]
  
  plot(y =  yvar,
       x =  xvar,
       pch = pointPch, cex = pointCex,
       ylab = paste0("AUC ", auc_type),
       #xlab = "mean least var (log10)",
	   xlab = paste0("Mean var. [log10] top ", nTopLast, " least variant genes\n(", exprTypeName, ")"),
       main = myTit,
cex.axis = cexAxis, cex.lab = cexLab
  )
  # text(x = log10(all_ds_geneVarDT$meanLeastVar),
  #      y = all_ds_geneVarDT[, auc_type],
  #      labels = all_ds_geneVarDT$dataset,
  #      pos=3, cex = 0.7)
  add_legend_narm(x = xvar, y = yvar, mypos="topleft", mymet="Pearson")  
  add_curv_fit(x = xvar, y = yvar, 
               withR2 = TRUE,
               col ="gray80")
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  myTit <- paste0(auc_type, ": mean variance most ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVar_log10_withFit_noLab.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  xvar <- log10(all_ds_geneVarDT$meanMostVar)
  yvar <- all_ds_geneVarDT[, auc_type]
  
  plot(y =  yvar,
       x =  xvar,
       pch = pointPch, cex = pointCex,
       ylab = paste0("AUC ", auc_type),
       #xlab = "mean most var (log10)",
	   xlab = paste0("Mean var. [log10] top ", nTopLast, " most variant genes\n(", exprTypeName, ")"),
       main = myTit,
cex.axis = cexAxis, cex.lab = cexLab
  )
  # text(x = log10(all_ds_geneVarDT$meanMostVar),
  #      y = all_ds_geneVarDT[, auc_type],
  #      labels = all_ds_geneVarDT$dataset,
  #      pos=3, cex = 0.7)
  add_legend_narm(x = xvar, y = yvar, mypos="topleft", mymet="Pearson")  
  add_curv_fit(x = xvar, y = yvar, 
               withR2 = TRUE,
               col ="gray80")
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  myTit <- paste0(auc_type, ": mean variance most/least ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVarleastVar_log10_ratio_withFit_noLab.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  xvar <- log10(all_ds_geneVarDT$meanMostVar)/log10(all_ds_geneVarDT$meanLeastVar)
  yvar <- all_ds_geneVarDT[, auc_type]
  plot(y =  yvar,
       x =  xvar,
       pch = pointPch, cex = pointCex,
       ylab = paste0("AUC ", auc_type),
       #xlab = "mean most/mean least var (log10)",
	   xlab = paste0("Mean most/mean least var. [log10] top ", nTopLast, " most variant genes\n(", exprTypeName, ")"),
       main = myTit,
cex.axis = cexAxis, cex.lab = cexLab
  )
  # text(x =  log10(all_ds_geneVarDT$meanMostVar)/log10(all_ds_geneVarDT$meanLeastVar),
  #      y = all_ds_geneVarDT[, auc_type],
  #      labels = all_ds_geneVarDT$dataset,
  #      pos=3, cex = 0.7)
  add_legend_narm(x = xvar, y = yvar, mypos="topleft", mymet="Pearson")  
  add_curv_fit(x = xvar, y = yvar, 
               withR2 = TRUE,
               col ="gray80")
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
