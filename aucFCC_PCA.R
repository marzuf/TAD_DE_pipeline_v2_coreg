
startTime <- Sys.time()
cat(paste0("> Rscript aucFCC_PCA.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggrepel, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggthemes, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript aucFCC_PCA.R

options(scipen=100)

buildTable <- TRUE

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 500, 8)
myWidth <- ifelse(plotType == "png", 500, 8)

source("set_dataset_colors.R")
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

### HARD CODED
caller <- "TopDom"
script8_name <- "8c_runAllDown"

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "AUCFCC_PCA")
system(paste0("mkdir -p ", outFold))

# logFile <- file.path(outFold, paste0("auccFCC_vs_aucCoexprDist_logFile.txt"))  
# system(paste0("rm -f ", logFile))

pipOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")
all_datasets <- list.files(pipOutFold)

cat(paste0("# of datasets found: ", length(all_datasets), "\n"))

if(buildTable) {

all_ratios <- foreach(curr_dataset = all_datasets) %dopar% {
  
  fcc_file <- file.path(pipOutFold, curr_dataset, script8_name, "all_obs_prodSignedRatio.Rdata")
  stopifnot(file.exists(fcc_file))
  tads_fcc <- eval(parse(text = load(fcc_file)))
  
  rd_file <- file.path(pipOutFold, curr_dataset, script8_name, "all_obs_ratioDown.Rdata")
  stopifnot(file.exists(rd_file))
  tads_rd <- eval(parse(text = load(rd_file)))
  
  dataset_tads <- intersect(names(tads_rd), names(tads_fcc))
  
  list(tads_fcc = tads_fcc, tads_rd = tads_rd, tads = dataset_tads)
  
}
names(all_ratios) <- all_datasets


outFile <- file.path(outFold, "all_ratios.Rdata")
save(all_ratios, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(outFold, "all_ratios.Rdata")
  all_ratios <- eval(parse(text = load(outFile)))
}


plot_pca <- function(prcomp_obj, PC_x = 1, PC_y = 2, ...) {
  var_exp <- prcomp_obj$sdev^2/sum(prcomp_obj$sdev^2)
  PC_x <- as.numeric(as.character(PC_x))
  PC_y <- as.numeric(as.character(PC_y))
  stopifnot(!is.na(PC_x))
  stopifnot(!is.na(PC_y))
  xlab <- paste0("PC", PC_x, " (", round(var_exp[PC_x]*100,2), "%)")
  ylab <- paste0("PC", PC_y, " (", round(var_exp[PC_y]*100,2), "%)")
  stopifnot(length(var_exp)  >= PC_x)
  stopifnot(length(var_exp)  >= PC_y)
  xpos <- prcomp_obj$x[,PC_x]
  ypos <- prcomp_obj$x[,PC_y]
  stopifnot(length(xpos) == length(ypos))
  
  plot(x=xpos,
       y=ypos,
       xlab=xlab,
       ylab=ylab,
       ...)
}


# retrieve all common tads
common_tads <- Reduce(intersect, lapply(all_ratios, function(x) as.character(x[["tads"]])))


iterate_ratios <- c("FCC", "ratioDown")
curr_ratio="FCC"

pcs_to_plot <- list(c(1,2), c(1,3))

for(curr_ratio in iterate_ratios) {
  
  cat("... START ", curr_ratio, "\n")
  
  for(pcs in pcs_to_plot) {
    stopifnot(length(pcs) == 2) 
    nPCa <- pcs[[1]]
    nPCb <- pcs[[2]]
    stopifnot(is.numeric(nPCa))
    stopifnot(is.numeric(nPCb))
    pcA <- paste0("PC", nPCa) 
    pcB <- paste0("PC", nPCb)
    cat("...... start plotting: ",pcA," and ", pcB, "\n")
    
    ##################
    ### plot the curr_ratio
    ##################
    if(curr_ratio == "FCC") {
      dataset_ratio <- lapply(all_ratios, function(x) x[[paste0("tads_", tolower(curr_ratio))]])  
    }else if(curr_ratio=="ratioDown") {
      dataset_ratio <- lapply(all_ratios, function(x) x[[paste0("tads_rd")]])    
    }
    
    dataset_ratio_filtered <- lapply( dataset_ratio, function(x) x[names(x) %in% common_tads])
    
    ratio_DT <- as.data.frame(dataset_ratio_filtered)
    
    pca_ratio <- prcomp(t(ratio_DT))
    
    mytit <- paste0("PCA - TADs ", curr_ratio)
    mysubtit <- paste0("# intersect TADs = ", length(common_tads), " - # datasets = ", nrow(pca_ratio$x))
    
    outFile <- file.path(outFold, paste0("PC", nPCa,"_PC", nPCb, "_", curr_ratio, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width = myWidth))
    plot_pca(prcomp_obj = pca_ratio, 
             PC_x=nPCa, 
             PC_y=nPCb,
             pch=16, 
             cex =0.9, 
             col = dataset_proc_colors[rownames(pca_ratio$x)])
    title(mytit)
    mtext(text = mysubtit, side=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    ##################
    ### with label
    ##################"
    outFile <- file.path(outFold, paste0("PC", nPCa, "_PC", nPCb, "_", curr_ratio, "_withLabels.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width = myWidth))
    plot_pca(prcomp_obj = pca_ratio, 
             PC_x=nPCa, 
             PC_y=nPCb,
             pch=16, 
             cex =0.9, 
             col = dataset_proc_colors[rownames(pca_ratio$x)])
    text(x = pca_ratio$x[,nPCa], y = pca_ratio$x[,nPCb], labels = rownames(pca_ratio$x), cex = 0.7)
    title(mytit)
    mtext(text = mysubtit, side=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ##################
    ### with label - ggrepel
    ##################
    
    dat <- as.data.frame(pca_ratio$x[,c(nPCa, nPCb)])
    dat$dataset <- rownames(dat)
    
    ggcols <-  dataset_proc_colors[dat$dataset]
    
    p <- ggplot(dat, aes_string(x=paste0(pcA), y=paste0(pcB), label = "dataset")) +
      ggtitle(label = mytit, subtitle = mysubtit)+
      geom_point(color = ggcols) + 
      geom_text_repel(size = 3) +
      theme_base()
    
    outFile <- file.path(outFold, paste0("PC", nPCa, "_PC", nPCb, "_", curr_ratio, "_withLabels_ggrepel.", plotType))
    ggsave(filename = outFile, plot = p, height = myHeight,width=myWidth)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
   
}


# # plot the FCC
# dataset_fcc <- lapply(all_ratios, function(x) x[["tads_fcc"]])  
# dataset_fcc_filtered <- lapply( dataset_fcc, function(x) x[names(x) %in% common_tads])
# 
# fcc_DT <- as.data.frame(dataset_fcc_filtered)
# 
# pca_fcc <- prcomp(t(fcc_DT))
# 
# outFile <- file.path(outFold, paste0("PC1_PC2_FCC.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_fcc, 
#          PC_x=1, 
#          PC_y=2,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# title("PCA - TADs FCC")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# outFile <- file.path(outFold, paste0("PC1_PC3_FCC.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_fcc, 
#          PC_x=1, 
#          PC_y=3,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# title("PCA - TADs FCC")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# 
# # with label
# outFile <- file.path(outFold, paste0("PC1_PC2_FCC_withLabels.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_fcc, 
#          PC_x=1, 
#          PC_y=2,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# text(x = pca_fcc$x[,1], y = pca_fcc$x[,2], labels = rownames(pca_fcc$x), cex = 0.7)
# title("PCA - TADs FCC")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# outFile <- file.path(outFold, paste0("PC1_PC3_FCC_withLabels.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_fcc, 
#          PC_x=1, 
#          PC_y=3,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# text(x = pca_fcc$x[,1], y = pca_fcc$x[,3], labels = rownames(pca_fcc$x), cex = 0.7)
# title("PCA - TADs FCC")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# #=================================================================================================================
# # plot the ratioDown
# dataset_rd <- lapply(all_ratios, function(x) x[["tads_rd"]])  
# dataset_rd_filtered <- lapply( dataset_rd, function(x) x[names(x) %in% common_tads])
# 
# rd_DT <- as.data.frame(dataset_rd_filtered)
# 
# pca_rd <- prcomp(t(rd_DT))
# 
# outFile <- file.path(outFold, paste0("PC1_PC2_ratioDown.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_rd, 
#          PC_x=1, 
#          PC_y=2,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# 
# title("PCA - TADs ratioDown")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# outFile <- file.path(outFold, paste0("PC1_PC3_ratioDown.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_rd, 
#          PC_x=1, 
#          PC_y=3,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# 
# title("PCA - TADs ratioDown")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# outFile <- file.path(outFold, paste0("PC1_PC2_ratioDown_withLabels.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_rd, 
#          PC_x=1, 
#          PC_y=2,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# title("PCA - TADs ratioDown")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# # with label
# text(x = pca_rd$x[,1], y = pca_rd$x[,2], labels = rownames(pca_rd$x), cex = 0.7)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# outFile <- file.path(outFold, paste0("PC1_PC3_ratioDown_withLabels.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width = myWidth))
# plot_pca(prcomp_obj = pca_rd, 
#          PC_x=1, 
#          PC_y=3,
#          pch=16, 
#          cex =0.9, 
#          col = dataset_proc_colors[rownames(pca_fcc$x)])
# title("PCA - TADs ratioDown")
# mtext(text = paste0("# intersect TADs = ", length(common_tads)), side=3)
# # with label
# text(x = pca_rd$x[,1], y = pca_rd$x[,3], labels = rownames(pca_rd$x), cex = 0.7)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))



#=======================================================================================
#=======================================================================================
#=======================================================================================
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
