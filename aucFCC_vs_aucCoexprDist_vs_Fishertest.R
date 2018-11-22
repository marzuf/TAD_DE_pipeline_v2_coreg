
startTime <- Sys.time()
cat(paste0("> Rscript aucFCC_vs_aucCoexprDist_vs_Fishertest.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

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


# Rscript aucFCC_vs_aucCoexprDist_vs_Fishertest.R

options(scipen=100)

buildTable <- F

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

coexprThresh <- 0

source("set_dataset_colors.R")
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

### HARD CODED
caller <- "TopDom"
script170_name <- "170_score_auc_pval_withShuffle"

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), paste0("AUCfcc_vs_AUCcoexprdist_vs_FT_WITHBARPLOT_", coexprThresh))
system(paste0("mkdir -p ", outFold))

# logFile <- file.path(outFold, paste0("auccFCC_vs_aucCoexprDist_logFile.txt"))  
# system(paste0("rm -f ", logFile))

pipOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")
all_datasets <- list.files(pipOutFold)

fisherFold <- file.path(setDir, 
                        paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/FISHER_SAMETAD_COEXPR_", coexprThresh))

cat(paste0("# of datasets found: ", length(all_datasets), "\n"))

if(buildTable) {

all_auc <- foreach(curr_dataset = all_datasets) %dopar% {
  
  aucFCC_file <- file.path(pipOutFold, curr_dataset, script170_name, "allratio_auc_pval.Rdata")
  stopifnot(file.exists(aucFCC_file))
  
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_dataset, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  
  ft_file <- file.path(fisherFold, curr_dataset, "resultsFT.Rdata")
  stopifnot(file.exists(ft_file))
  
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  
  ft_results <-  eval(parse(text = load(ft_file)))
  
  list(aucFCC = aucFCC, aucCoexprDist = aucCoexprDist, 
       pvalFT = as.numeric(ft_results["p_value"]), oddsFT = as.numeric(ft_results["odds_ratio"]))
  
}
names(all_auc) <- all_datasets

all_auc_FCC <- sapply(all_auc, function(x) x[["aucFCC"]])
all_auc_CoexprDist <- sapply(all_auc, function(x) x[["aucCoexprDist"]])
stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist))


all_ft_pval <- sapply(all_auc, function(x) x[["pvalFT"]])
all_ft_odds <- sapply(all_auc, function(x) x[["oddsFT"]])
stopifnot(names(all_auc_FCC) == names(all_ft_pval))
stopifnot(names(all_auc_FCC) == names(all_ft_odds))

outFile <- file.path(outFold, "all_auc_FCC.Rdata")
save(all_auc_FCC, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, "all_auc_CoexprDist.Rdata")
save(all_auc_CoexprDist, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "all_ft_pval.Rdata")
save(all_ft_pval, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "all_ft_odds.Rdata")
save(all_ft_odds, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(outFold, "all_auc_FCC.Rdata")
  all_auc_FCC <- eval(parse(text = load(outFile)))
  outFile <- file.path(outFold, "all_auc_CoexprDist.Rdata")
  all_auc_CoexprDist <- eval(parse(text = load(outFile)))
  outFile <- file.path(outFold, "all_ft_odds.Rdata")
  all_ft_odds <- eval(parse(text = load(outFile)))
  outFile <- file.path(outFold, "all_ft_pval.Rdata")
  all_ft_pval <- eval(parse(text = load(outFile)))
}

##############################################################################################
############################################## SCATTERPLOT ALL DATASETS
##############################################################################################

stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist) )
stopifnot(names(all_auc_FCC) == names(all_ft_pval) )
stopifnot(names(all_auc_FCC) == names(all_ft_odds) )
stopifnot(names(all_auc_FCC) %in% names(dataset_proc_colors) )

#curr_colors <- dataset_proc_colors[names(all_auc_FCC)]

plotRatio <- c("auc_FCC", "auc_CoexprDist", "ft_pval", "ft_odds")

ratioName <- c("auc_FCC" = "AUC FCC",
                "auc_CoexprDist" = "AUC coexpr. dist.", 
               "ft_pval" = "FT p-val", 
               "ft_odds"= "FT odds ratio")

all_ft_pval <- -log10(all_ft_pval)

for(i1 in 1:(length(plotRatio)-1)) {

  curr_ratio1 <- plotRatio[i1]  
  ratioVect1 <- eval(parse(text = paste0("all_", curr_ratio1)))

  
  auc_DT <- data.frame(dataset=names(ratioVect1),
                       auc_fcc = ratioVect1,
                       stringsAsFactors = FALSE)
  rownames(auc_DT) <- NULL
  
  auc_DT <- auc_DT[order(auc_DT$auc_fcc, decreasing = TRUE),]
  
  auc_DT_m <- melt(auc_DT, id=c("dataset"))
  auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))
  
  stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]
  
  
  
  p_AUC <- ggplot(auc_DT_m, aes(x = dataset, y = value, fill = variable)) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(label="", subtitle=paste0("(coexpr. ", coexprThresh, ")")) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 5)) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
    # facet_grid(~dataset, switch="x") +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = c(auc_fcc = "dodgerblue4"),
                      labels = c(auc_fcc = paste0(ratioName[curr_ratio1]) ))+
    labs(fill  = "") +
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.margin = unit(c(1, 1, 2, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(size = 6),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    )+
    geom_hline(yintercept = 1, linetype = 2)
  
  if(SSHFS) p_AUC
  
  outFile <- file.path(outFold, paste0(curr_ratio1, "_all_datasets_barplot_order", curr_ratio1, ".", plotType))
  ggsave(p_AUC, filename = outFile, height = 6, width=12)
  cat(paste0("... written: ", outFile, "\n"))








  for(i2 in i1:length(plotRatio)) {  


  curr_ratio2 <- plotRatio[i2] 
  

  ratioVect2 <- eval(parse(text = paste0("all_", curr_ratio2))) # all_auc_FCC

  curr_colors <- dataset_proc_colors[names(ratioVect2)]
  
  outFile <- file.path(outFold, paste0(curr_ratio2, "_vs_", curr_ratio1, "_all_datasets.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width =myHeight)) # square
  plot(x = ratioVect2,
       y = ratioVect1,
       xlim = range(ratioVect2),#*c(0.95, 1.15),
       ylim = range(ratioVect1),#*c(0.95,1.15),
       xlab = ratioName[curr_ratio2],
       ylab = ratioName[curr_ratio1],
       pch = 16, cex=0.7,
       col = curr_colors,
       main = paste0(ratioName[curr_ratio1], " vs. ", ratioName[curr_ratio2]))
  mtext(text = paste0(caller, " - # of datasets = ", length(ratioVect2), " - coexpr. ", coexprThresh), side = 3)
  text(x = ratioVect2,
       y = ratioVect1,
       labels = names(ratioVect2),
       col = curr_colors,
       pos=3, cex = 0.7)
  corTest <- cor.test(ratioVect2, ratioVect1)
  legTxt <- paste0("PCC = ", round(corTest$estimate,2), "\n(p-val =  ", sprintf("%1.2e", corTest$p.value), ")")
  legend("topleft", legend = legTxt, bty="n")
  # abline(lm(ratioVect1 ~ ratioVect2), col="grey", lty=2)
  add_curv_fit(x = ratioVect2, y=ratioVect1, withR2 = TRUE, R2shiftX = 0.01, R2shiftY = 0, col="grey", lty=2)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ##############################################################################################
  ############################################################################################## BARPLOT WITH BOTH AUC RATIOS
  ##############################################################################################
  
  stopifnot(names(ratioVect2) == names(ratioVect1) )
  
  auc_DT <- data.frame(dataset=names(ratioVect2),
                       auc_fcc = ratioVect2,
                       auc_coexpr = ratioVect1,
                       stringsAsFactors = FALSE)
  rownames(auc_DT) <- NULL
  
  auc_DT <- auc_DT[order(auc_DT$auc_fcc, decreasing = TRUE),]
  
  auc_DT_m <- melt(auc_DT, id=c("dataset"))
  auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))
  
  stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]
  
  p_AUC <- ggplot(auc_DT_m, aes(x = variable, y = value, fill = variable)) +
    ggtitle(label="",subtitle=paste0("(coexpr. ", coexprThresh, ")")) +
    geom_bar(stat="identity")+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 5))+
    facet_grid(~dataset, switch="x") +
    scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                      labels = c(auc_fcc = paste0(ratioName[curr_ratio2]), auc_coexpr = paste0(ratioName[curr_ratio1])))+
    labs(fill  = "") +
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      #plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.margin = unit(c(1, 1, 2, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      strip.text.x = element_text(size = 6, color = curr_colors),
      axis.text.x = element_blank(),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    )
  
  if(SSHFS) p_AUC
  
  
  p_AUC <- ggplot(auc_DT_m, aes(x = dataset, y = value, fill = variable)) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(label="", subtitle=paste0("(coexpr. ", coexprThresh, ")")) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 5)) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
    # facet_grid(~dataset, switch="x") +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                      labels = c(auc_fcc = paste0(ratioName[curr_ratio2]), auc_coexpr = paste0(ratioName[curr_ratio1])))+
    labs(fill  = "") +
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.margin = unit(c(1, 1, 2, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(size = 6),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    )+
    geom_hline(yintercept = 1, linetype = 2)
  
  if(SSHFS) p_AUC
  
  outFile <- file.path(outFold, paste0(curr_ratio2, "_", curr_ratio1, "_all_datasets_barplot_order", curr_ratio2, ".", plotType))
  ggsave(p_AUC, filename = outFile, height = 6, width=12)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  auc_DT <- auc_DT[order(auc_DT$auc_coexpr, decreasing = TRUE),]
  auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))
  
  
  p_AUC <- ggplot(auc_DT_m, aes(x = dataset, y = value, fill = variable)) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle(label="", subtitle=paste0("(coexpr. ", coexprThresh, ")")) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(""),
                       breaks = scales::pretty_breaks(n = 5)) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
    # facet_grid(~dataset, switch="x") +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                      labels = c(auc_fcc = paste0(ratioName[curr_ratio2]), auc_coexpr = paste0(ratioName[curr_ratio1])))+
    labs(fill  = "") +
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.margin = unit(c(1, 1, 2, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(size = 6),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    )+
    geom_hline(yintercept = 1, linetype = 2)
  
  if(SSHFS) p_AUC
  
  outFile <- file.path(outFold, paste0(curr_ratio2, "_", curr_ratio1, "_all_datasets_barplot_orderCoexpr.", plotType))
  ggsave(p_AUC, filename = outFile, height = 6, width=12)
  cat(paste0("... written: ", outFile, "\n"))  
}
}




