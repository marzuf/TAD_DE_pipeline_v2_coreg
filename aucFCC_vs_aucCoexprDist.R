
startTime <- Sys.time()
cat(paste0("> Rscript aucFCC_vs_aucCoexprDist.R\n"))

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




# Rscript aucFCC_vs_aucCoexprDist.R

options(scipen=100)

buildTable <- T

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

source("set_dataset_colors.R")
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

### HARD CODED
caller <- "TopDom"
script170_name <- "170_score_auc_pval_withShuffle"

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "AUCfcc_vs_AUCcoexprdist")
system(paste0("mkdir -p ", outFold))

# logFile <- file.path(outFold, paste0("auccFCC_vs_aucCoexprDist_logFile.txt"))  
# system(paste0("rm -f ", logFile))

pipOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")
all_datasets <- list.files(pipOutFold)

cat(paste0("# of datasets found: ", length(all_datasets), "\n"))

if(buildTable) {

all_auc <- foreach(curr_dataset = all_datasets) %dopar% {
  
  aucFCC_file <- file.path(pipOutFold, curr_dataset, script170_name, "allratio_auc_pval.Rdata")
  stopifnot(file.exists(aucFCC_file))
  
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_dataset, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  
  list(aucFCC = aucFCC, aucCoexprDist = aucCoexprDist)
  
}
names(all_auc) <- all_datasets

all_auc_FCC <- sapply(all_auc, function(x) x[["aucFCC"]])
all_auc_CoexprDist <- sapply(all_auc, function(x) x[["aucCoexprDist"]])
stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist))

outFile <- file.path(outFold, "all_auc_FCC.Rdata")
save(all_auc_FCC, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, "all_auc_CoexprDist.Rdata")
save(all_auc_CoexprDist, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_auc_FCC.Rdata")
  all_auc_FCC <- eval(parse(text = load(outFile)))
  outFile <- file.path(outFold, "all_auc_CoexprDist.Rdata")
  all_auc_CoexprDist <- eval(parse(text = load(outFile)))
}

##############################################################################################
############################################## SCATTERPLOT ALL DATASETS
##############################################################################################

stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist) )

stopifnot(names(all_auc_FCC) %in% names(dataset_proc_colors) )

curr_colors <- dataset_proc_colors[names(all_auc_FCC)]

outFile <- file.path(outFold, paste0("aucFCC_vs_aucCoexprDist_all_datasets.", plotType))
do.call(plotType, list(outFile, height = myHeight, width =myWidth))
plot(x = all_auc_FCC,
     y = all_auc_CoexprDist,
     xlim = range(all_auc_FCC),#*c(0.95, 1.15),
     ylim = range(all_auc_CoexprDist),#*c(0.95,1.15),
     xlab = "AUC FCC",
     ylab = "AUC coexpr. dist.",
     pch = 16, cex=0.7,
     col = curr_colors,
     main = "AUC coexpr. dist. vs. AUC FCC")
mtext(text = paste0(caller, " - # of datasets = ", length(all_auc_FCC)), side = 3)
text(x = all_auc_FCC,
     y = all_auc_CoexprDist,
     labels = names(all_auc_FCC),
     col = curr_colors,
     pos=3, cex = 0.7)
corTest <- cor.test(all_auc_FCC, all_auc_CoexprDist)
legTxt <- paste0("PCC = ", round(corTest$estimate,2), "\n(p-val =  ", sprintf("%1.2e", corTest$p.value), ")")
legend("topleft", legend = legTxt, bty="n")
# abline(lm(all_auc_CoexprDist ~ all_auc_FCC), col="grey", lty=2)
add_curv_fit(x = all_auc_FCC, y=all_auc_CoexprDist, withR2 = TRUE, R2shiftX = 0.01, R2shiftY = 0, col="grey", lty=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


##############################################################################################
############################################################################################## BARPLOT WITH BOTH AUC RATIOS
##############################################################################################

stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist) )

auc_DT <- data.frame(dataset=names(all_auc_FCC),
                     auc_fcc = all_auc_FCC,
                     auc_coexpr = all_auc_CoexprDist,
                     stringsAsFactors = FALSE)
rownames(auc_DT) <- NULL

auc_DT <- auc_DT[order(auc_DT$auc_fcc, decreasing = TRUE),]

auc_DT_m <- melt(auc_DT, id=c("dataset"))
auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))

stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]

p_AUC <- ggplot(auc_DT_m, aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat="identity")+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("AUC ratios"),
                     breaks = scales::pretty_breaks(n = 5))+
  facet_grid(~dataset, switch="x") +
  scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."))+
  labs(fill  = "AUC ratio") +
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
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
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("AUC ratio"),
                     breaks = scales::pretty_breaks(n = 5)) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."))+
  labs(fill  = "AUC ratio") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
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

outFile <- file.path(outFold, paste0("aucFCC_aucCoexprDist_all_datasets_barplot_orderFCC.", plotType))
ggsave(p_AUC, filename = outFile, height = 6, width=12)
cat(paste0("... written: ", outFile, "\n"))



auc_DT <- auc_DT[order(auc_DT$auc_coexpr, decreasing = TRUE),]
auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))


p_AUC <- ggplot(auc_DT_m, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("AUC ratio"),
                     breaks = scales::pretty_breaks(n = 5)) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."))+
  labs(fill  = "AUC ratio") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
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

outFile <- file.path(outFold, paste0("aucFCC_aucCoexprDist_all_datasets_barplot_orderCoexpr.", plotType))
ggsave(p_AUC, filename = outFile, height = 6, width=12)
cat(paste0("... written: ", outFile, "\n"))



##############################################################################################
############################################################################################## SCATTERPLOT 64 DATASETS (FOR THE POSTER)
##############################################################################################


all_auc_FCC_64 <- all_auc_FCC[names(all_auc_FCC) != "TCGAskcm_skcm_mutCTNNB1" & names(all_auc_FCC) != "TCGAstad_EBVneg_EBVpos"]
all_auc_CoexprDist_64 <- all_auc_CoexprDist[names(all_auc_CoexprDist) != "TCGAskcm_skcm_mutCTNNB1" & names(all_auc_CoexprDist) != "TCGAstad_EBVneg_EBVpos"]

stopifnot(names(all_auc_FCC_64) == names(all_auc_CoexprDist_64) )

stopifnot(names(all_auc_FCC_64) %in% names(dataset_proc_colors) )

curr_colors_64 <- dataset_proc_colors[names(all_auc_FCC_64)]



outFile <- file.path(outFold, paste0("aucFCC_vs_aucCoexprDist_64_datasets.", plotType))
do.call(plotType, list(outFile, height = 6, width =6))
plot(x = all_auc_FCC_64,
     y = all_auc_CoexprDist_64,
     xlim = range(all_auc_FCC_64),#*c(0.95, 1.15),
     ylim = range(all_auc_CoexprDist_64),#*c(0.95,1.15),
     xlab = "AUC FCC",
     ylab = "AUC coexpression",
     pch = 16, cex=1, cex.lab=1.2,
     col = curr_colors_64,
     main = "AUC coexpression vs. AUC FCC")
mtext(text = paste0(caller, " - # of datasets = ", length(all_auc_FCC_64)), side = 3)
# text(x = all_auc_FCC_64,
#      y = all_auc_CoexprDist_64,
#      labels = names(all_auc_FCC_64),
#      pos=3, cex = 0.7)
corTest <- cor.test(all_auc_FCC_64, all_auc_CoexprDist_64)
legTxt <- paste0("PCC = ", round(corTest$estimate,2), "\n(p-val =  ", sprintf("%1.2e", corTest$p.value), ")")
legend("topleft", legend = legTxt, bty="n")
#abline(lm(all_auc_CoexprDist_64 ~ all_auc_FCC_64), col="grey", lty=2)
add_curv_fit(x = all_auc_FCC_64, y=all_auc_CoexprDist_64, withR2 = TRUE, R2shiftX = 0.05, R2shiftY = 0, col="grey", lty=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


