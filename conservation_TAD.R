require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTableTAD <- TRUE
buildTableCorr <- TRUE
buildTableAllTADs <- TRUE

tolRad <- 2*40000

registerDoMC(ifelse(SSHFS, 2, 40))

cat(paste0("> START with\tbuildTableTAD = ", buildTableTAD, "\n"))
cat(paste0("> START with\tbuildTableCorr = ", buildTableCorr, "\n"))
cat(paste0("> START with\tbuildTableAllTADs = ", buildTableAllTADs, "\n"))
cat(paste0("> START with\ttolRad = ", tolRad, "\n"))

plotType <- "svg"
myHeight <- 6
myWidth <- 12

source("set_dataset_colors.R")
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)



startTime <- Sys.time()

outFold <- "TAD_CONSERV"
system(paste0("mkdir -p ", outFold))

TADassign_file <- file.path(setDir, 
    "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    

mainPipDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER" )

curr_dataset <- "TCGAcrc_msi_mss"
all_datasets <- list.files(mainPipDir)

if(buildTableTAD) {
  all_files <- list.files(file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom/domains_TopDom/"), full.names = T)
  
  all_tissues_TADs_DT <- foreach(i_file = all_files, .combine='rbind') %dopar% {
    curr_dt <- read.delim(i_file, header=F, col.names=c("chromo", "start", "end"))
    curr_dataset <- gsub("(.+)_chr.+", "\\1", basename(i_file))
    curr_dt$dataset <- curr_dataset
    curr_dt[,c("dataset", "chromo", "start", "end")]
  }
  outFile <- file.path(outFold, "all_tissues_TADs_DT.Rdata")
  save(all_tissues_TADs_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  all_tissues_TADs_DT <- eval(parse(text = load(paste0("all_tissues_TADs_DT.Rdata"))))
  
  tads_DT <- read.delim(TADassign_file, header=F, stringsAsFactors = F, col.names = c("chromo", "region", "start", "end"))
  tads_DT <- tads_DT[grepl("_TAD", tads_DT$region),]
  nTissues <- length(unique(all_tissues_TADs_DT$dataset))
  
  # for each TAD, get conservation
  ratioConserv <- foreach(i_tad = 1:nrow(tads_DT), .combine='c') %dopar% {
    curr_chromo <- tads_DT$chromo[i_tad]
    curr_start <- tads_DT$start[i_tad]
    curr_end <- tads_DT$end[i_tad]
    tissues_matchStart <-  all_tissues_TADs_DT$dataset[
      all_tissues_TADs_DT$chromo == curr_chromo &
        abs(all_tissues_TADs_DT$start-curr_start) <= tolRad]
    tissues_matchEnd <-  all_tissues_TADs_DT$dataset[
      all_tissues_TADs_DT$chromo == curr_chromo &
        abs(all_tissues_TADs_DT$end-curr_end) <= tolRad]
    length(intersect(tissues_matchStart, tissues_matchEnd))/nTissues
  }
  tads_DT$ratioConserv <- ratioConserv
  tads_DT_withConserv <- tads_DT
  outFile <- file.path(outFold, "tads_DT_withConserv.Rdata")
  save(tads_DT_withConserv, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "tads_DT_withConserv.Rdata")
}

tads_DT_withConserv <- eval(parse(text = load(outFile)))
tads_conserv <- setNames(tads_DT_withConserv$ratioConserv, tads_DT_withConserv$region)

cat(paste0("range TAD conserv: ", paste0(round(range(tads_conserv),2), collapse=" - "),  "\n"))


if(buildTableCorr) {
  dataset_conserv_corr_DT <- foreach(curr_dataset = all_datasets, .combine='rbind') %dopar% {
    cat(paste0("... start ", curr_dataset, "\n"))
    fcc_file <- file.path( 
      mainPipDir,
      curr_dataset, 
      "8c_runAllDown",
      "all_obs_prodSignedRatio.Rdata")
    if(!file.exists(fcc_file))
      return(NULL)
    #cat(fcc_file, "\n")
    tads_fcc <- eval(parse(text=load(fcc_file)))
    my_tads_fcc <- intersect(names(tads_fcc), tads_DT_withConserv$region)
    # plot(tads_fcc[my_tads] ~ tads_conserv[my_tads])
    corTest_fcc <- cor.test(tads_conserv[my_tads_fcc], tads_fcc[my_tads_fcc], method = "pearson")
    
    rd_file <- file.path( 
      mainPipDir,
      curr_dataset, 
      "8c_runAllDown",
      "all_obs_ratioDown.Rdata")
    if(!file.exists(rd_file))
      return(NULL)
    #cat(rd_file, "\n")
    tads_rd <- eval(parse(text=load(rd_file)))
    my_tads_rd <- intersect(names(tads_rd), tads_DT_withConserv$region)
    # plot(tads_rd[my_tads] ~ tads_conserv[my_tads])
    corTest_rd <- cor.test(tads_conserv[my_tads_rd], tads_rd[my_tads_rd], method = "pearson")
    
    data.frame(dataset = curr_dataset,
               fcc_PCC = corTest_fcc$estimate,
               fcc_pval = corTest_fcc$p.value,
               rd_PCC = corTest_rd$estimate,
               rd_pval = corTest_rd$p.value,
               stringsAsFactors=F
    )
  }
  
  outFile <- file.path(outFold, "dataset_conserv_corr_DT.Rdata")
  save(dataset_conserv_corr_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFold, "dataset_conserv_corr_DT.Rdata")
}

dataset_conserv_corr_DT <- eval(parse(text = load(outFile)))
rownames(dataset_conserv_corr_DT) <- NULL

dataset_conserv_corr_DT <- dataset_conserv_corr_DT[order(dataset_conserv_corr_DT$fcc_PCC),]

dataset_conserv_corr_DT_m <- melt(dataset_conserv_corr_DT, id=c("dataset")) 


dataset_conserv_corr_DT_m$dataset <- factor(as.character(dataset_conserv_corr_DT_m$dataset), levels = as.character(dataset_conserv_corr_DT$dataset))

stopifnot(as.character(dataset_conserv_corr_DT_m$dataset)  %in% names(dataset_proc_colors) )
curr_colors <- dataset_proc_colors[as.character(levels(dataset_conserv_corr_DT_m$dataset))]


p_Corr <- ggplot(dataset_conserv_corr_DT_m[dataset_conserv_corr_DT_m$variable %in% c("fcc_PCC", "rd_PCC"),], aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("PCC"),
                     breaks = scales::pretty_breaks(n = 5)) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(fcc_PCC = "dodgerblue4", rd_PCC = "darkorange2"),
                    labels = c(fcc_PCC = "FCC", rd_PCC = "ratioDown"))+
  labs(fill  = "Pearson's corr. coeff.") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 1.5, 1), "lines"),
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
  )

if(SSHFS) p_Corr


outFile <- file.path(outFold, paste0("FCC_rD_PCC_with_ratioConserv.", plotType))
ggsave(p_Corr, filename = outFile, height = myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################################################################

  
if(buildTableAllTADs) {
  all_tad_conserv_DT <- foreach(curr_dataset = all_datasets, .combine='rbind') %dopar% {
    cat(paste0("... start ", curr_dataset, "\n"))
    fcc_file <- file.path( 
      mainPipDir,
      curr_dataset, 
      "8c_runAllDown",
      "all_obs_prodSignedRatio.Rdata")
    if(!file.exists(fcc_file))
      return(NULL)
    #cat(fcc_file, "\n")
    tads_fcc <- eval(parse(text=load(fcc_file)))
    
    
    rd_file <- file.path( 
      mainPipDir,
      curr_dataset, 
      "8c_runAllDown",
      "all_obs_ratioDown.Rdata")
    if(!file.exists(rd_file))
      return(NULL)
    #cat(rd_file, "\n")
    tads_rd <- eval(parse(text=load(rd_file)))
    
    meanCorr_file <-  file.path( 
      mainPipDir,
      curr_dataset, 
      "4_runMeanTADCorr",
      "all_meanCorr_TAD.Rdata")
    if(!file.exists(meanCorr_file))
      return(NULL)
    #cat(rd_file, "\n")
    tads_meanCorr <- eval(parse(text=load(meanCorr_file)))
    
    meanFC_file <-  file.path( 
      mainPipDir,
      curr_dataset, 
      "3_runMeanTADLogFC",
      "all_meanLogFC_TAD.Rdata")
    if(!file.exists(meanFC_file))
      return(NULL)
    #cat(rd_file, "\n")
    tads_meanFC <- eval(parse(text=load(meanFC_file)))
    
    
    
    data.frame(dataset = curr_dataset,
              TAD = names(tads_conserv),
               TAD_conserv = tads_conserv,
               TAD_FCC = as.numeric(tads_fcc[names(tads_conserv)]),
               TAD_rD = as.numeric(tads_rd[names(tads_conserv)]),
              TAD_meanFC = as.numeric(tads_meanFC[names(tads_conserv)]),
              TAD_meanCorr = as.numeric(tads_meanCorr[names(tads_conserv)]),
               stringsAsFactors=F
    )
  }
  
  outFile <- file.path(outFold, "all_tad_conserv_DT.Rdata")
  save(all_tad_conserv_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
} else {
  outFile <- file.path(outFold, "all_tad_conserv_DT.Rdata")
  
}

all_tad_conserv_DT <- eval(parse(text = load(outFile)))

all_tad_conserv_DT_noDataset <- all_tad_conserv_DT
all_tad_conserv_DT_noDataset$dataset <- NULL
mean_tad_conserv_DT <- aggregate(. ~  TAD, data=all_tad_conserv_DT_noDataset, FUN=mean, na.rm=T)


densplot <- function(x,y, pch=19, cex=1, ...){
  df <- data.frame(x,y)
  d <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(d)[1,] + 1L
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  df <- df[order(df$dens),]
  plot(df$x,df$y, pch=pch, col=df$col, ...)
}


outFile <- file.path(outFold, paste0("FCC_ratioConserv_all_datasets_meanTADs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myHeight))
densplot(x = mean_tad_conserv_DT$TAD_FCC, 
         y = mean_tad_conserv_DT$TAD_conserv, cex = 0.7, 
         xlab = "TAD FCC (mean across datasets)",
         ylab ="TAD conserv. ratio")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, paste0("rD_ratioConserv_all_datasets_meanTADs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myHeight))
densplot(x = mean_tad_conserv_DT$TAD_rD, 
         y = mean_tad_conserv_DT$TAD_conserv, cex = 0.7, 
         xlab = "TAD ratioDown (mean across datasets)",
         ylab ="TAD conserv. ratio")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



outFile <- file.path(outFold, paste0("meanFC_ratioConserv_all_datasets_meanTADs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myHeight))
densplot(x = mean_tad_conserv_DT$TAD_meanFC, 
         y = mean_tad_conserv_DT$TAD_conserv, cex = 0.7, 
         xlab = "TAD meanFC (mean across datasets)",
         ylab ="TAD conserv. ratio")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, paste0("meanCorr_ratioConserv_all_datasets_meanTADs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myHeight))
densplot(x = mean_tad_conserv_DT$TAD_meanCorr, 
         y = mean_tad_conserv_DT$TAD_conserv, cex = 0.7, 
         xlab = "TAD meanCorr (mean across datasets)",
         ylab ="TAD conserv. ratio")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

mean_tad_conserv_DT$TAD_meanFC_transf <- log10(abs(mean_tad_conserv_DT$TAD_meanFC) + 0.1)
outFile <- file.path(outFold, paste0("meanFCtransf_ratioConserv_all_datasets_meanTADs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width=myHeight))
densplot(x = mean_tad_conserv_DT$TAD_meanFC_transf, 
         y = mean_tad_conserv_DT$TAD_conserv, cex = 0.7, 
         xlab = "TAD log10(abs(meanFC)+0.1) (mean across datasets)",
         ylab ="TAD conserv. ratio")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



###=============================================================================================================================================
###=============================================================================================================================================
###=============================================================================================================================================




###===============================================
###===============================================
###===============================================

cat(paste0(startTime, "\n", Sys.time(), "\n"))







