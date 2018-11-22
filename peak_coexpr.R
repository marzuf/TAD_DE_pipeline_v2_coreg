SSHFS <- F

options(scipen=13)

require(dplyr)

setDir <- ifelse(SSHFS, "/media/electron", "")

distMin <- 100*10^3
distMax <- 200*10^3
justRange <- 50*10^3
nTop_inPeak_barplot <- 10
nTop_inPeak_geneCoexpr <- 100

plotType <- "svg"
myWidth <- ifelse(plotType == "png", 400,  7)
myHeight <- ifelse(plotType == "png", 300,  5)

dataType <- "all"
dataType <- "sameTAD"
dataType <- "sameFam"
dataType <- "sameTAD_sameFam"

curr_dataset <- "TCGAstad_EBVpos_EBVneg"
curr_dataset <- "stomach_ebv"

curr_family <- "hgnc"

# Rscript peak_coexpr.R sameTAD_sameFam TCGAstad_EBVneg_EBVpos
# Rscript peak_coexpr.R sameTAD_sameFam TCGAcrc_msi_mss
# Rscript peak_coexpr.R sameTAD_sameFam TCGAstad_EBVneg_EBVpos
# Rscript peak_coexpr.R sameTAD_sameFam GSE102073_stic_nostic
# Rscript peak_coexpr.R sameTAD_sameFam GSE65540_before_after
# Rscript peak_coexpr.R sameTAD_sameFam GSE79209_dysp_nodysp

args <- commandArgs(trailingOnly = TRUE)
dataType <- args[1]
curr_dataset <- args[2]
stopifnot(dataType %in% c("all", "sameTAD","sameFam", "sameTAD_sameFam"))

# settingsTxt <- paste0("distMin", distMin/1000, "kb_distMax", distMax/1000, "kb_justRange", justRange/1000, "kb")
settingsTxt <- paste0(dataType, " - distMin", distMin/1000, "kb - distMax", distMax/1000, "kb - justRange", justRange/1000, "kb")


outFold <- file.path("PEAK_COEXPR",
                     gsub(" - ", "_",  settingsTxt),
                     paste0(curr_dataset)
                    )
system(paste0("mkdir -p ", outFold))

#****************************************************************************************************
#**************************************************************************************************** STAD EBV DATA
#****************************************************************************************************

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/stomach_ebv")
source("../stomach_ebv/stomach_ebv_utils.R")


# DEfile <- file.path(setDir, 
#                     "/mnt/ed4/marie/scripts/stomach_ebv/DE_analysis",
#                     "DE_topTable.Rdata")
# topTableDT <- eval(parse(text = load(file.path(DEfile))))

if(curr_dataset == "stomach_ebv") {
  
  dataFile <- file.path(setDir, 
                        "/mnt/ed4/marie/scripts/stomach_ebv/COEXPR_DIST_v3_PLOT_SORTNODUP_500kb",
                        "hgnc_family_allData_dt.Rdata")

} else {
  dataFile <- file.path(setDir, 
                              "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3_SORTNODUP",
                              paste0(curr_dataset, "_", curr_family),
                              "hgnc_family_allData_dt.Rdata")
}



dataDT <- eval(parse(text = load(file.path(dataFile))))
nrow(dataDT)

if(dataType == "sameTAD") {
  dataDT <- dataDT[dataDT$sameTAD == 1,]   
} else if(dataType == "sameFam") {
  dataDT <- dataDT[dataDT$sameFamily == 1,]   
} else if(dataType == "sameTAD_sameFam") {
  dataDT <- dataDT[dataDT$sameTAD == 1 & dataDT$sameFamily == 1,]   
}  else {
  stopifnot(dataType == "all")
}

beforeDT <- dataDT[dataDT$dist < distMin,]
nrow(beforeDT)

afterDT <- dataDT[dataDT$dist > distMax,]
nrow(afterDT)

peakDT <- dataDT[dataDT$dist >= distMin & dataDT$dist <= distMax,]
nrow(peakDT)

# before-peak-after

dataDT$position <- ifelse(dataDT$dist < distMin, "before", 
                          ifelse(dataDT$dist > distMax, "after", "peak"))

dataDT$position <- factor(dataDT$position, levels = c("before", "peak", "after"))

outFile <- file.path(outFold, paste0("coexpr_density_3catego.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot_multiDens(
  list(before = dataDT$coexpr[dataDT$position == "before"],
       peak = dataDT$coexpr[dataDT$position == "peak"],
       after = dataDT$coexpr[dataDT$position == "after"]),
  plotTit = paste0(curr_dataset),
  my_xlab = paste0("gene coexpression")
)
mtext(text = settingsTxt, side=3)
foo <- dev.off()
checkAndPrint(outFile)


# before-justBefore-peak-justAfter-after
dataDT$position <- ifelse(dataDT$dist < distMin-justRange, "before", 
                          ifelse(dataDT$dist >= distMin-justRange & dataDT$dist < distMin, "justBefore", 
                                 ifelse(dataDT$dist >  distMax & dataDT$dist <= distMax+justRange, "justAfter", 
                                        ifelse(dataDT$dist > distMax+justRange, "after", "peak"))))

dataDT$position <- factor(dataDT$position, levels = c("before", "justBefore", "peak", "justAfter", "after"))



outFile <- file.path(outFold, paste0("coexpr_density_5catego.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot_multiDens(
  list(before = dataDT$coexpr[as.character(dataDT$position) == "before"],
       justBefore = dataDT$coexpr[as.character(dataDT$position) == "justBefore"],
       peak = dataDT$coexpr[as.character(dataDT$position) == "peak"],
       justAfter = dataDT$coexpr[as.character(dataDT$position) == "justAfter"],
       after = dataDT$coexpr[as.character(dataDT$position) == "after"]),
  plotTit = paste0(curr_dataset),
  my_xlab = paste0("gene coexpression")
)
mtext(text=settingsTxt, side=3)
foo <- dev.off()
checkAndPrint(outFile)


outFile <- file.path(outFold, paste0("coexpr_boxplot_5catego.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(dataDT$coexpr ~  dataDT$position,
        main = paste0(curr_dataset)
        )
mtext(text = settingsTxt, side=3)
foo <- dev.off()
checkAndPrint(outFile)

subDT <- dataDT[ as.character(dataDT$position) != "before" & as.character(dataDT$position) != "after", ]
subDT$position <- factor(subDT$position, levels = c("justBefore", "peak", "justAfter"))

outFile <- file.path(outFold, paste0("coexpr_boxplot_3catego.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(subDT$coexpr ~  
          subDT$position,
        main =  paste0(curr_dataset))
mtext(text = settingsTxt, side=3)
foo <- dev.off()
checkAndPrint(outFile)

###### WHAT MAKES THE INCREASE ? MORE PAIRS WITH HIGH EXPRESSION VALUES OR LESS WITH LOW EXPRESSION VALUES ?

if(any(subDT$coexpr > 0 & subDT$coexpr <= 0.5)) {
  outFile <- file.path(outFold, paste0("coexpr_boxplot_coexpr0_05.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  boxplot(subDT[subDT$coexpr > 0 & subDT$coexpr <= 0.5,]$coexpr ~  
            subDT[subDT$coexpr > 0 & subDT$coexpr <= 0.5,]$position,
          main = paste0(curr_dataset)
  )
  mtext(text = paste0(settingsTxt, " - coexpr > 0 & <= 0.5"), side=3)
  foo <- dev.off()
  checkAndPrint(outFile)
}

if(any(subDT$coexpr >  0.5)) {
outFile <- file.path(outFold, paste0("coexpr_boxplot_coexpr_bigger05.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(subDT[subDT$coexpr > 0.5,]$coexpr ~  
          subDT[subDT$coexpr > 0.5,]$position,
        main = paste0(curr_dataset)
        )
mtext(text = paste0(settingsTxt, " - coexpr > 0.5"), side=3)
foo <- dev.off()
checkAndPrint(outFile)
}

if(any(subDT$coexpr >  0)) {
outFile <- file.path(outFold, paste0("coexpr_boxplot_coexpr_pos.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(subDT[subDT$coexpr > 0,]$coexpr ~  
          subDT[subDT$coexpr > 0,]$position,
        main = paste0(curr_dataset)
        )
mtext(text = paste0(settingsTxt, " - coexpr > 0"), side=3)
foo <- dev.off()
checkAndPrint(outFile)
}

if(any(subDT$coexpr <  0)) {
outFile <- file.path(outFold, paste0("coexpr_boxplot_coexpr_neg.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(subDT[subDT$coexpr < 0,]$coexpr ~  
          subDT[subDT$coexpr <0,]$position,
        main = paste0(curr_dataset)
)
mtext(text = paste0(settingsTxt, " - coexpr < 0"), side=3)
foo <- dev.off()
checkAndPrint(outFile)
}

if(any(subDT$coexpr <  -0.3)) {
outFile <- file.path(outFold, paste0("coexpr_boxplot_coexpr_smaller03.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(subDT[subDT$coexpr < -0.3,]$coexpr ~  
          subDT[subDT$coexpr < -0.3,]$position,
        main = paste0(curr_dataset)
        )
mtext(text = paste0(settingsTxt, " - coexpr < -0.3"), side=3)
foo <- dev.off()
checkAndPrint(outFile)
}
# plot_multiDens(
#   list(before = dataDT$dist[as.character(dataDT$position) == "before"],
#   justBefore = dataDT$dist[as.character(dataDT$position) == "justBefore"],
#   peak = dataDT$dist[as.character(dataDT$position) == "peak"],
#   justAfter = dataDT$dist[as.character(dataDT$position) == "justAfter"],
#   after = dataDT$dist[as.character(dataDT$position) == "after"])
# )



#*****************************  1st barplot -> without considering coexpr values
inPeakDT <- dataDT[dataDT$dist >= distMin & dataDT$dist <= distMax,]
inPeakDT_sameFam <- inPeakDT[inPeakDT$sameFamily == 1,]
nrow(inPeakDT)
nrow(inPeakDT_sameFam)
inPeakDT_sameFam <- inPeakDT_sameFam[order(inPeakDT_sameFam$coexpr,decreasing = T),]
head(inPeakDT_sameFam)

caller <- "TopDom"
familyData <- "hgnc"
inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)

inPeakDT_sameFam <- left_join(inPeakDT_sameFam, familyDT[,c("entrezID", "hgnc_family_short")], by = c("gene1" = "entrezID"))
colnames(inPeakDT_sameFam)[colnames(inPeakDT_sameFam) == "hgnc_family_short"] <- "family_gene1"

inPeakDT_sameFam <- left_join(inPeakDT_sameFam, familyDT[,c("entrezID", "hgnc_family_short")], by = c("gene2" = "entrezID"))
colnames(inPeakDT_sameFam)[colnames(inPeakDT_sameFam) == "hgnc_family_short"] <- "family_gene2"

stopifnot(inPeakDT_sameFam$family_gene1 == inPeakDT_sameFam$family_gene2)

inPeakDT_sameFam <- inPeakDT_sameFam[order(inPeakDT_sameFam$coexpr,decreasing = T),]
head(inPeakDT_sameFam)

#first select the top coexpressed genes before looking at their family to which they belong to
# inPeakDT_sameFam <- inPeakDT_sameFam[1:nTop_inPeak_geneCoexpr,]

inPeakDT_sameFam_agg <- data.frame(
  family = names(table(inPeakDT_sameFam$family_gene2)),
  nGenes =as.numeric(table(inPeakDT_sameFam$family_gene2)),
  stringsAsFactors = FALSE
) 
inPeakDT_sameFam_agg <-  inPeakDT_sameFam_agg[order(inPeakDT_sameFam_agg$nGenes,decreasing = T),]
inPeakDT_sameFam_nTop_agg <- inPeakDT_sameFam_agg[1:nTop_inPeak_barplot,]

# barplot(inPeakDT_sameFam_nTop_agg$nGenes,names.arg = inPeakDT_sameFam_nTop_agg$family,las=2)
settingSave <- par()$oma

outFile <- file.path(outFold, paste0("peak_topFamily_barplot.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
par(oma = settingSave + c(0,5,0,0))
barplot(rev(inPeakDT_sameFam_nTop_agg$nGenes),
        names.arg = rev(inPeakDT_sameFam_nTop_agg$family),
        cex.names = 0.8,
        horiz=T,las=2)
title(main = paste0(curr_dataset," - nTop = ", nTop_inPeak_barplot, "/", nrow(inPeakDT_sameFam_agg)))
mtext(paste0("in peak: distMin = ", distMin, " - distMax = ", distMax, "; all genes in peak"), side =3)
foo <- dev.off()
checkAndPrint(outFile)

par(oma = settingSave)

#*****************************  2nd barplot -> with considering coexpr values

inPeakDT <- dataDT[dataDT$dist >= distMin & dataDT$dist <= distMax,]
inPeakDT_sameFam <- inPeakDT[inPeakDT$sameFamily == 1,]
nrow(inPeakDT)
nrow(inPeakDT_sameFam)
# inPeakDT_sameFam <- inPeakDT_sameFam[order(inPeakDT_sameFam$coexpr,decreasing = T),]
head(inPeakDT_sameFam)

caller <- "TopDom"
familyData <- "hgnc"
inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)

inPeakDT_sameFam <- left_join(inPeakDT_sameFam, familyDT[,c("entrezID", "hgnc_family_short")], by = c("gene1" = "entrezID"))
colnames(inPeakDT_sameFam)[colnames(inPeakDT_sameFam) == "hgnc_family_short"] <- "family_gene1"

inPeakDT_sameFam <- left_join(inPeakDT_sameFam, familyDT[,c("entrezID", "hgnc_family_short")], by = c("gene2" = "entrezID"))
colnames(inPeakDT_sameFam)[colnames(inPeakDT_sameFam) == "hgnc_family_short"] <- "family_gene2"

stopifnot(inPeakDT_sameFam$family_gene1 == inPeakDT_sameFam$family_gene2)

inPeakDT_sameFam <- inPeakDT_sameFam[order(inPeakDT_sameFam$coexpr,decreasing = T),]
head(inPeakDT_sameFam)

#first select the top coexpressed genes before looking at their family to which they belong to
inPeakDT_sameFam <- inPeakDT_sameFam[1:nTop_inPeak_geneCoexpr,]

inPeakDT_sameFam_agg <- data.frame(
  family = names(table(inPeakDT_sameFam$family_gene2)),
  nGenes =as.numeric(table(inPeakDT_sameFam$family_gene2)),
  stringsAsFactors = FALSE
) 
inPeakDT_sameFam_agg <-  inPeakDT_sameFam_agg[order(inPeakDT_sameFam_agg$nGenes,decreasing = T),]
inPeakDT_sameFam_nTop_agg <- inPeakDT_sameFam_agg[1:nTop_inPeak_barplot,]

# barplot(inPeakDT_sameFam_nTop_agg$nGenes,names.arg = inPeakDT_sameFam_nTop_agg$family,las=2)
settingSave <- par()$oma

outFile <- file.path(outFold, paste0("peak_topFamily_barplot_topCoexprGenes.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
par(oma = settingSave + c(0,5,0,0))
barplot(rev(inPeakDT_sameFam_nTop_agg$nGenes),
        names.arg = rev(inPeakDT_sameFam_nTop_agg$family),
        cex.names = 0.8,
        horiz=T,las=2)
title(main = paste0(curr_dataset," - nTop = ", nTop_inPeak_barplot, "/", nrow(inPeakDT_sameFam_agg)))
mtext(paste0("in peak: distMin = ", distMin, " - distMax = ", distMax, "; nTopGenes = ", nTop_inPeak_geneCoexpr), side =3)
foo <- dev.off()
checkAndPrint(outFile)

par(oma = settingSave)



#****************************************************************************************************
#**************************************************************************************************** pairwise coexpr. comparison
#****************************************************************************************************

########################################
######################################## CRC MSI MSS VS. STOMACH_EBV
########################################

crcMS_dataFile <- file.path(setDir, 
                            "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3_SORTNODUP/TCGAcrc_msi_mss_hgnc",
                            "hgnc_family_allData_dt.Rdata")

crcMS_dataDT <- eval(parse(text = load(file.path(crcMS_dataFile))))
nrow(crcMS_dataDT)
head(crcMS_dataDT)

stadEBV_dataFile <- file.path(setDir, 
                              "/mnt/ed4/marie/scripts/stomach_ebv/COEXPR_DIST_v3_PLOT_SORTNODUP_500kb",
                              "hgnc_family_allData_dt.Rdata")

stadEBV_dataDT <- eval(parse(text = load(file.path(stadEBV_dataFile))))
nrow(stadEBV_dataDT)
head(stadEBV_dataDT)


crc_stad_dataDT <- inner_join(crcMS_dataDT, stadEBV_dataDT, suffix=c("_crc", "_stad"), by=c("gene1", "gene2"))
stopifnot(crc_stad_dataDT$dist_crc == crc_stad_dataDT$dist_stad)

outFile <- file.path(outFold, paste0("crc_stomachEbv_coexprCmp.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = crc_stad_dataDT$coexpr_crc,
     y = crc_stad_dataDT$coexpr_stad,
     xlab = "coexpr CRC",
     ylab = "coexpr stomach_ebv",
        main = paste0("stomach_ebv vs. CRC"),
     pch = 16, cex = 0.7)
points(x = crc_stad_dataDT$coexpr_crc[crc_stad_dataDT$dist_stad >= distMin & crc_stad_dataDT$dist_stad <= distMax],
       y = crc_stad_dataDT$coexpr_stad[crc_stad_dataDT$dist_stad >= distMin & crc_stad_dataDT$dist_stad <= distMax],
       col="red",
       pch = 16, cex = 0.7)
legend("topleft", bty="n", lty=-1, text.col = "red", legend = "in peak")
foo <- dev.off()
checkAndPrint(outFile)

crc_stad_dataDT$delta_coexpr <- crc_stad_dataDT$coexpr_stad - crc_stad_dataDT$coexpr_crc

outFile <- file.path(outFold, paste0("crc_stomachEbv_deltaCoexprDist.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = crc_stad_dataDT$dist_stad/1000,
     y = crc_stad_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb)",
     ylab = "delta coexpr stomachEBV-CRC",
        main = paste0("delta stomachEBV-CRC vs. dist"),
     pch = 16, cex = 0.7)
abline(v = distMin/1000, col="gray", lty=2 )
abline(v = distMax/1000, col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)

outFile <- file.path(outFold, paste0("crc_stomachEbv_deltaCoexprDist_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = log10(crc_stad_dataDT$dist_stad/1000),
     y = crc_stad_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb) [log10]",
     ylab = "delta coexpr stomachEBV-CRC",
        main = paste0("delta stomachEBV-CRC vs. dist [log10]"),
     pch = 16, cex = 0.7)
abline(v = log10(distMin/1000), col="gray", lty=2 )
abline(v = log10(distMax/1000), col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)


########################################
######################################## STAD EBVNEG EBVPOS VS. STOMACH_EBV
########################################

stadPip_dataFile <- file.path(setDir, 
                            "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3_SORTNODUP/TCGAstad_EBVneg_EBVpos_hgnc",
                            "hgnc_family_allData_dt.Rdata")

stadPip_dataDT <- eval(parse(text = load(file.path(stadPip_dataFile))))
nrow(stadPip_dataDT)
head(stadPip_dataDT)


stad_stadPip_dataDT <- inner_join(stadPip_dataDT, stadEBV_dataDT, suffix=c("_stadPip", "_stad"), by=c("gene1", "gene2"))
stopifnot(stad_stadPip_dataDT$dist_stadPip == stad_stadPip_dataDT$dist_stad)

outFile <- file.path(outFold, paste0("stadPip_stomachEbv_coexprCmp.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = stad_stadPip_dataDT$coexpr_stadPip,
     y = stad_stadPip_dataDT$coexpr_stad,
     xlab = "coexpr stadPip",
     ylab = "coexpr stomach_ebv",
     pch = 16, cex = 0.7)

points(x = stad_stadPip_dataDT$coexpr_stadPip[stad_stadPip_dataDT$dist_stad >= distMin & stad_stadPip_dataDT$dist_stad <= distMax],
       y = stad_stadPip_dataDT$coexpr_stad[stad_stadPip_dataDT$dist_stad >= distMin & stad_stadPip_dataDT$dist_stad <= distMax],
       col="red",
       pch = 16, cex = 0.7)
foo <- dev.off()
checkAndPrint(outFile)

stad_stadPip_dataDT$delta_coexpr <- stad_stadPip_dataDT$coexpr_stad - stad_stadPip_dataDT$coexpr_stadPip

outFile <- file.path(outFold, paste0("stadPip_stomachEbv_deltaCoexprDist.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = stad_stadPip_dataDT$dist_stad/1000,
     y = stad_stadPip_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb)",
     ylab = "delta coexpr stadPip - stomachEBV",
     pch = 16, cex = 0.7)
abline(v = distMin/1000, col="gray", lty=2 )
abline(v = distMax/1000, col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)

outFile <- file.path(outFold, paste0("stadPip_stomachEbv_deltaCoexprDist_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = log10(stad_stadPip_dataDT$dist_stad/1000),
     y = stad_stadPip_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb) [log10]",
     ylab = "delta coexpr stadPip - stomachEBV",
     pch = 16, cex = 0.7)
abline(v = log10(distMin/1000), col="gray", lty=2 )
abline(v = log10(distMax/1000), col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)


########################################
######################################## STAD EBVNEG EBVPOS VS. CRC MSI MSS
########################################


stadPip_crc_dataDT <- inner_join(stadPip_dataDT, crcMS_dataDT, suffix=c("_stadPip", "_crc"), by=c("gene1", "gene2"))
stopifnot(stadPip_crc_dataDT$dist_stadPip == stadPip_crc_dataDT$dist_stad)

outFile <- file.path(outFold, paste0("stadPip_crc_coexprCmp.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = stadPip_crc_dataDT$coexpr_stadPip,
     y = stadPip_crc_dataDT$coexpr_crc,
     xlab = "coexpr stadPip",
     ylab = "coexpr crc",
     pch = 16, cex = 0.7)

points(x = stadPip_crc_dataDT$coexpr_stadPip[stadPip_crc_dataDT$dist_crc >= distMin & stadPip_crc_dataDT$dist_crc <= distMax],
       y = stadPip_crc_dataDT$coexpr_crc[stadPip_crc_dataDT$dist_crc >= distMin & stadPip_crc_dataDT$dist_crc <= distMax],
       col="red",
       pch = 16, cex = 0.7)
foo <- dev.off()
checkAndPrint(outFile)


stadPip_crc_dataDT$delta_coexpr <- stadPip_crc_dataDT$coexpr_crc - stadPip_crc_dataDT$coexpr_stadPip

outFile <- file.path(outFold, paste0("stadPip_crc_deltaCoexprDist_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = stadPip_crc_dataDT$dist_crc/1000,
     y = stadPip_crc_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb)",
     ylab = "delta coexpr stadPip - crc",
     pch = 16, cex = 0.7)
abline(v = distMin/1000, col="gray", lty=2 )
abline(v = distMax/1000, col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)

outFile <- file.path(outFold, paste0("stadPip_crc_deltaCoexprDist_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = log10(stadPip_crc_dataDT$dist_crc/1000),
     y = stadPip_crc_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb) [log10]",
     ylab = "delta coexpr stadPip - crc",
     pch = 16, cex = 0.7)
abline(v = log10(distMin/1000), col="gray", lty=2 )
abline(v = log10(distMax/1000), col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)


########################################
######################################## CRC MSI MSS VS. GSE102073
########################################

crcMS_dataFile <- file.path(setDir, 
                            "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3_SORTNODUP/TCGAcrc_msi_mss_hgnc",
                            "hgnc_family_allData_dt.Rdata")

crcMS_dataDT <- eval(parse(text = load(file.path(crcMS_dataFile))))
nrow(crcMS_dataDT)
head(crcMS_dataDT)

GSE102073_dataFile <- file.path(setDir, 
                              "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3_SORTNODUP/GSE102073_stic_nostic_hgnc",
                              "hgnc_family_allData_dt.Rdata")

GSE102073_dataDT <- eval(parse(text = load(file.path(GSE102073_dataFile))))
nrow(GSE102073_dataDT)
head(GSE102073_dataDT)


crc_GSE102073_dataDT <- inner_join(crcMS_dataDT, GSE102073_dataDT, suffix=c("_crc", "_GSE102073"), by=c("gene1", "gene2"))
stopifnot(crc_GSE102073_dataDT$dist_crc == crc_GSE102073_dataDT$dist_GSE102073)

outFile <- file.path(outFold, paste0("GSE102073_crc_coexprCmp.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = crc_GSE102073_dataDT$coexpr_crc,
     y = crc_GSE102073_dataDT$coexpr_GSE102073,
     xlab = "coexpr crc",
     ylab = "coexpr GSE102073",
     pch = 16, cex = 0.7)
points(x = crc_GSE102073_dataDT$coexpr_crc[crc_GSE102073_dataDT$dist_GSE102073 >= distMin & crc_GSE102073_dataDT$dist_GSE102073 <= distMax],
       y = crc_GSE102073_dataDT$coexpr_GSE102073[crc_GSE102073_dataDT$dist_GSE102073 >= distMin & crc_GSE102073_dataDT$dist_GSE102073 <= distMax],
       col="red",
       pch = 16, cex = 0.7)
foo <- dev.off()
checkAndPrint(outFile)


crc_GSE102073_dataDT$delta_coexpr <- crc_GSE102073_dataDT$coexpr_GSE102073 - crc_GSE102073_dataDT$coexpr_crc

outFile <- file.path(outFold, paste0("GSE102073_crc_deltaCoexprDist.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = crc_GSE102073_dataDT$dist_GSE102073/1000,
     y = crc_GSE102073_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb)",
     ylab = "delta coexpr GSE102073 - crc",
     pch = 16, cex = 0.7)
abline(v = distMin/1000, col="gray", lty=2 )
abline(v = distMax/1000, col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)

outFile <- file.path(outFold, paste0("GSE102073_crc_deltaCoexprDist_log10.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = log10(crc_GSE102073_dataDT$dist_GSE102073/1000),
     y = crc_GSE102073_dataDT$delta_coexpr,
     xlab = "dist. between genes (kb) [log10]",
     ylab = "delta coexpr GSE102073 - crc",
     pch = 16, cex = 0.7)
abline(v = log10(distMin/1000), col="gray", lty=2 )
abline(v = log10(distMax/1000), col="gray", lty=2 )
foo <- dev.off()
checkAndPrint(outFile)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# 
# #****************************************************************************************************
# #**************************************************************************************************** CRC MSI MSS DATA
# #****************************************************************************************************
# 
# 
# # if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/stomach_ebv")
# # source("stomach_ebv_utils.R")
# 
# # DEfile <- file.path(setDir, 
# #                     "/mnt/ed4/marie/scripts/stomach_ebv/DE_analysis",
# #                     "DE_topTable.Rdata")
# # topTableDT <- eval(parse(text = load(file.path(DEfile))))
# 
# 
# crcMS_dataFile <- file.path(setDir, 
#                               "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3_SORTNODUP/TCGAcrc_msi_mss_hgnc",
#                               "hgnc_family_allData_dt.Rdata")
# 
# crcMS_dataDT <- eval(parse(text = load(file.path(crcMS_dataFile))))
# nrow(crcMS_dataDT)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# if(dataType == "sameTAD") {
#   crcMS_dataDT <- crcMS_dataDT[crcMS_dataDT$sameTAD == 1,]   
# }  else if(dataType == "sameFam") {
#   crcMS_dataDT <- crcMS_dataDT[crcMS_dataDT$sameFamily == 1,]   
# }  else if(dataType == "sameTAD_sameFam") {
#   crcMS_dataDT <- crcMS_dataDT[crcMS_dataDT$sameTAD == 1 & crcMS_dataDT$sameFamily == 1,]   
# } else {
#   stopifnot(dataType == "all")
# }
# 
# crcMS_beforeDT <- crcMS_dataDT[crcMS_dataDT$dist < distMin,]
# nrow(crcMS_beforeDT)
# 
# crcMS_afterDT <- crcMS_dataDT[crcMS_dataDT$dist > distMax,]
# nrow(crcMS_afterDT)
# 
# crcMS_peakDT <- crcMS_dataDT[crcMS_dataDT$dist >= distMin & crcMS_dataDT$dist <= distMax,]
# nrow(crcMS_peakDT)
# 
# # before-peak-after
# 
# crcMS_dataDT$position <- ifelse(crcMS_dataDT$dist < distMin, "before", 
#                                   ifelse(crcMS_dataDT$dist > distMax, "after", "peak"))
# 
# crcMS_dataDT$position <- factor(crcMS_dataDT$position, levels = c("before", "peak", "after"))
# 
# plot_multiDens(
#   list(before = crcMS_dataDT$coexpr[crcMS_dataDT$position == "before"],
#        peak = crcMS_dataDT$coexpr[crcMS_dataDT$position == "peak"],
#        after = crcMS_dataDT$coexpr[crcMS_dataDT$position == "after"])
# )
# 
# # before-justBefore-peak-justAfter-after
# crcMS_dataDT$position <- ifelse(crcMS_dataDT$dist < distMin-justRange, "before", 
#                                   ifelse(crcMS_dataDT$dist >= distMin-justRange & crcMS_dataDT$dist < distMin, "justBefore", 
#                                          ifelse(crcMS_dataDT$dist >  distMax & crcMS_dataDT$dist <= distMax+justRange, "justAfter", 
#                                                 ifelse(crcMS_dataDT$dist > distMax+justRange, "after", "peak"))))
# 
# crcMS_dataDT$position <- factor(crcMS_dataDT$position, levels = c("before", "justBefore", "peak", "justAfter", "after"))
# 
# boxplot(crcMS_dataDT$coexpr ~  crcMS_dataDT$position)
# 
# crcMS_subDT <- crcMS_dataDT[ as.character(crcMS_dataDT$position) != "before" & as.character(crcMS_dataDT$position) != "after", ]
# crcMS_subDT$position <- factor(crcMS_subDT$position, levels = c("justBefore", "peak", "justAfter"))
# boxplot(crcMS_subDT$coexpr ~  
#           crcMS_subDT$position)
# 
# ###### WHAT MAKES THE INCREASE ? MORE PAIRS WITH HIGH EXPRESSION VALUES OR LESS WITH LOW EXPRESSION VALUES ?
# 
# boxplot(crcMS_subDT[crcMS_subDT$coexpr > 0 & crcMS_subDT$coexpr <= 0.5,]$coexpr ~  
#           crcMS_subDT[crcMS_subDT$coexpr > 0 & crcMS_subDT$coexpr <= 0.5,]$position)
# 
# boxplot(crcMS_subDT[crcMS_subDT$coexpr > 0.5,]$coexpr ~  
#           crcMS_subDT[crcMS_subDT$coexpr > 0.5,]$position)
# 
# 
# boxplot(crcMS_subDT[crcMS_subDT$coexpr < 0,]$coexpr ~  
#           crcMS_subDT[crcMS_subDT$coexpr <0,]$position)
# 
# boxplot(crcMS_subDT[crcMS_subDT$coexpr < -0.3,]$coexpr ~  
#           crcMS_subDT[crcMS_subDT$coexpr < -0.3,]$position)
# 
# # plot_multiDens(
# #   list(before = crcMS_dataDT$dist[as.character(crcMS_dataDT$position) == "before"],
# #   justBefore = crcMS_dataDT$dist[as.character(crcMS_dataDT$position) == "justBefore"],
# #   peak = crcMS_dataDT$dist[as.character(crcMS_dataDT$position) == "peak"],
# #   justAfter = crcMS_dataDT$dist[as.character(crcMS_dataDT$position) == "justAfter"],
# #   after = crcMS_dataDT$dist[as.character(crcMS_dataDT$position) == "after"])
# # )
# plot_multiDens(
#   list(before = crcMS_dataDT$coexpr[as.character(crcMS_dataDT$position) == "before"],
#        justBefore = crcMS_dataDT$coexpr[as.character(crcMS_dataDT$position) == "justBefore"],
#        peak = crcMS_dataDT$coexpr[as.character(crcMS_dataDT$position) == "peak"],
#        justAfter = crcMS_dataDT$coexpr[as.character(crcMS_dataDT$position) == "justAfter"],
#        after = crcMS_dataDT$coexpr[as.character(crcMS_dataDT$position) == "after"])
# )
# 
# 
# ########################################################################################## FAMILY DATA
# 
# cat(paste0("... load FAMILY data\t", Sys.time(), "\t"))
# caller <- "TopDom"
# i_fam <- "hgnc_family"
# familyFile <- file.path(setDir, 
#             paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_FAMILY_SORTNODUP/", i_fam, "_all_family_pairs.Rdata"))
# familyDT <- eval(parse(text = load(familyFile)))
# cat(paste0(Sys.time(), "\n"))
# head(familyDT)
# nrow(familyDT)
# familyDT$gene1 <- as.character(familyDT$gene1)
# familyDT$gene2 <- as.character(familyDT$gene2)
# stopifnot(familyDT$gene1 < familyDT$gene2)
# 
# 
# crc_peakFamilyDT <- familyDT[familyDT$gene1 %in% c(crcMS_peakDT$gene1, crcMS_peakDT$gene2) &
#                                familyDT$gene2 %in% c(crcMS_peakDT$gene1, crcMS_peakDT$gene2), ]
# nrow(crc_peakFamilyDT)
# length(unique(familyDT$family))
# # 1297
# length(unique(crc_peakFamilyDT$family))
# # 999
# 
# peakFamilyDT <- familyDT[familyDT$gene1 %in% c(peakDT$gene1, peakDT$gene2) &
#                                familyDT$gene2 %in% c(peakDT$gene1, peakDT$gene2), ]
# nrow(peakFamilyDT)
# length(unique(familyDT$family))
# # 1297
# length(unique(peakFamilyDT$family))
# # 1225
# 
# 
# ########################################################################################## CMP COEXPR
# 
# 
# # is the correlation of coexpression CRC vs. STAD stronger for gene pairs in the peak ???
# 
# crcMS_dataDT$genePairs <- paste0(crcMS_dataDT$gene1, "_", crcMS_dataDT$gene2)
# dataDT$genePairs <- paste0(dataDT$gene1, "_", dataDT$gene2)
# 
# crc_genePairs <- crcMS_dataDT$genePairs
# stad_genePairs <- dataDT$genePairs
# 
# 
# interDT <- inner_join(crcMS_dataDT, dataDT, by=c("gene1", "gene2", "genePairs"), suffix = c("_crc", "_stad"))
# 
# stopifnot(interDT$dist_crc == interDT$dist_stad)
# 
# plot(interDT$coexpr_crc ~ interDT$coexpr_stad, pch=16, cex=0.7)
# nrow(interDT)
# cor.test(x=interDT$coexpr_crc, y=interDT$coexpr_stad)
# nrow(interDT[interDT$dist_crc > distMin & interDT$dist_crc < distMax,])
# plot(interDT[interDT$dist_crc > distMin & interDT$dist_crc < distMax,]$coexpr_crc ~ 
#        interDT[interDT$dist_crc > distMin & interDT$dist_crc < distMax,]$coexpr_stad, pch=16, cex=0.7)
# 
# cor.test(x=interDT[interDT$dist_crc > distMin & interDT$dist_crc < distMax,]$coexpr_crc ,
#          y=interDT[interDT$dist_crc > distMin & interDT$dist_crc < distMax,]$coexpr_stad)
# 
# inter_genePairs <- intersect(crc_genePairs, stad_genePairs)
# length(inter_genePairs)
# 
# 
# ########################################################################################## 
# 
# # COEXPR AND RELATIVE POSITION IN tad
# # -> HIGHER COEXPR WHEN ONE OF THE GENE IS AT START OR END ??
# 
# # (midPos-tadStart)/(tadEnd-tadStart)
