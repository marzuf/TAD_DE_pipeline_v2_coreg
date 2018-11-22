histMark <- "H3K4me1"

# HARD CODED MAIN FOLDER
mainFold <- file.path(setDir, 
                      "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017")

refPeakFile <- file.path(mainFold, histMark, paste0(histMark, "_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"))
stopifnot(file.exists(refPeakFile))

refPeaks_DT <- read.delim(refPeakFile, stringsAsFactors = FALSE, col.names=c("peak_chromo", "peak_start", "peak_end", "peak_name"), header=FALSE)
# check the peaks are not overlapping
peak_ends <- refPeaks_DT$end[-length(refPeaks_DT$end)]
peak_starts <- refPeaks_DT$start[-1]
stopifnot(peak_ends <= peak_starts)

refPeaks_DT$size <- refPeaks_DT$peak_end - refPeaks_DT$peak_start + 1

outFile <- file.path(mainFold, histMark, paste0(histMark, "_all_peaks_filtered_merged_named_sizeDensity.svg"))
svg(outFile, height=7, width=8)
plot(density(log10(refPeaks_DT$size)), main = paste0(histMark, " - merged peak size [log10]"))
foo <- dev.off()

outFile  <- file.path(mainFold, histMark, paste0(histMark, "_all_peaks_filtered_merged_named_sizeStat.txt"))
sink(outFile)
print(summary(refPeaks_DT$size))
print(nrow(refPeaks_DT))
sink()