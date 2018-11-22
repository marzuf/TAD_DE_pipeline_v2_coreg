all_hist_marks = c("H3K27ac", "H3K4me1")
all_tads = c("chr1_TAD150", "chr6_TAD58", "chr12_TAD81")

hist_mark = "H3K27ac"
curr_tad = "chr12_TAD81"

for(hist_mark in all_hist_marks){
  for(curr_tad in all_tads){
    
    fc_thresh = 2
    
    diff_peak_file = file.path("HIST_PEAKS_LISTRATIO",
                               "TCGAcrc_msi_mss",
                                paste0("histpeaks_list_ratio_logFile_", hist_mark, "_", curr_tad, ".txt"))
    
    reg_feat_file = file.path("REGELEMENTS_HISTPEAKS",
                              "TCGAcrc_msi_mss",
                              paste0(curr_tad, "_overlap_histPeaks_regFeatures.txt"))
    
    outFold = file.path("TOP_PEAKS_ANNOT", "TCGAcrc_msi_mss")
    system(paste0("mkdir -p ", outFold))
    
    
    diffPeak_DT = read.delim(diff_peak_file, header=T, stringsAsFactors = FALSE)
    regFeat_DT = read.delim(reg_feat_file, header=T, stringsAsFactors = FALSE)
    
    diffPeak_FCfilter_DT = diffPeak_DT[diffPeak_DT$ratioBin_FC >= fc_thresh,]
    
    diffPeak_FCfilter_DT = diffPeak_FCfilter_DT[order(diffPeak_FCfilter_DT$ratioBin_FC, decreasing=T),]
    
    
    
    if(nrow(diffPeak_FCfilter_DT) > 0) {
      
      diffPeaks = unique(diffPeak_FCfilter_DT$peak_name)  
      
      peak_regFeat_DT = regFeat_DT[regFeat_DT$peak_name %in% diffPeaks, c("peak_name", "DS", "my_desc")]
      
      unique(peak_regFeat_DT$peak_name  )
      
      peak_regFeat_DT = peak_regFeat_DT[match(diffPeaks[diffPeaks %in% peak_regFeat_DT$peak_name], peak_regFeat_DT$peak_name),]
      
      unique(peak_regFeat_DT$peak_name  )
      
      outFile = file.path(outFold, paste0(curr_tad, "_", hist_mark, "_FCthresh", fc_thresh, ".txt"))
      
      write.table(peak_regFeat_DT, 
                  file = outFile,
                  row.names = F,
                  col.names = T,
                  quote = F,
                  sep="\t")
      
      cat(paste0("... written: ", outFile, "\n"))
      
      # for(i in 1:nrow(diffPeak_FCfilter_DT)) {
      #   curr_peak = diffPeak_FCfilter_DT$peak_name[i]
      #   peak_regFeat_DT = regFeat_DT[regFeat_DT$peak_name == curr_peak,c("DS", "my_desc")]
      # }
      
    }


  }
}
