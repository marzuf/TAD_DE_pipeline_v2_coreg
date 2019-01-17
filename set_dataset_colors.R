library(RColorBrewer)

annotFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg", "dataset_annot.csv")

cat(paste0("... set datasetcolors using annotFile:\t", annotFile, "\n"))

score_DT <- read.delim(annotFile, header=T, stringsAsFactors = F, sep=",")
score_DT <- score_DT[,c("prefix", "dataset", "process", "cond1", "cond2", "prodSignedRatio_auc_permGenes")]

score_DT <- score_DT[order(score_DT$prodSignedRatio_auc_permGenes),]


score_DT$label <- paste0(score_DT$cond1, " vs.\n", score_DT$cond2)

shorten_process <- c("cancer", "infection", "psychiatric disorder", "muscular disease", "heart disease",
                     "lung disease", "neurologic disorder")

score_DT$process_short <- score_DT$process

for(i in shorten_process) {
  score_DT$process_short[grepl(i, score_DT$process_short)] <- i  
}
score_DT$process_short[grepl("lung disorder", score_DT$process_short)] <- "lung disease"

# score_DT$process[grepl("cancer", score_DT$process)] <- "cancer"

other_processes <- names(table(score_DT$process_short))[as.numeric(table(score_DT$process_short)) <= 3 ]
score_DT$process_short[score_DT$process_short %in% other_processes] <- "other"

table(score_DT$process_short)
# cancer       differentiation embryonic development             infection                 other  psychiatric disorder 
# 27                     8                     4                     8                    13                     5 

my_colors <- c(cancer = "green4",
               differentiation= "chocolate1",
               infection = "red",
               `psychiatric disorder`= "blue3", #"dodgerblue",
               `embryonic development` = "yellow",
               other = "orchid1"
)

#my_colors <- c(cancer = "yellow1",
#               `psychiatric disorder`= "wheat", #"dodgerblue",
#               infection = "snow2",
#               differentiation= "slategray1",
#               other = "thistle1",
#               `embryonic development` = "lightgreen"    
#)



all_proc_short <- unique(score_DT$process_short)

# score_DT$proc_col <- unlist(sapply(score_DT$process_short, function(x)
#   brewer.pal(length(all_proc_short), "Dark2")[which(all_proc_short == x)]
#   ))
score_DT$proc_col <- unlist(sapply(score_DT$process_short, function(x)
  my_colors[x]
))



dataset_colors <- c(noTCGA_cancerDS = "darkorange", TCGA_cancerDS = "firebrick2", no_cancerDS="dodgerblue2")
