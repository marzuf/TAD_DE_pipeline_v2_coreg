SSFHS <- TRUE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path(setDir, "/mnt/ed4/marie/scripts/SV_project/SV_utils.R"))

obsTSS <- eval(parse(
  text = load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/TSS_DIST_TAD_ALL_GENES/all_TAD_seqDist_DT.Rdata"))))

obs_mean_dist_dna <- obsTSS$mean_dist_dna
obs_mean_dist_alignment <- obsTSS$mean_dist_alignment

randomGenesTSS <- eval(parse(
  text = load(file.path(setDir, "//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/TSS_DIST_TAD_RANDOM_GENES_ALL_GENES/randomizeGenes_randomList_all_TAD_seqDist_DT.Rdata"))))

nRandomGenes <- length(randomGenesTSS)

randomGenes_dist_dna <- unlist(lapply(randomGenesTSS, function(x) x$mean_dist_dna))
randomGenes_dist_alignment <- unlist(lapply(randomGenesTSS, function(x) x$mean_dist_alignment))

randomTADsTSS <- eval(parse(
  text = load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/TSS_DIST_TAD_RANDOM_TADs_ALL_GENES_nRandom2/randomizeTADs_randomList_all_TAD_seqDist_DT.Rdata"))))

nRandomTADs <- length(randomTADsTSS)

randomTADs_dist_dna <- unlist(lapply(randomTADsTSS, function(x) x$mean_dist_dna))
randomTADs_dist_alignment <- unlist(lapply(randomTADsTSS, function(x) x$mean_dist_alignment))


# outFile2 <- file.path(outFold, paste0("pooledChromo_", type1, "Breakpoints_", type2, caller, "TADs_density.", plotType))
# do.call(plotType, list(outFile2, height = myHeight_dens, width=myWidth_dens))
plot_multiDens(list(
  observed = obs_mean_dist_dna,
  randomGenes =randomGenes_dist_dna,
  randomTADs = randomTADs_dist_dna
),
my_xlab = "sequence distance",
plotTit = paste0(mytit)
)
# foo <- dev.off()
# cat(paste0("... written: ", outFile2, "\n"))
