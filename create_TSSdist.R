suppressPackageStartupMessages(library(Biostrings, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(DECIPHER, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(fossil, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2,20))

caller <- "TopDom"

# Rscript create_TSSdist.R

windowBeforeTSS <- 1500
windowAfterTSS <- 500

nCpu <- ifelse(SSHFS, 2, 30)

genePairsDT <- eval(parse(
  text = load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_COMBN/all_gene_pairs_DT.Rdata"))))

# genePairsDT <- genePairsDT[1:5,]

stopifnot(colnames(genePairsDT) == c("gene1", "gene2"))
genePairsDT$gene1 <- as.character(genePairsDT$gene1)
genePairsDT$gene2 <- as.character(genePairsDT$gene2)

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_TSS_DIST") )
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "create_TSS_DIST.txt")
system(paste0("rm -f ", logFile))

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
entrezDT$entrezID <- as.character(entrezDT$entrezID)


geneList_DT <- entrezDT
geneList_DT$gene_start <- ifelse(geneList_DT$strand == "+", geneList_DT$start, geneList_DT$end)


#****************************************************************************************************
#**************************************************************************************************** sequence alignment for each TAD
#****************************************************************************************************

all_pair_seqDist_DT <- foreach( i = 1:nrow(genePairsDT), .combine = 'rbind') %dopar% {
  
  cat(paste0("... gene pair ", i, "/", nrow(genePairsDT), "\n"))
  
  gene1 <- genePairsDT$gene1[i]
  gene2 <- genePairsDT$gene2[i]
  
  curr_genes <- c(gene1, gene2)
  
  pair_geneList_DT <- geneList_DT[geneList_DT$entrezID %in% curr_genes,]
  
  if(nrow(pair_geneList_DT) == 1) {
    return(
      data.frame(
        gene1 = gene1,
        gene2 = gene2,
        mean_dist_dna = NA,
        mean_dist_alignment = NA,
        # cor_dna_alignment = NA,
        stringsAsFactors = FALSE
      )
    )
    
  }
  
  pair_seqStringSet <- getSeq(x = Hsapiens, 
                             names = pair_geneList_DT$chromo,
                             start = pair_geneList_DT$gene_start - windowBeforeTSS, 
                             end = pair_geneList_DT$gene_start + windowAfterTSS,
                             strand = pair_geneList_DT$strand)
  
  # align the sequences
  pair_alignedSeq <- AlignSeqs(pair_seqStringSet)
  # processors = ifelse(SSHFS, 2, nCpu)) # run dopar!
  names(pair_alignedSeq) <- pair_geneList_DT$entrezID
  
  # distance between the aligned sequences: ape::dist.dna, applied on as.DNAbin
  pair_DNAbin <- as.DNAbin(pair_alignedSeq)
  pair_dist_dna <- dist.dna(pair_DNAbin, model="raw")
  
  # distance: seqinr::dist.alignment, applied on alignment object
  pair_alignment <- seqinr::as.alignment(
    nb = length(pair_alignedSeq),
    nam = names(pair_alignedSeq),
    seq = as.character(pair_alignedSeq)
  )
  pair_dist_alignment <- dist.alignment(pair_alignment)
  
  data.frame(
    gene1 = gene1,
    gene2 = gene2,
    mean_dist_dna = mean(pair_dist_dna, na.rm=T),
    mean_dist_alignment = mean(pair_dist_alignment, na.rm=T),
    # cor_dna_alignment = cor(pair_dist_alignment, pair_dist_dna),
    stringsAsFactors = FALSE
  )
}

# names(all_pair_seqDist_DT) <- regionList # THIS IS A DT !!!
outFile <- file.path(outFold, "all_pair_seqDist_DT.Rdata")
save(all_pair_seqDist_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
cat("*** DONE\n")
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)