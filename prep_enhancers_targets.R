suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

genehancer_DT <- read.delim("genehancer.csv", header=T, stringsAsFactors = F, sep=",")
nrow(genehancer_DT)
genehancer_DT <- genehancer_DT[grepl("enhancer", tolower(genehancer_DT$feature.name)),]
nrow(genehancer_DT)

head(genehancer_DT)

source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/main_settings.R"))
source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R"))

outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/GeneHancer_GENES_MAPPED")
system(paste0("mkdir -p ", outFold))

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)


# prep_enhancer_DT <- foreach(i = 1:1000, .combine="rbind") %do% {
  prep_enhancer_DT <- foreach(i = 1:nrow(genehancer_DT), .combine="rbind") %dopar% {
  
  cat(paste0("... start enhancer: ", i, "/", nrow(genehancer_DT), "\n"))
  
  curr_att <- genehancer_DT$attributes[i]
  all_att <- strsplit(x=curr_att,split=";")[[1]]
  curr_enhancer <- all_att[grep("genehancer_id", all_att)]
  stopifnot(length(curr_enhancer) == 1)
  curr_enhancer <- gsub("genehancer_id=", "", curr_enhancer)
  
  mapped_genes <- all_att[grep("connected_gene", all_att)]
  mapped_genes <- gsub("connected_gene=", "", mapped_genes)
  
  mapped_genes_symbol <- get_geneList_fromSymbol(refList=mapped_genes, g2t=gene2tadDT, symbDT_file=symbolDT_file)
  mapped_genes_entrez <- get_geneList_fromEntrez(refList=mapped_genes, g2t=gene2tadDT, histDT_file=historyDT_file)
  mapped_genes_ensembl <- get_geneList_fromEnsembl(refList=mapped_genes, g2t=gene2tadDT, ensDT_file=ensemblDT_file)  
  mapped_genes_renamed <- c(mapped_genes_symbol,mapped_genes_entrez, mapped_genes_ensembl)
  
  if(length(mapped_genes_renamed) == 0) {
    data.frame(
      chromo = genehancer_DT$chrom[i],
      start = genehancer_DT$start[i],
      end = genehancer_DT$end[i],
      enhancer = curr_enhancer,
      connected_gene = mapped_genes,
      connected_gene_renamed = NA,
      connected_gene_chromo = NA,
      connected_gene_start = NA,
      connected_gene_end = NA,
      stringsAsFactors = F
    )
  }else {
    stopifnot(all(mapped_genes_renamed %in% gene2tadDT$entrezID))
    mapped_genes_renamed_chromo <- unlist(sapply(mapped_genes_renamed, function(x) gene2tadDT$chromo[gene2tadDT$entrezID == x]))
    mapped_genes_renamed_start <- unlist(sapply(mapped_genes_renamed, function(x) gene2tadDT$start[gene2tadDT$entrezID == x]))
    mapped_genes_renamed_end <- unlist(sapply(mapped_genes_renamed, function(x) gene2tadDT$end[gene2tadDT$entrezID == x]))
    # chromo = genehancer_DT$chrom[i]
    # start = genehancer_DT$start[i]
    # end = genehancer_DT$end[i]
    # enhancer = curr_enhancer
    # connected_gene = mapped_genes
    # connected_gene_renamed = as.character(mapped_genes_renamed[mapped_genes])
    # connected_gene_chromo = as.character(mapped_genes_renamed_chromo[mapped_genes])
    # connected_gene_start = as.numeric(mapped_genes_renamed_start[mapped_genes])
    # connected_gene_end = as.numeric(mapped_genes_renamed_end[mapped_genes])
    data.frame(
      chromo = genehancer_DT$chrom[i],
      start = genehancer_DT$start[i],
      end = genehancer_DT$end[i],
      enhancer = curr_enhancer,
      connected_gene = mapped_genes,
      connected_gene_renamed = as.character(mapped_genes_renamed[mapped_genes]),
      connected_gene_chromo = as.character(mapped_genes_renamed_chromo[mapped_genes]),
      connected_gene_start = as.numeric(mapped_genes_renamed_start[mapped_genes]),
      connected_gene_end = as.numeric(mapped_genes_renamed_end[mapped_genes]),
      stringsAsFactors = F
    )    
  }

  
}
outFile <- file.path(outFold, "prep_enhancer_DT.Rdata")
save(prep_enhancer_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
