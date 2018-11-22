curr_TAD <- "chr1_TAD150"

test_fams <- c("unassigned", "RNA binding motif containing")

tad_g2t <- gene2tadDT[gene2tadDT$region == curr_TAD & gene2tadDT$entrezID %in% pipeline_geneList,]


tad_g2t$family <- unlist(sapply(tad_g2t$entrezID, function(x) 
  if(x %in% familyDT$entrezID) { 
    return(familyDT$hgnc_family[familyDT$entrezID == x])
  } else {
    return("unassigned")  
  }
  ))

tad_g2t

for(i_fam in test_fams) {
  cat(i_fam, ":")
  cat(sum(tad_g2t$family == i_fam))
  cat("\n")
}
cat("tot n genes: ")
nrow(tad_g2t)
cat("\n")






mod1Members <- ( as.character(gene2mod1_topTADs) == mod1Labels[i_mod1] )
mod2Members <- ( as.character(gene2mod2_topTADs)  == mod2Labels[i_mod2])
pval <- -log10(fisher.test(mod1Members, mod2Members, alternative = "greater")$p.value);

mod1 <- c(
  geneA = "chr1_TAD1",
  geneB = "chr1_TAD1", 
  geneC = "chr1_TAD2",
  geneD = "chr1_TAD3",
  geneE = "chr1_TAD4",
  geneF = "chr1_TAD1"
)

mod2 <- c(
  geneA = "fam1",
  geneB = "fam1", 
  geneC = "fam1",
  geneD = "fam2",
  geneE = "fam3",
  geneF = "fam4"
)

i_mod1 = "chr1_TAD1"
i_mod2 = "fam1"

mod1Members <- ( as.character(mod1) == i_mod1 )
mod2Members <- ( as.character(mod2) == i_mod2 )

fisher.test(mod1Members, mod2Members)

testMat <- matrix(c(2,0,1,3), nrow=2,
                  byrow = T, dimnames = list(c("sameTAD","diffTAD"), c("sameFam","diffFam")))
fisher.test(testMat)


testMat <- matrix(c(2,1,1,2), nrow=2,
                  byrow = T, dimnames = list(c("sameTAD","diffTAD"), c("sameFam","diffFam")))
fisher.test(testMat)


