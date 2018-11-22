#=> check that all genes are present in gene1+gene2

# coexpr

cat("start coexpr\n")

v1 <- eval(parse(text = load("//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_COEXPR_SORTNODUP/TCGAcrc_msi_mss_pearson/coexprDT.Rdata")))
v0 <- eval(parse(text = load("//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_COEXPR_DONT_USE_ME_NOT_SORTED/TCGAcrc_msi_mss_pearson/coexprDT.Rdata")))

v1$gene1 <- as.character(v1$gene1)
v1$gene2 <- as.character(v1$gene2)

v0$gene1 <- as.character(v0$gene1)
v0$gene2 <- as.character(v0$gene2)

v1_set <- union(v1$gene1, v1$gene2)
v2_set <- union(v0$gene1, v0$gene2)

stopifnot(setequal(v1_set, v2_set))

all(v1$gene1 == v0$gene1)

# fam
cat("start fam\n")
v1 <- eval(parse(text = load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_SAME_FAMILY_SORTNODUP/hgnc_family_all_family_pairs.Rdata")))
v0 <- eval(parse(text = load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_SAME_FAMILY_DONT_USE_ME_NOT_SORTED/hgnc_family_all_family_pairs.Rdata")))


v1$gene1 <- as.character(v1$gene1)
v1$gene2 <- as.character(v1$gene2)

v0$gene1 <- as.character(v0$gene1)
v0$gene2 <- as.character(v0$gene2)

v1_set <- union(v1$gene1, v1$gene2)
v2_set <- union(v0$gene1, v0$gene2)

stopifnot(setequal(v1_set, v2_set))

all(v1$gene1 == v0$gene1)

cat("start dist\n")
# dist


v1 <- eval(parse(text = load("//mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_DIST_SORTNODUP/all_dist_pairs.Rdata")))
v0 <- eval(parse(text = load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_DIST_DONT_USE_ME_NOT_SORTED/all_dist_pairs.Rdata")))


v1$gene1 <- as.character(v1$gene1)
v1$gene2 <- as.character(v1$gene2)

v0$gene1 <- as.character(v0$gene1)
v0$gene2 <- as.character(v0$gene2)

v1_set <- union(v1$gene1, v1$gene2)
v2_set <- union(v0$gene1, v0$gene2)

stopifnot(setequal(v1_set, v2_set))

# tad
cat("start tad\n")

v1 <- eval(parse(text = load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_SAME_TAD_SORTNODUP/all_TAD_pairs.Rdata")))
v0 <- eval(parse(text = load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_SAME_TAD_DONT_USE_ME_NOT_SORTED/all_TAD_pairs.Rdata")))

v1$gene1 <- as.character(v1$gene1)
v1$gene2 <- as.character(v1$gene2)

v0$gene1 <- as.character(v0$gene1)
v0$gene2 <- as.character(v0$gene2)

v1_set <- union(v1$gene1, v1$gene2)
v2_set <- union(v0$gene1, v0$gene2)

all(v1$gene1 == v0$gene1)

stopifnot(setequal(v1_set, v2_set))

