#!/usr/bin/bash

#Rscript coexpr_dist_v3.R TCGAcrc_msi_mss hgnc
Rscript coexpr_dist_v3.R GSE74927_neg_pos hgnc
Rscript coexpr_dist_v3.R GSE102073_stic_nostic hgnc   

Rscript coexpr_dist_v3.R GSE65540_before_after hgnc
Rscript coexpr_dist_v3.R GSE84231_lhb_rhb hgnc
Rscript coexpr_dist_v3.R GSE86356_tibMD1_tibNorm hgnc


