#!/usr/bin/bash

# ./all_create_coexpr_sortNoDup.sh

# Rscript create_coexpr_sortNoDup.R GSE81046_noninf_list
# Rscript create_coexpr_sortNoDup.R GSE81046_noninf_salm
# Rscript create_coexpr_sortNoDup.R GSE81046_salm_list
# Rscript create_coexpr_sortNoDup.R TCGAprad_norm_prad


start_time=$(date -R)   

script_name="create_coexpr_sortNoDup.R"
 
#############################################################

all_data=( "TCGAprad_norm_prad" )

### TCGA datasets 08.12.2018
##all_data=(
#"TCGAacc_wt_mutCTNNB1"
#"TCGAblca_norm_blca"
##"TCGAbrca_lum_bas lum bas"
#"TCGAcesc_adeno_squam"
##"TCGAcoad_msi_mss msi mss"
#"TCGAgbm_classical_neural"
#"TCGAgbm_classical_proneural"
#"TCGAgbm_classical_mesenchymal"
#"TCGAhnsc_HPVneg_HPVpos"
#"TCGAkich_norm_kich"
#"TCGAlaml_wt_mutFLT3"
#"TCGAlgg_IDHwt_IDHmutnc"
#"TCGAlihc_wt_mutCTNNB1"
#"TCGAluad_mutKRAS_mutEGFR"
#"TCGAluad_nonsmoker_smoker nonsmoker smoker"
#"TCGAluad_wt_mutKRAS"
#"TCGAlusc_norm_lusc"
#"TCGApaad_wt_mutKRAS"
#"TCGAsarc_ddlps_lms"
#"TCGAsarc_ddlps_mfs"
#"TCGAsarc_lms_mfs"
#"TCGAskcm_lowInf_highInf"
#"TCGAskcm_wt_mutBRAF"
#"TCGAskcm_wt_mutCTNNB1"
#"TCGAthca_wt_mutBRAF"
#"TCGAstad_EBVpos_gs"
#"TCGAstad_msi_gs"
#"TCGAstad_norm_gs"
#"TCGAtgct_sem_nonsem"
#"TCGAthca_mut.RAS_mutBRAF"
#"TCGAucec_msi_cnl"
#)

#all_data=(
            #"TCGAacc_wt_mutCTNNB1"
            #"TCGAblca_norm_blca"
            ##"TCGAbrca_lum_bas lum bas"
            #"TCGAcesc_adeno_squam"
            ##"TCGAcoad_msi_mss msi mss"
            #"TCGAgbm_classical_neural"
                #"TCGAgbm_classical_proneural"
                #"TCGAgbm_classical_mesenchymal"
                #"TCGAhnsc_HPVneg_HPVpos"
                #"TCGAkich_norm_kich"
                #"TCGAlaml_wt_mutFLT3"
                #"TCGAlgg_IDHwt_IDHmutnc"
                #"TCGAlihc_wt_mutCTNNB1"
                #"TCGAluad_mutKRAS_mutEGFR"
#                        "TCGAluad_nonsmoker_smoker"
                        #"TCGAluad_wt_mutKRAS"
                        #"TCGAlusc_norm_lusc"
                        #"TCGApaad_wt_mutKRAS"
                        #"TCGAsarc_ddlps_lms"
                        #"TCGAsarc_ddlps_mfs"
                        #"TCGAsarc_lms_mfs"
                            #"TCGAskcm_lowInf_highInf"
                            #"TCGAskcm_wt_mutBRAF"
                            #"TCGAskcm_wt_mutCTNNB1"
                            #"TCGAthca_wt_mutBRAF"
                            #"TCGAstad_EBVpos_gs"
                            #"TCGAstad_msi_gs"
#"TCGAstad_norm_gs"
#"TCGAtgct_sem_nonsem"
#"TCGAthca_mut.RAS_mutBRAF"
#"TCGAucec_msi_cnl"
#)



## new datasets
#all_data=(
##"TCGAcesc_adeno_squam" 
#"TCGAhnsc_HPVneg_HPVpos" 
#"TCGAlgg_IDHwt_IDHmutnc" 
#"TCGAsarc_ddlps_lms" 
#"TCGAsarc_ddlps_mfs" 
#"TCGAsarc_lms_mfs" 
#"TCGAtgct_sem_nonsem" 
#)
## new datasets
#all_data=(
#"TCGAblca_norm_blca" 
#"TCGAkich_norm_kich" 
#"TCGAlusc_norm_lusc" 
#"TCGAstad_norm_gs" 
#)
## DESeq2
#all_data=(
#"GSE77509_normal_tumor" 
#"GSE77509_normal_ptt" 
#"GSE77509_ptt_tumor" 
#"GSE68719_norm_park" 
#"GSE64810_control_carrier" 
#"GSE101521_control_mdd" 
#"GSE77314_normal_tumor" 
#)

#all_data=(
#"GSE101521_control_mdd"
#"GSE102073_stic_nostic" # -> ok
#"GSE40419_normal_cancer"
#"GSE48166_control_ICM"
#"GSE51799_control_carrier"
#"GSE52166_prePf_postPf"
#"GSE57148_normal_COPD"
#"GSE58135_ERpos_adjERpos"
#"GSE58135_ERpos_tripleNeg"
#"GSE58135_tripleNeg_adjTripleNeg"
#"GSE61476_unaffectNPCday11_unaffectNPCday31"
#"GSE64810_control_carrier"
#"GSE64813_prePTSD_postPTSD"
#"GSE65540_before_after" # -> ok
#"GSE66306_before_after"
#"GSE67528_ASDiPSC_ASDneuron"
#"GSE67528_ASDiPSC_ASDnpc"
#"GSE67528_CTRLiPSC_CTRLneuron"
#"GSE67528_CTRLiPSC_CTRLnpc"
#"GSE68719_norm_park"
#"GSE71119_dediffSM_MFSM"
#"GSE71119_undiffSM_LMSM"
#"GSE73765_noninf_list"
#"GSE73765_noninf_salm"
#"GSE73765_salm_list"
# "GSE74927_neg_pos" # -> ok
#"GSE77314_normal_tumor"
#"GSE77509_normal_ptt"
#"GSE77509_normal_tumor"
#"GSE77509_ptt_tumor"
#"GSE78936_ba11Ctrl_ba11BD"
#"GSE78936_ba11Ctrl_ba11Sz"
#"GSE79209_dysp_nodysp"
#"GSE79362_prog_nonprog"
#"GSE81046_noninf_list"
#"GSE81046_noninf_salm"
#"GSE81046_salm_list"
#"GSE81089_normal_nsclc"
#"GSE84231_lhb_lsc"
#"GSE84231_lhb_rhb" # -> ok
#"GSE84231_lsc_rsc"
#"GSE84231_rhb_rsc"
#"GSE86356_quadMD1_quadNorm"
#"GSE86356_tibMD1_quadMD1"
#"GSE86356_tibMD1_tibNorm" # -> ok
#"GSE86422_myoPre_myoPost"
#"GSE87194_control_schi"
#"GSE87340_ad_nl"
#"GSE90749_iPSC_adipo"
#"GSE90749_iPSC_hepato"
#"GSE92592_control_ipf"
#"GSE94631_adMono_adDC"
#"GSE94736_old_young"
#"TCGAacc_acc_mutCTNNB1"
#"TCGAbrca_lum_bas"
#"TCGAcrc_msi_mss" # -> ok
#"TCGAlaml_laml_mutFLT3"
#"TCGAlihc_lihc_mutCTNNB1"
#"TCGAluad_luad_mutKRAS"
#"TCGApaad_paad_mutKRAS"
#"TCGAskcm_skcm_mutBRAF"
#"TCGAskcm_skcm_mutCTNNB1"
#"TCGAstad_msi_gs"
#"TCGAthca_thca_mutBRAF"
#"TCGAucec_msi_cnl"
#)
#############################################################
for data in "${all_data[@]}"; do
    echo "> START for $data"
	echo Rscript $script_name $data
	Rscript $script_name $data
done




###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

