##################################################################################################################################
# Runs scripts that run analyses in Trachsel et al. 2020 Quaternary research 
##################################################################################################################################

#define working directory
wd <- '~/r_code_calibration_paper/'


source(paste0(wd,'expert_elicitation/R/Elicitation Plots.R'))  #works


# load the data that is not in neotoma
source(paste0(wd,'calibration/R/source_data_not_in_neotoma.R'))
source(paste0(wd,'calibration/helper_funs/load_vegetation_data_mean.R'))
source(paste0(wd,'calibration/R/prep_data_calibration_only_abies_changed_other_hardwood_final.R'))
#from here calibration is run in a bash script or on condor




#----------------------------------------------------------------------------------------------------------------------------------
#prepare data for spatial process
#----------------------------------------------------------------------------------------------------------------------------------
source(paste0(wd,'vegetation/R/prepare_veg_data_township.R'))
#----------------------------------------------------------------------------------------------------------------------------------
# prepare data for prediction
#----------------------------------------------------------------------------------------------------------------------------------
source(paste0(wd,'calibration/R/prep_data_for_prediction.R'))
source(paste0(wd,'prediction/R/prepare_data_settlement_era_final_var110_december.R'))
#----------------------------------------------------------------------------------------------------------------------------------
#analyze predicted data
#----------------------------------------------------------------------------------------------------------------------------------
source(paste0(wd,'prediction/R/analysis _of_posterior_nb_13_taxa_combine_runs.R'))
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
#run REVEALS
#----------------------------------------------------------------------------------------------------------------------------------
source(paste0(wd,'Reveals_NEUS/R/eval_ppes.r'))
source(paste0(wd,'Reveals_NEUS/R/reveals.r'))
source(paste0(wd,'Reveals_NEUS/R/gridded_REVEALS.R'))
source(paste0(wd,'prediction/R/Chord_distance_observed_predicted_vegetation.R'))
#----------------------------------------------------------------------------------------------------------------------------------
#analyses/Figures shown in paper
#----------------------------------------------------------------------------------------------------------------------------------
source(paste0(wd,'figures/R/Figure1.R'))
source(paste0(wd,'figures/R/inset_map_figure1.R'))
source(paste0(wd,'figures/R/Figure_2.r'))
source(paste0(wd,'figures/R/Figure5.R'))
source(paste0(wd,'figures/R/Figure_6.R'))



source(paste0(wd,'figures/R/Figure7.R'))




