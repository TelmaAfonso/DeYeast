
'''
Case 7 simulations

Author: Telma Afonso
'''

import warnings
from phenomenaly.io import load_yeast_76
from case_7 import *

warnings.filterwarnings('ignore')

#Initialization
# case7 = Case7(cobra_model = load_yeast_76()) #Load model using phenomenally
# case7.saveObjectToFile(case7.model, 'model_yeast_76.sav')
case7 = Case7()
case7.model = case7.loadObjectFromFile('model_yeast_76.sav')
case7.model.solver = 'optlang-cplex'
case7.setMedium('MINIMAL_CASE7')
case7.dictsForCase7()

# General datasets
exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2 = case7.case7Pipeline(plot = False, makeFigs = False, type = 'pfba', res_exists = True)

# FBA
res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors = case7.fbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, plotGenes = False, plotReacts = False, saveGenesPlot = False, res_exists = True)

# pFBA
res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors = case7.pfbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, saveGenesPlot = False, plotReacts = False, plotGenes = False, res_exists = True)

# FVA
res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim = case7.fvaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, res_exists = True)

# LMOMA
res_lmoma, res_lmoma_df, wt_lmoma_df, df_lmoma_exp_sim, df_lmoma_exp_sim_errors = case7.lmomaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, reference_dict = res_pfba, plotGenes = False, plotReacts = False, saveGenesPlot = False, res_exists = True)

# Create xlsx with results
# case7.convertPandasDFToExcel(reactions, df_fba_exp_sim_errors, filename = 'Results/Case 7/fba_results_case7_test.xlsx', imgFolder = 'Results/Case 7/FBA_figs')
# case7.convertPandasDFToExcel(reactions, df_pfba_exp_sim_errors, title = 'pFBA Results Case 7', filename = 'Results/Case 7/pfba_results_case7_test.xlsx', imgFolder = 'Results/Case 7/pFBA_figs', type = 'pfba')
# case7.convertPandasDFToExcelFVA(reactions, df_fva_exp_sim, title = 'FVA Results Case 7', filename = 'Results/Case 7/fva_results_case7.xlsx')
# case7.convertPandasDFToExcel(reactions, df_lmoma_exp_sim_errors, title = 'LMOMA Results Case 7', filename = 'Results/Case 7/lmoma_results_case7.xlsx', imgFolder = 'Results/Case 7/LMOMA_figs', type = 'lmoma')


#TESTS
genes_res = case7.createResultsDictByGene(df_fba_exp_sim_errors, df_pfba_exp_sim_errors, df_lmoma_exp_sim_errors, df_fva_exp_sim)
wt_res = case7.createResultsDataframeWT(reactions, wt_fba_df, wt_pfba_df, wt_lmoma_df, wt_fva_df)


#FIXING BIOMASS UB

# General datasets

# # FBA
# res_fba_biomub, res_fba_df_biomub, wt_fba_df_biomub, df_fba_exp_sim_biomub, df_fba_exp_sim_errors_biomub = case7.fbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, plotGenes = False, plotReacts = False, saveGenesPlot = False, res_exists = False, fname = 'res_fba_biomub_case7.sav', biom_ub = True)
#
# case7.plotReactExpVsSim(df_pfba_exp_sim_biomub, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'r_2116', reaction = 'r_2116')
#
# # pFBA
# res_pfba_biomub, res_pfba_df_biomub, wt_pfba_df_biomub, df_pfba_exp_sim_biomub, df_pfba_exp_sim_errors_biomub = case7.pfbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, saveGenesPlot = False, plotReacts = False, plotGenes = False, res_exists = False, fname = 'res_pfba_biomub_case7.sav', biom_ub = True)
#
# # FVA
# res_fva_biomub, res_fva_df_biomub, wt_fva_df_biomub, df_fva_exp_sim_biomub = case7.fvaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, res_exists = False, fname = 'res_fva_biomub_case7.sav', biom_ub = True)
#
# # LMOMA
# res_lmoma_biomub, res_lmoma_df_biomub, wt_lmoma_df_biomub, df_lmoma_exp_sim_biomub, df_lmoma_exp_sim_errors_biomub = case7.lmomaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, reference_dict = res_pfba, plotGenes = False, plotReacts = False, saveGenesPlot = False, res_exists = False, fname = 'res_lmoma_biomub_case7.sav', biom_ub = True)

res = case7.singleSimulation('ADH3', fluxes_O2 = fluxes_O2, gl_lb = None, o2_lb = None, type = 'fba')
res.fluxes.loc['r_4041']
res.fluxes.loc['r_2116']




