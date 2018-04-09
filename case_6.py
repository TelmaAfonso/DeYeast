
'''
Functions for case 6 simulations

Author: Telma Afonso
'''

from phenomenaly.simulation import fba, fva, pfba, lmoma
from phenomenaly.variables import Media
from types import *
import pickle
import string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.backends.backend_pdf import PdfPages

import warnings
warnings.filterwarnings('ignore')

from yeastpack_test import PhenomenalySim


class Case6 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case6, self).__int__(self, cobra_model)

    def loadExperimentalRes (self, filename_exp, sep = ';'):
        exp_fluxes = pd.read_csv(filename_exp, sep = sep)
        exp_fluxes.set_index('Yeast7_ID', inplace = True)
        reactions = exp_fluxes.iloc[:,0]
        sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if 'Exp Flux' not in col], 1)
        sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.')).astype('float')
        sub_exp_fluxes.columns = ['MAE1 Real Flux', 'WT Real Flux']

        #Add Biomass Values (r_4041)
        reactions['r_4041'] = 'Biomass'
        sub_exp_fluxes.loc['r_4041'] = [0.1, 0.1] #growth rate

        return sub_exp_fluxes, reactions

    def singleSimulation (self, gl_lb = -1.5, geneko = None, o2_lb = None, type = 'fba'):
        with self.model as m:
            m.set_carbon_source('r_1714', lb = gl_lb)
            if geneko is not None:
                m.set_environmental_conditions(gene_knockout = geneko)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb)

            if type == 'fba':
                res = fba(m)
            elif type == 'pfba':
                res = pfba(m)
            elif type == 'fva':
                res = fva(m, reaction_list = m.reactions, fix_biomass = True)
            elif type == 'lmoma':
                r = pfba(m)
                res = lmoma(m, r)

        return res

    def getDFWithoutExtremeFluxes (self, dataframe, column_index = 1, val = 900):
        df = dataframe[(dataframe.ix[:, column_index] > -val) & (dataframe.ix[:, column_index] < val)]

        return df.astype('float')

    def createDatasetExpVsSimul (self, exp_fluxes, sim_fluxes):
        dsim = sim_fluxes.copy()
        df = pd.concat([exp_fluxes, sim_fluxes], axis = 1, join = 'inner')

        return df

    def createDatasetWithAbsRelError (self, dataset_exp_vs_sim):
        df = dataset_exp_vs_sim.copy()
        ae = self.absoluteError(df.ix[:, 0], df.ix[:, 1])
        df.insert(loc = 2, column = 'Abs Error', value = ae)

        re = self.relativeError(df.ix[:, 0], df.ix[:, 1])
        df.insert(loc = 3, column = 'Rel Error', value = re)

        return df

    def simulationPipeline (self, exp_dataset, gl_lb = -1.5, geneko = None, o2_lb = None, type = 'fba', res_exists = False):
        if geneko is not None:
            path = 'Results/Case 6/res_{}_case6.sav'.format('mae1_' + type)
        else:
            path = 'Results/Case 6/res_{}_case6.sav'.format('wt_' + type)

        if res_exists:
            res = self.loadObjectFromFile(path)
        else:
            res = self.singleSimulation (gl_lb = gl_lb, geneko = geneko, o2_lb = o2_lb, type = type)
            self.saveObjectToFile(res, path)

        if type == 'lmoma':
            fluxes = pd.DataFrame(list(res.x_dict.items())).set_index(0) #legacy solution
        elif type == 'fva':
            fluxes = res #pandas df
        else:
            fluxes = res.fluxes

        # Dataset with experimental vs simulated fluxes
        df_exp_sim = self.createDatasetExpVsSimul(exp_dataset, self.correctReversedReactions(fluxes))
        if type != 'fva':
            df_exp_sim = df_exp_sim.rename(columns = {1: 'Sim Flux'})

        # Dataset with absolute and realtive errors
        df_exp_sim_errors = self.createDatasetWithAbsRelError(df_exp_sim)

        return res, df_exp_sim, df_exp_sim_errors


if __name__ == '__main__':

    #Initialization
    case6 = Case6()
    case6.model = case6.loadObjectFromFile('model_yeast_76.sav')
    case6.model.solver = 'optlang-cplex'
    case6.setMedium('MINIMAL')

    #General datasets
    exp_dataset, reactions = case6.loadExperimentalRes('Results/Case 6/case6_experimental_fluxes.csv')
    # NO EXPERIMENTAL ETHANOL FLUXES TO ADJUST O2 LB

    # ====== WILD TYPE ======
    #FBA
    wt_fba_res, wt_fba_exp_sim, wt_fba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'fba', res_exists = True)
    case6.plotExpVsSim(wt_fba_exp_sim_errors, save_fig_path = 'Results/Case 6/wt_fba_exp_sim_plot.png', title = 'FBA Wild Type')
    plt.close('all')

    #pFBA
    wt_pfba_res, wt_pfba_exp_sim, wt_pfba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'pfba', res_exists = True)
    case6.plotExpVsSim(wt_pfba_exp_sim_errors, save_fig_path = 'Results/Case 6/wt_pfba_exp_sim_plot.png', title = 'pFBA Wild Type')
    plt.close('all')

    # case6.getListOfMetabolitesSummary(wt_pfba_res)
    # case6.getMetaboliteSummaryWithNames('s_0373', wt_pfba_res)

    #FVA
    wt_fva_res, wt_fva_exp_sim, _ = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'fva', res_exists = True)


    # ====== MAE1 DELETION ======

    mae1 = case6.convertStdToSyst(['MAE1'])['MAE1']

    #FBA
    mae1_fba_res, mae1_fba_exp_sim, mae1_fba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'fba', res_exists = True)
    mae1_fba_exp_sim_errors = case6.getDFWithoutExtremeFluxes(mae1_fba_exp_sim_errors)
    case6.plotExpVsSim(mae1_fba_exp_sim_errors, save_fig_path = 'Results/Case 6/mae1_fba_exp_sim_plot.png', title = 'FBA MAE1 Del')
    plt.close('all')

    #pFBA
    mae1_pfba_res, mae1_pfba_exp_sim, mae1_pfba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'pfba', res_exists = True)
    case6.plotExpVsSim(mae1_pfba_exp_sim_errors, save_fig_path = 'Results/Case 6/mae1_pfba_exp_sim_plot.png', title = 'pFBA MAE1 Del')
    plt.close('all')

    #LMOMA
    mae1_lmoma_res, mae1_lmoma_exp_sim, mae1_lmoma_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'lmoma', res_exists = True)
    case6.plotExpVsSim(mae1_lmoma_exp_sim_errors, save_fig_path = 'Results/Case 6/mae1_lmoma_exp_sim_plot.png', title = 'LMOMA MAE1 Del')
    plt.close('all')

    # case6.getListOfMetabolitesSummary(mae1_pfba_res)
    # case6.getMetaboliteSummaryWithNames('s_0373', wt_pfba_res)

    #FVA
    mae1_fva_res, mae1_fva_exp_sim, _ = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'fva', res_exists = True)



    # case6.checkReaction('r_0302')
    # case6.checkReactionLB('r_1992')


    # =========================================
    #    Save all results into a binary file
    # =========================================

    all_res = {'d6_wt_fba': wt_fba_exp_sim, 'd6_wt_pfba': wt_pfba_exp_sim, 'd6_wt_fva': wt_fva_exp_sim, 'd6_reactions': reactions,
               'd6_mae1_fba': mae1_fba_exp_sim, 'd6_mae1_pfba': mae1_pfba_exp_sim, 'd6_mae1_fva': mae1_fva_exp_sim, 'd6_mae1_lmoma': mae1_lmoma_exp_sim}

    case6.saveObjectToFile(all_res, 'Results/case6_all_res.sav')









