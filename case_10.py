
'''
Functions for case 10 simulations

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
pd.set_option('display.max_colwidth', -1)

from yeastpack_test import PhenomenalySim


class Case10 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case10, self).__int__(self, cobra_model)

    def dictsForCase10 (self):
        # Carbon source lbs
        #cs_lb = {'glucose': -1.5, 'maltose': -1, 'ethanol': -5, 'acetate': -5} #Shophia's values
        cs_lb = {'glucose': -1.15, 'maltose': -0.61, 'ethanol': -3.78, 'acetate': -5.89} #Author's values

        # Oxygen lbs
        o2_lb = {'glucose': -2.74, 'maltose': -3.05, 'ethanol': -6.87, 'acetate': -7.4}

        # Carbon source exchange reactions
        cs_reaction = {'glucose': 'r_1714', 'maltose': 'r_1931', 'ethanol': 'r_1761', 'acetate': 'r_1634'}

        self.cs_lb, self.o2_lb, self.cs_reaction = cs_lb, o2_lb, cs_reaction

    def loadExperimentalRes (self, filename_exp, sep = ';'):
        exp_fluxes = pd.read_csv(filename_exp, sep = sep)
        exp_fluxes.set_index('Yeast7_ID', inplace = True)
        reactions = exp_fluxes.iloc[:,0]
        sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if 'exp flux' not in col], 1)
        sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.')).astype('float')
        #sub_exp_fluxes.columns = ['MAE1 Real Flux', 'WT Real Flux']

        #Add Biomass Values (r_4041)
        # reactions['r_4041'] = 'Biomass'
        # sub_exp_fluxes.loc['r_4041'] = [0.1, 0.1] #growth rate

        return sub_exp_fluxes, reactions

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

    def simulationPipeline (self, exp_dataset, cs = 'glucose', geneko = None, type = 'fba', res_exists = False, fname = None):
        if res_exists:
            res = self.loadObjectFromFile(fname)
        else:
            res = self.singleSimulation(carbon_source = self.cs_reaction[cs], cs_lb = self.cs_lb[cs], geneko = geneko, o2_lb = self.o2_lb[cs], type = type)
            self.saveObjectToFile(res, fname)

        if type == 'fva':
            fluxes = res #pandas df
        elif hasattr(res, 'x_dict'): #if legacy solution
            fluxes = pd.DataFrame(list(res.x_dict.items())).set_index(0)
        else:
            fluxes = res.fluxes

        # Dataset with experimental vs simulated fluxes
        df_exp_sim = self.createDatasetExpVsSimul(exp_dataset, self.correctReversedReactions(fluxes))
        if type != 'fva':
            df_exp_sim = df_exp_sim.rename(columns = {1: 'Sim Flux'})

        # Dataset with absolute and realtive errors
        df_exp_sim_errors = self.createDatasetWithAbsRelError(df_exp_sim)

        return res, df_exp_sim, df_exp_sim_errors

    def plotExpVsSim (self, absRelErrorDataset, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'Wild Type', label_adjust = 0.05, save_fig_path = None):
        plt.rcParams["figure.figsize"] = (10,5)

        x = absRelErrorDataset.ix[:,0]
        y = absRelErrorDataset.ix[:,1]
        react_IDs = list(absRelErrorDataset.index)
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        meanRelErr = absRelErrorDataset.ix[:,3].mean()
        corr = x.corr(y)

        plt.plot(x, y, 'o', x, line)
        for ind, react_ID in enumerate(react_IDs):
            plt.annotate(react_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + label_adjust, y[ind] + label_adjust))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.plot([], [], ' ') # To show correlation in legend
        plt.plot([], [], ' ') # To show mean relative error in legend
        plt.legend(['Reactions', 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

        if save_fig_path is not None:
            plt.savefig(save_fig_path)


if __name__ == '__main__':

    #Initialization
    case10 = Case10()
    case10.model = case10.loadObjectFromFile('model_yeast_76.sav')
    case10.model.solver = 'optlang-cplex'
    case10.setMedium('MINIMAL')
    case10.dictsForCase10()

    #General datasets
    exp_dataset, reactions = case10.loadExperimentalRes('Results/Case 10/case10_experimental_fluxes.csv')

    # ====== CS: GLUCOSE ======
    #FBA
    g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,0], cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_glucose_case10.sav')
    g_fba_exp_sim_errors = case10.getDFWithoutExtremeFluxes(g_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case10.plotExpVsSim(g_fba_exp_sim_errors, save_fig_path = 'Results/Case 10/g_fba_exp_sim_plot.png', title = 'FBA GLucose Carbon Source')
    plt.close('all')

    #pFBA
    g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,0], cs = 'glucose', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_glucose_case10.sav')
    case10.plotExpVsSim(g_pfba_exp_sim_errors, save_fig_path = 'Results/Case 10/g_pfba_exp_sim_plot.png', title = 'pFBA GLucose Carbon Source')
    plt.close('all')

    #FVA
    g_fva_res, g_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,0], cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_glucose_case10.sav')


    # ====== CS: MALTOSE ======
    # Maltose only participates in two reactions (transport and exchange)
    # Consider maltose as glucose?? (maltose <==> D-Glucose, relative flux of 100)
    #FBA
    m_fba_res, m_fba_exp_sim, m_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,1], cs = 'maltose', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_maltose_case10.sav')
    m_fba_exp_sim_errors = case10.getDFWithoutExtremeFluxes(m_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case10.plotExpVsSim(m_fba_exp_sim_errors, save_fig_path = 'Results/Case 10/m_fba_exp_sim_plot.png', title = 'FBA Maltose Carbon Source')
    plt.close('all')

    #pFBA
    m_pfba_res, m_pfba_exp_sim, m_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,1], cs = 'maltose', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_maltose_case10.sav')
    case10.plotExpVsSim(m_pfba_exp_sim_errors, save_fig_path = 'Results/Case 10/m_pfba_exp_sim_plot.png', title = 'pFBA Maltose Carbon Source')
    plt.close('all')

    #FVA
    m_fva_res, m_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,1], cs = 'maltose', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_maltose_case10.sav')


    # ====== CS: ETHANOL ======
    #FBA
    e_fba_res, e_fba_exp_sim, e_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,2], cs = 'ethanol', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_ethanol_case10.sav')
    e_fba_exp_sim_errors = case10.getDFWithoutExtremeFluxes(e_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case10.plotExpVsSim(e_fba_exp_sim_errors, save_fig_path = 'Results/Case 10/e_fba_exp_sim_plot.png', title = 'FBA Ethanol Carbon Source')
    plt.close('all')

    #pFBA
    e_pfba_res, e_pfba_exp_sim, e_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,2], cs = 'ethanol', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_ethanol_case10.sav')
    case10.plotExpVsSim(e_pfba_exp_sim_errors, save_fig_path = 'Results/Case 10/e_pfba_exp_sim_plot.png', title = 'pFBA Ethanol Carbon Source')
    plt.close('all')

    #FVA
    e_fva_res, e_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,2], cs = 'ethanol', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_ethanol_case10.sav')


    # ====== CS: ACETATE ======
    #FBA
    a_fba_res, a_fba_exp_sim, a_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,3], cs = 'acetate', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_acetate_case10.sav')
    a_fba_exp_sim_errors = case10.getDFWithoutExtremeFluxes(a_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case10.plotExpVsSim(a_fba_exp_sim_errors, save_fig_path = 'Results/Case 10/a_fba_exp_sim_plot.png', title = 'FBA Acetate Carbon Source')
    plt.close('all')

    #pFBA
    a_pfba_res, a_pfba_exp_sim, a_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,3], cs = 'acetate', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_acetate_case10.sav')
    case10.plotExpVsSim(a_pfba_exp_sim_errors, save_fig_path = 'Results/Case 10/a_pfba_exp_sim_plot.png', title = 'pFBA Acetate Carbon Source')
    plt.close('all')

    #FVA
    a_fva_res, a_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,3], cs = 'acetate', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_acetate_case10.sav')



    # TESTS
    # r_0717, r_0982 inverse
    #
    # case10.checkReaction(case6.convertKeggID('R02035'))
    # case10.getReactionInfo('r_0982')
    #
    # case10.checkReaction('r_0567')
    # r = pd.DataFrame(reactions)
    # r.index
    #
    # for i in r.index:
    #     print(i, str(i) + '_reversed' in res.keys())
    #
    # a_fba_exp_sim_errors.loc['r_1049']
    # res = a_pfba_res.x_dict
    # res.keys()
    #
    # a = case10.model.reactions.get_by_id('r_1049')
    # vars(a)
    #
    # reactions = ['r_0962', 'r_0300', 'r_1022', 'r_0454', 'r_1054', 'r_0452', 'r_0892', 'r_0893', 'r_1049']
    # for r in reactions:
    #     case10.checkReaction(r)

