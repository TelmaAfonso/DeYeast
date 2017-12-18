
'''
Functions for case 13 simulations

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


class Case13 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case13, self).__int__(self, cobra_model)

    def dictsForCase13 (self):
        # Carbon source lb
        cs_lb = {'g_oxidative': -1.56, 'g_resp_fermentative': -4.9, 'g_fermentative': -8.23}

        # Carbon source exchange reaction
        cs_reaction = {'g_oxidative': 'r_1714', 'g_resp_fermentative': 'r_1714', 'g_fermentative': 'r_1714'}

        self.cs_lb, self.cs_reaction = cs_lb, cs_reaction

    def loadExperimentalRes (self, filename_exp, sep = ';'):
        exp_fluxes = pd.read_csv(filename_exp, sep = sep)
        exp_fluxes.set_index('Yeast7_ID', inplace = True)
        reactions = exp_fluxes.iloc[:,0]
        sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if 'exp flux' not in col], 1)
        sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.'))
        #sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.astype('float') if str(x).isdigit() else None)

        #Add Biomass Values (r_4041)
        reactions['r_4041'] = 'Biomass'
        sub_exp_fluxes.loc['r_4041'] = [0.15, 0.30, 0.40] #growth rate


        return sub_exp_fluxes, reactions

    def getColumnWithoutNAs (self, dataframe, column_index, na = 'x'):
        df = dataframe[dataframe.ix[:, column_index] != na].ix[:, column_index]

        return df.astype('float')

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

    def simulationPipeline (self, exp_dataset, cs = 'g_oxidative', o2_lb = None, geneko = None, type = 'fba', res_exists = False, fname = None):
        if res_exists:
            res = self.loadObjectFromFile(fname)
        else:
            res = self.singleSimulation(carbon_source = self.cs_reaction[cs], cs_lb = self.cs_lb[cs], geneko = geneko, o2_lb = o2_lb, type = type)
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


#if __name__ == '__main__':

#Initialization
case13 = Case13()
case13.model = case13.loadObjectFromFile('model_yeast_76.sav')
case13.model.solver = 'optlang-cplex'
case13.setMedium('MINERAL')
case13.dictsForCase13()

#General datasets
exp_dataset, reactions = case13.loadExperimentalRes('Results/Case 13/case13_experimental_fluxes.csv')

# ====== OXIDATIVE GROWTH ======
o_exp_df = case13.getColumnWithoutNAs(exp_dataset, 0, 'X')
# NO EtOH fluxes available for O2 flux estimation

#FBA
o_fba_res, o_fba_exp_sim, o_fba_exp_sim_errors = case13.simulationPipeline(o_exp_df, cs = 'g_oxidative', type = 'fba', res_exists = True, fname = 'Results/Case 13/res_fba_oxidative_case13.sav')
o_fba_exp_sim_errors = case13.getDFWithoutExtremeFluxes(o_fba_exp_sim_errors) #without extreme fluxes (for plotting)
case13.plotExpVsSim(o_fba_exp_sim_errors, save_fig_path = 'Results/Case 13/o_fba_exp_sim_plot.png', title = 'FBA Oxidative Growth')
plt.close('all')

#pFBA
o_pfba_res, o_pfba_exp_sim, o_pfba_exp_sim_errors = case13.simulationPipeline(o_exp_df, cs = 'g_oxidative',type = 'pfba', res_exists = True, fname = 'Results/Case 13/res_pfba_oxidative_case13.sav')
case13.plotExpVsSim(o_pfba_exp_sim_errors, save_fig_path = 'Results/Case 13/o_pfba_exp_sim_plot.png', title = 'pFBA Oxidative Growth')
plt.close('all')

#FVA
o_fva_res, o_fva_exp_sim, _ = case13.simulationPipeline(o_exp_df, cs = 'g_oxidative', type = 'fva', res_exists = True, fname = 'Results/Case 13/res_fva_oxidative_case13.sav')


# ====== RESPIRO-FERMENTATIVE GROWTH ======
rf_exp_df = case13.getColumnWithoutNAs(exp_dataset, 1, 'X')
# NO EtOH fluxes available for O2 flux estimation

#FBA
rf_fba_res, rf_fba_exp_sim, rf_fba_exp_sim_errors = case13.simulationPipeline(rf_exp_df, cs = 'g_resp_fermentative', type = 'fba', res_exists = True, fname = 'Results/Case 13/res_fba_resp_fermentative_case13.sav')
rf_fba_exp_sim_errors = case13.getDFWithoutExtremeFluxes(rf_fba_exp_sim_errors) #without extreme fluxes (for plotting)
case13.plotExpVsSim(rf_fba_exp_sim_errors, save_fig_path = 'Results/Case 13/rf_fba_exp_sim_plot.png', title = 'FBA Respiro-Fermentative Growth')
plt.close('all')

#pFBA
rf_pfba_res, rf_pfba_exp_sim, rf_pfba_exp_sim_errors = case13.simulationPipeline(rf_exp_df, cs = 'g_resp_fermentative',type = 'pfba', res_exists = True, fname = 'Results/Case 13/res_pfba_resp_fermentative_case13.sav')
case13.plotExpVsSim(rf_pfba_exp_sim_errors, save_fig_path = 'Results/Case 13/rf_pfba_exp_sim_plot.png', title = 'pFBA Respiro-Fermentative Growth')
plt.close('all')

#FVA
rf_fva_res, rf_fva_exp_sim, _ = case13.simulationPipeline(rf_exp_df, cs = 'g_resp_fermentative', type = 'fva', res_exists = True, fname = 'Results/Case 13/res_fva_resp_fermentative_case13.sav')


# ====== FERMENTATIVE GROWTH ======
f_exp_df = case13.getColumnWithoutNAs(exp_dataset, 2, 'X')
# NO EtOH fluxes available for O2 flux estimation

#FBA
f_fba_res, f_fba_exp_sim, f_fba_exp_sim_errors = case13.simulationPipeline(f_exp_df, cs = 'g_fermentative', type = 'fba', res_exists = True, fname = 'Results/Case 13/res_fba_fermentative_case13.sav')
f_fba_exp_sim_errors = case13.getDFWithoutExtremeFluxes(f_fba_exp_sim_errors) #without extreme fluxes (for plotting)
case13.plotExpVsSim(f_fba_exp_sim_errors, save_fig_path = 'Results/Case 13/f_fba_exp_sim_plot.png', title = 'FBA Fermentative Growth')
plt.close('all')

#pFBA
f_pfba_res, f_pfba_exp_sim, f_pfba_exp_sim_errors = case13.simulationPipeline(f_exp_df, cs = 'g_fermentative',type = 'pfba', res_exists = True, fname = 'Results/Case 13/res_pfba_fermentative_case13.sav')
case13.plotExpVsSim(f_pfba_exp_sim_errors, save_fig_path = 'Results/Case 13/f_pfba_exp_sim_plot.png', title = 'pFBA Fermentative Growth')
plt.close('all')

#FVA
f_fva_res, f_fva_exp_sim, _ = case13.simulationPipeline(f_exp_df, cs = 'g_fermentative', type = 'fva', res_exists = True, fname = 'Results/Case 13/res_fva_fermentative_case13.sav')


