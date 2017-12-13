
'''
Functions for case 3 simulations

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


class Case3 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case3, self).__int__(self, cobra_model)

    def dictsForCase3 (self):
        # Carbon source lbs
        cs_lb = {'glucose': -1.5, 'mannose': -1.5, 'galactose': -1.5, 'pyruvate': -4.5} #Shophia's values

        #NO ETHANOL FLUX VALUES TO ESTIMATE O2 FLUXES

        # Carbon source exchange reactions
        cs_reaction = {'glucose': 'r_1714', 'mannose': 'r_1715', 'galactose': 'r_1710', 'pyruvate': 'r_2033'}

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
        sub_exp_fluxes.loc['r_4041'] = [0.33, 0.32, 0.20, 0.10] #growth rate

        return sub_exp_fluxes, reactions

    def getColumnWithoutNAs (self, dataframe, column_index, na = 'x'):
        df = dataframe[dataframe.ix[:, column_index] != na].ix[:, column_index]

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
            res = self.singleSimulation(carbon_source = self.cs_reaction[cs], cs_lb = self.cs_lb[cs], geneko = geneko, type = type)
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
case3 = Case3()
case3.model = case3.loadObjectFromFile('model_yeast_76.sav')
case3.model.solver = 'optlang-cplex'
case3.setMedium('MINIMAL')
case3.dictsForCase3()

#General datasets
exp_dataset, reactions = case3.loadExperimentalRes('Results/Case 3/case3_experimental_fluxes.csv')

# ====== CS: GLUCOSE ======
g_exp_df = case3.getColumnWithoutNAs(exp_dataset, 0)

#FBA
g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case3.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_glucose_case3.sav')
case3.plotExpVsSim(g_fba_exp_sim_errors, save_fig_path = 'Results/Case 3/g_fba_exp_sim_plot.png', title = 'FBA GLucose Carbon Source')
plt.close('all')

#pFBA
g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case3.simulationPipeline(g_exp_df, cs = 'glucose', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_glucose_case3.sav')
case3.plotExpVsSim(g_pfba_exp_sim_errors, save_fig_path = 'Results/Case 3/g_pfba_exp_sim_plot.png', title = 'pFBA GLucose Carbon Source')
plt.close('all')

#FVA
g_fva_res, g_fva_exp_sim, _ = case3.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_glucose_case3.sav')


# ====== CS: MANNOSE ======
m_exp_df = case3.getColumnWithoutNAs(exp_dataset, 1)

#FBA
m_fba_res, m_fba_exp_sim, m_fba_exp_sim_errors = case3.simulationPipeline(m_exp_df, cs = 'mannose', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_mannose_case10.sav')
case3.plotExpVsSim(m_fba_exp_sim_errors, save_fig_path = 'Results/Case 3/m_fba_exp_sim_plot.png', title = 'FBA Mannose Carbon Source')
plt.close('all')

#pFBA
m_pfba_res, m_pfba_exp_sim, m_pfba_exp_sim_errors = case3.simulationPipeline(m_exp_df, cs = 'mannose', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_mannose_case10.sav')
case3.plotExpVsSim(m_pfba_exp_sim_errors, save_fig_path = 'Results/Case 3/m_pfba_exp_sim_plot.png', title = 'pFBA Mannose Carbon Source')
plt.close('all')

#FVA
m_fva_res, m_fva_exp_sim, _ = case3.simulationPipeline(m_exp_df, cs = 'mannose', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_mannose_case10.sav')


# ====== CS: GALACTOSE ======
gal_exp_df = case3.getColumnWithoutNAs(exp_dataset, 2)

#FBA
gal_fba_res, gal_fba_exp_sim, gal_fba_exp_sim_errors = case3.simulationPipeline(gal_exp_df, cs = 'galactose', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_galactose_case3.sav')
case3.plotExpVsSim(gal_fba_exp_sim_errors, save_fig_path = 'Results/Case 3/gal_fba_exp_sim_plot.png', title = 'FBA Galactose Carbon Source')
plt.close('all')

#pFBA
gal_pfba_res, gal_pfba_exp_sim, gal_pfba_exp_sim_errors = case3.simulationPipeline(gal_exp_df, cs = 'galactose', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_galactose_case3.sav')
case3.plotExpVsSim(gal_pfba_exp_sim_errors, save_fig_path = 'Results/Case 3/gal_pfba_exp_sim_plot.png', title = 'pFBA Galactose Carbon Source')
plt.close('all')

#FVA
gal_fva_res, gal_fva_exp_sim, _ = case3.simulationPipeline(gal_exp_df, cs = 'galactose', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_galactose_case3.sav')


# ====== CS: PYRUVATE ======
p_exp_df = case3.getColumnWithoutNAs(exp_dataset, 3)

#FBA
p_fba_res, p_fba_exp_sim, p_fba_exp_sim_errors = case3.simulationPipeline(p_exp_df, cs = 'pyruvate', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_pyruvate_case3.sav')
case3.plotExpVsSim(p_fba_exp_sim_errors, save_fig_path = 'Results/Case 3/p_fba_exp_sim_plot.png', title = 'FBA Pyruvate Carbon Source')
plt.close('all')

#pFBA
p_pfba_res, p_pfba_exp_sim, p_pfba_exp_sim_errors = case3.simulationPipeline(p_exp_df, cs = 'pyruvate', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_pyruvate_case3.sav')
case3.plotExpVsSim(p_pfba_exp_sim_errors, save_fig_path = 'Results/Case 3/p_pfba_exp_sim_plot.png', title = 'pFBA Pyruvate Carbon Source')
plt.close('all')

#FVA
p_fva_res, p_fva_exp_sim, _ = case3.simulationPipeline(p_exp_df, cs = 'pyruvate', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_pyruvate_case3.sav')




# Acetyl-CoA  <==> Acetyl-CoA-mit (transport) no correspondence
