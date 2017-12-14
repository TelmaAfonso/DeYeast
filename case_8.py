
'''
Functions for case 8 simulations

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


class Case8 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case8, self).__int__(self, cobra_model)

    def dictsForCase8 (self):
        # Carbon source lbs
        cs_lb = {'glucose': -1.5, 'galactose': -1.5, 'glycerol': -2.5, 'ethanol': -5} #Shophia's values

        # Carbon source exchange reactions
        cs_reaction = {'glucose': 'r_1714', 'galactose': 'r_1710', 'glycerol': 'r_1808', 'ethanol': 'r_1761'}

        self.cs_lb, self.cs_reaction = cs_lb, cs_reaction

    def loadExperimentalRes (self, filename_exp, sep = ';'):
        exp_fluxes = pd.read_csv(filename_exp, sep = sep)
        exp_fluxes.set_index('Yeast7_ID', inplace = True)
        reactions = exp_fluxes.iloc[:,0]
        sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if 'exp flux' not in col], 1)
        sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.'))
        #sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.astype('float') if str(x).isdigit() else None)

        #Add Biomass Values (r_4041)
        # reactions['r_4041'] = 'Biomass'
        # sub_exp_fluxes.loc['r_4041'] = [0.33, 0.32, 0.20, 0.10] #growth rate

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

    def simulationPipeline (self, exp_dataset, cs = 'glucose', o2_lb = None, geneko = None, type = 'fba', res_exists = False, fname = None):
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

    def testO2EthanolProd (self, g_knockout = None, cs = 'r_1714', cs_lb = -1.5, range_o2 = list(np.arange(-20, 0, 2))):
        loading_bars = 40*'='
        res = {}
        for i in range_o2:
            with self.model as m:
                m.set_carbon_source(cs, lb = cs_lb)
                m.reactions.get_by_id('r_2115').lower_bound = float(i)
                if g_knockout is not None:
                    m.set_environmental_conditions(gene_knockout = g_knockout)
                r = pfba(m)
                if hasattr(r, 'x_dict'): #if legacy solution
                    fluxes = pd.DataFrame(list(r.x_dict.items())).set_index(0)
                else:
                    fluxes = r.fluxes
                res[str(i)] = fluxes.loc['r_0163']
        for key, val in sorted(res.items()): print(key, '\t', val)
        return res

    def plotO2vsEtOH (self, dict_EtOH_res, real_EtOH_flux = 0, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', legend = 'Wild Type', fname = None):
        plt.figure(figsize = (10, 5))
        x = sorted([float(x) for x in dict_EtOH_res.keys()])
        y = [float(dict_EtOH_res[str(key)]) for key in x]
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        real_O2 = lambda x0: (y0 - intercept) / slope
        y0 = real_EtOH_flux
        plt.plot(x, y, 'o', x, line)
        plt.axhline(y = real_EtOH_flux, ls = 'dashed')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.legend([legend, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f (O2 flux of %.2f)' % (real_EtOH_flux, real_O2(y0))])
        plt.show()

        if fname is not None:
            plt.savefig(fname)

        return round(real_O2(y0), 4)


#if __name__ == '__main__':

#Strain FY4 is derived from S288C (very similar)

#Initialization
case8 = Case8()
case8.model = case8.loadObjectFromFile('model_yeast_76.sav')
case8.model.solver = 'optlang-cplex'
case8.setMedium('MINIMAL')
case8.dictsForCase8()

#S. Cerevisiae BY4741 deltas
genes = ['HIS3', 'LEU2', 'MET17', 'URA3']
genes = list(case8.convertStdToSyst(genes).values())

#General datasets
exp_dataset, reactions = case8.loadExperimentalRes('Results/Case 8/case8_experimental_fluxes.csv')

# ====== CS: GLUCOSE ======
g_exp_df = case8.getColumnWithoutNAs(exp_dataset, 0, 'X')

# O2 FLUX ESTIMATION - ALL ZERO!
# g_etOH = case5.testO2EthanolProd(range_o2 = list(np.arange(-10, 0, 1)))
# case5.saveObjectToFile(g_etOH, 'Results/Case 5/g_dict_etOH_O2_fluxes.sav')
# g_etOH = case5.loadObjectFromFile('Results/Case 5/g_dict_etOH_O2_fluxes.sav')
# g_o2_lb = case5.plotO2vsEtOH(g_etOH, real_EtOH_flux = -0.2478, fname = 'Results/Case 5/g_etOH_plot.png')
# plt.close('all')

#FBA
g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case5.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 5/res_fba_glucose_case5.sav')
g_fba_exp_sim_errors = case5.getDFWithoutExtremeFluxes(g_fba_exp_sim_errors) #without r_0302, for plotting
case5.plotExpVsSim(g_fba_exp_sim_errors, save_fig_path = 'Results/Case 5/g_fba_exp_sim_plot.png', title = 'FBA Glucose Carbon Source')
plt.close('all')

#pFBA
g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case5.simulationPipeline(g_exp_df, cs = 'glucose', type = 'pfba', res_exists = True, fname = 'Results/Case 5/res_pfba_glucose_case5.sav')
case5.plotExpVsSim(g_pfba_exp_sim_errors, save_fig_path = 'Results/Case 5/g_pfba_exp_sim_plot.png', title = 'pFBA Glucose Carbon Source')
plt.close('all')

#FVA
g_fva_res, g_fva_exp_sim, _ = case5.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 5/res_fva_glucose_case5.sav')


# R03668 - KEGGID Wrong in cecafdb dataset
# D-Glyceraldehyde-3-phosphate <==> Glycerol (transport) not present in model
# Acetyl-CoA <==> Acetyl-CoA-mit	(transport) not present in model


case5.checkReaction(case5.convertKeggID('R03668'))

