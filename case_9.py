
'''
Functions for case 9 simulations

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


class Case9 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case9, self).__int__(self, cobra_model)

    def dictsForCase9 (self):
        # Carbon source lb
        # cs_lb = {'glucose': -1.5} #Shophia's values
        cs_lb = {'batch': -15.9, 'chemostat': -1.17} #Authors's values

        # Carbon source exchange reaction
        cs_reaction = {'batch': 'r_1714', 'chemostat': 'r_1714'}

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
        sub_exp_fluxes.loc['r_4041'] = [0.37, 0.1] #growth rate

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

    def testO2EthanolProd (self, g_knockout = None, cs = 'glucose', range_o2 = list(np.arange(-20, 0, 2))):
        res = {}
        for i in range_o2:
            with self.model as m:
                m.set_carbon_source(self.cs_reaction[cs], lb = self.cs_lb[cs])
                m.reactions.get_by_id('r_1992').lower_bound = float(i)
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


if __name__ == '__main__':

    #Initialization
    case9 = Case9()
    case9.model = case9.loadObjectFromFile('model_yeast_76.sav')
    case9.model.solver = 'optlang-cplex'
    case9.setMedium('MINIMAL')
    case9.dictsForCase9()

    #General datasets
    exp_dataset, reactions = case9.loadExperimentalRes('Results/Case 9/case9_experimental_fluxes.csv')

    # ====== BATCH ======
    b_exp_df = case9.getColumnWithoutNAs(exp_dataset, 0, 'X')

    # O2 FLUX ESTIMATION - ALL ZERO
    # b_etOH = case9.testO2EthanolProd(cs = 'batch', range_o2 = list(np.arange(-2, 0, 0.1)))
    # case9.saveObjectToFile(b_etOH, 'Results/Case 9/b_dict_etOH_O2_fluxes.sav')
    # b_etOH = case9.loadObjectFromFile('Results/Case 9/b_dict_etOH_O2_fluxes.sav')
    # b_o2_lb = case9.plotO2vsEtOH(b_etOH, real_EtOH_flux = 2.385, legend = 'Batch Culture', fname = 'Results/Case 9/b_etOH_plot.png')
    # plt.close('all')

    #FBA
    b_fba_res, b_fba_exp_sim, b_fba_exp_sim_errors = case9.simulationPipeline(b_exp_df, cs = 'batch', type = 'fba', res_exists = True, fname = 'Results/Case 9/res_fba_batch_case9.sav')
    b_fba_exp_sim_errors = case9.getDFWithoutExtremeFluxes(b_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case9.plotExpVsSim(b_fba_exp_sim_errors, save_fig_path = 'Results/Case 9/b_fba_exp_sim_plot.png', title = 'FBA Batch')
    plt.close('all')

    #pFBA
    b_pfba_res, b_pfba_exp_sim, b_pfba_exp_sim_errors = case9.simulationPipeline(b_exp_df, cs = 'batch',type = 'pfba', res_exists = True, fname = 'Results/Case 9/res_pfba_batch_case9.sav')
    case9.plotExpVsSim(b_pfba_exp_sim_errors, save_fig_path = 'Results/Case 9/b_pfba_exp_sim_plot.png', title = 'pFBA Batch')
    plt.close('all')

    #FVA
    b_fva_res, b_fva_exp_sim, _ = case9.simulationPipeline(b_exp_df, cs = 'batch', type = 'fva', res_exists = True, fname = 'Results/Case 9/res_fva_batch_case9.sav')


    # ====== CHEMOSTAT ======
    c_exp_df = case9.getColumnWithoutNAs(exp_dataset, 1, 'X')

    # O2 FLUX ESTIMATION - ALL ZERO
    # c_etOH = case9.testO2EthanolProd(cs = 'chemostat', range_o2 = list(np.arange(-2, 0, 0.1)))
    # case9.saveObjectToFile(c_etOH, 'Results/Case 9/c_dict_etOH_O2_fluxes.sav')
    # c_etOH = case9.loadObjectFromFile('Results/Case 9/c_dict_etOH_O2_fluxes.sav')
    # c_o2_lb = case9.plotO2vsEtOH(c_etOH, real_EtOH_flux = 2.385, legend = 'Batch Culture', fname = 'Results/Case 9/c_etOH_plot.png')
    # plt.close('all')

    #FBA
    c_fba_res, c_fba_exp_sim, c_fba_exp_sim_errors = case9.simulationPipeline(c_exp_df, cs = 'chemostat', type = 'fba', res_exists = True, fname = 'Results/Case 9/res_fba_chemostat_case9.sav')
    c_fba_exp_sim_errors = case9.getDFWithoutExtremeFluxes(c_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case9.plotExpVsSim(c_fba_exp_sim_errors, save_fig_path = 'Results/Case 9/c_fba_exp_sim_plot.png', title = 'FBA Chemostat')
    plt.close('all')

    #pFBA
    c_pfba_res, c_pfba_exp_sim, c_pfba_exp_sim_errors = case9.simulationPipeline(c_exp_df, cs = 'chemostat',type = 'pfba', res_exists = True, fname = 'Results/Case 9/res_pfba_chemostat_case9.sav')
    case9.plotExpVsSim(c_pfba_exp_sim_errors, save_fig_path = 'Results/Case 9/c_pfba_exp_sim_plot.png', title = 'pFBA Chemostat')
    plt.close('all')

    #FVA
    c_fva_res, c_fva_exp_sim, _ = case9.simulationPipeline(c_exp_df, cs = 'chemostat', type = 'fva', res_exists = True, fname = 'Results/Case 9/res_fva_chemostat_case9.sav')
