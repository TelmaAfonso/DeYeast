
'''
Functions for case 14 simulations

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


class Case14 (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(Case14, self).__int__(self, cobra_model)


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

    def simulationPipeline (self, exp_dataset, cs = 'r_1714', cs_lb = -1.5, o2_lb = None, geneko = None, type = 'fba', res_exists = False, fname = None):
        if res_exists:
            res = self.loadObjectFromFile(fname)
        else:
            res = self.singleSimulation(carbon_source = cs, cs_lb = cs_lb, geneko = geneko, o2_lb = o2_lb, type = type)
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

    def testO2EthanolProd (self, g_knockout = None, cs = 'r_1714', cs_lb = -1.5, range_o2 = list(np.arange(-20, 0, 2))):
        loading_bars = 40*'='
        res = {}
        for i in range_o2:
            with self.model as m:
                m.set_carbon_source(cs, lb = cs_lb)
                m.reactions.get_by_id('r_1992').lower_bound = float(i)
                if g_knockout is not None:
                    m.set_environmental_conditions(gene_knockout = g_knockout)
                r = pfba(m)
                if hasattr(r, 'x_dict'): #if legacy solution
                    fluxes = pd.DataFrame(list(r.x_dict.items())).set_index(0)
                else:
                    fluxes = r.fluxes
                res[str(i)] = fluxes.loc['r_2115']
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
    case14 = Case14()
    case14.model = case14.loadObjectFromFile('model_yeast_76.sav')
    case14.model.solver = 'optlang-cplex'
    case14.setMedium('MINIMAL')

    genes = ['REG1', 'MIG1', 'MIG2', 'GRR1', 'HXK2'] #Knockouts in this study
    HXK2 = case14.convertStdToSyst(genes)['HXK2'] # Gene match only for HXK2 gene
    # no reaction in the model for Acetyl-CoA <==> Acetyl-CoA-mit (transport) and Pyruvate-mit + CO2-mit <==> Oxaloacetate-mit (R00217-M)

    #General datasets
    exp_dataset, reactions = case14.loadExperimentalRes('Results/Case 14/case14_experimental_fluxes.csv')


    # ====== WILD TYPE ======

    # O2 FLUX ESTIMATION
    # wt_etOH = case14.testO2EthanolProd(range_o2 = list(np.arange(-2, 0, 0.2)))
    # case14.saveObjectToFile(wt_etOH2, 'Results/Case 14/wt_dict_etOH_O2_fluxes.sav')
    wt_etOH = case14.loadObjectFromFile('Results/Case 14/wt_dict_etOH_O2_fluxes.sav')
    wt_o2_lb = case14.plotO2vsEtOH(wt_etOH, real_EtOH_flux = 2.1080, fname = 'Results/Case 14/wt_etOH_plot.png')
    plt.close('all')


    #FBA
    wt_fba_res, wt_fba_exp_sim, wt_fba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,0], o2_lb = wt_o2_lb, type = 'fba', res_exists = True, fname = 'Results/Case 14/res_fba_wt_case14.sav')
    wt_fba_exp_sim_errors = case14.getDFWithoutExtremeFluxes(wt_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case14.plotExpVsSim(wt_fba_exp_sim_errors, save_fig_path = 'Results/Case 14/wt_fba_exp_sim_plot.png', title = 'FBA Wild Type')
    plt.close('all')

    #pFBA
    wt_pfba_res, wt_pfba_exp_sim, wt_pfba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,0], o2_lb = wt_o2_lb, type = 'pfba', res_exists = True, fname = 'Results/Case 14/res_pfba_wt_case14.sav')
    case14.plotExpVsSim(wt_pfba_exp_sim_errors, save_fig_path = 'Results/Case 14/wt_pfba_exp_sim_plot.png', title = 'pFBA Wild Type')

    plt.close('all')

    #FVA
    wt_fva_res, wt_fva_exp_sim, _ = case14.simulationPipeline(exp_dataset.ix[:,0], o2_lb = wt_o2_lb, type = 'fva', res_exists = True, fname = 'Results/Case 14/res_fva_wt_case14.sav')



    # ====== HXK2 DELETION ======

    # O2 FLUX ESTIMATION
    # hxk2_etOH = case14.testO2EthanolProd(g_knockout = HXK2,range_o2 = list(np.arange(-2, 0, 0.2)))
    # case14.saveObjectToFile(hxk2_etOH, 'Results/Case 14/hxk2_dict_etOH_O2_fluxes.sav')
    hxk2_etOH = case14.loadObjectFromFile('Results/Case 14/hxk2_dict_etOH_O2_fluxes.sav')
    hxk2_02_lb = case14.plotO2vsEtOH(hxk2_etOH, real_EtOH_flux = 1.4221, fname = 'Results/Case 14/hxk2_etOH_plot.png')
    plt.close('all')


    #FBA
    hxk2_fba_res, hxk2_fba_exp_sim, hxk2_fba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, type = 'fba', res_exists = True, fname = 'Results/Case 14/res_fba_hxk2_case14.sav')
    hxk2_fba_exp_sim_errors = case14.getDFWithoutExtremeFluxes(hxk2_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case14.plotExpVsSim(hxk2_fba_exp_sim_errors, save_fig_path = 'Results/Case 14/hxk2_fba_exp_sim_plot.png', title = 'FBA HXK2 Del')
    plt.close('all')

    #pFBA
    hxk2_pfba_res, hxk2_pfba_exp_sim, hxk2_pfba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, type = 'pfba', res_exists = True, fname = 'Results/Case 14/res_pfba_hxk2_case14.sav')
    case14.plotExpVsSim(hxk2_pfba_exp_sim_errors, save_fig_path = 'Results/Case 14/hxk2_pfba_exp_sim_plot.png', title = 'pFBA HXK2 Del')
    plt.close('all')

    #FVA
    hxk2_fva_res, hxk2_fva_exp_sim, _ = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, type = 'fva', res_exists = True, fname = 'Results/Case 14/res_fva_hxk2_case14.sav')

    #LMOMA
    hxk2 = case14.convertStdToSyst(['HXK2'])['HXK2']
    hxk2_lmoma_res, hxk2_lmoma_exp_sim, hxk2_lmoma_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, geneko = hxk2, type = 'lmoma', res_exists = True, fname = 'Results/Case 14/res_lmoma_hxk2_case14.sav')
    case14.plotExpVsSim(hxk2_lmoma_exp_sim_errors, save_fig_path = 'Results/Case 14/hxk2_lmoma_exp_sim_plot.png', title = 'LMOMA hxk2 Del')
    plt.close('all')


    # SEE r_0962 signal in exp dataset

    # error = case14.rmse(wt_pfba_exp_sim)
    # r_sqrd = case14.r_squared(wt_pfba_exp_sim)


    # =========================================
    #    Save all results into a binary file
    # =========================================

    all_res = {'d14_wt_fba': wt_fba_exp_sim, 'd14_wt_pfba': wt_pfba_exp_sim, 'd14_wt_fva': wt_fva_exp_sim, 'd14_reactions': reactions,
               'd14_hxk2_fba': hxk2_fba_exp_sim, 'd14_hxk2_pfba': hxk2_pfba_exp_sim, 'd14_hxk2_fva': hxk2_fva_exp_sim, 'd14_hxk2_lmoma': hxk2_lmoma_exp_sim}

    case14.saveObjectToFile(all_res, 'Results/case14_all_res.sav')
    
    

# TESTS
# import time
# import sys
#
# toolbar_width = 40
#
# # setup toolbar
# sys.stdout.write("[%s]" % (" " * toolbar_width))
# sys.stdout.flush()
# sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
#
# for i in range(toolbar_width):
#     time.sleep(0.1) # do real work here
#     # update the bar
#     sys.stdout.write("=")
#     sys.stdout.flush()
#
# sys.stdout.write("\n")
#
#
# import time
#
# from tqdm import *
#
# for i in tqdm(range(10)):
#     time.sleep(3)
#
#
# from progress.bar import Bar
#
# bar = Bar('Processing', max=20)
# for i in range(20):
#     print('Do stuff')
#     bar.next()
# bar.finish()
# #
#
#
# import progressbar
# import time
#
# progress = progressbar.ProgressBar()
# for i in progress(range(80)):
#     time.sleep(0.01)
#
#
#
# import pyprind
# import sys
# import time
#
# n = 100
# bar = pyprind.ProgBar(n, stream=sys.stdout)
# for i in range(n):
#     time.sleep(0.1)
# bar.update()
#
