
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
        super(Case10, self).__int__(self, cobra_model)


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

    def testO2EthanolProd (self, g_knockout = None, cs = 'r_1714', cs_lb = -1.5, range_o2 = list(np.arange(-1000, 0, 100))):
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



# if __name__ == '__main__':

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
wt_etOH = case14.testO2EthanolProd()

# ====== WILD TYPE ======
#FBA
wt_fba_res, wt_fba_exp_sim, wt_fba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,0], type = 'fba', res_exists = False, fname = 'Results/Case 14/res_fba_wt_case14.sav')
case14.plotExpVsSim(wt_fba_exp_sim_errors, save_fig_path = 'Results/Case 14/wt_fba_exp_sim_plot.png', title = 'FBA Wild Type')
plt.close('all')

#pFBA
wt_pfba_res, wt_pfba_exp_sim, wt_pfba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,0], type = 'pfba', res_exists = False, fname = 'Results/Case 14/res_pfba_wt_case14.sav')
case14.plotExpVsSim(wt_pfba_exp_sim_errors, save_fig_path = 'Results/Case 14/wt_pfba_exp_sim_plot.png', title = 'pFBA Wild Type')
plt.close('all')

#FVA
wt_fva_res, wt_fva_exp_sim, _ = case14.simulationPipeline(exp_dataset.ix[:,0], type = 'fva', res_exists = False, fname = 'Results/Case 14/res_fva_wt_case14.sav')



#
#
# # TESTS
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
#
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
