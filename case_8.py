
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

    # def plotExpVsSim (self, absRelErrorDataset, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'Wild Type', label_adjust = 0.05, save_fig_path = None):
    #     plt.rcParams["figure.figsize"] = (10,5)
    #
    #     x = absRelErrorDataset.ix[:,0]
    #     y = absRelErrorDataset.ix[:,1]
    #     react_IDs = list(absRelErrorDataset.index)
    #     slope, intercept, r_value, p_value, std_err = linregress(x, y)
    #     line = [slope * x + intercept for x in x]
    #     meanRelErr = absRelErrorDataset.ix[:,3].mean()
    #     corr = x.corr(y)
    #
    #     plt.plot(x, y, 'o', x, line)
    #     for ind, react_ID in enumerate(react_IDs):
    #         plt.annotate(react_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + label_adjust, y[ind] + label_adjust))
    #
    #     plt.ylabel(ylab)
    #     plt.xlabel(xlab)
    #     plt.title(title)
    #     plt.plot([], [], ' ') # To show correlation in legend
    #     plt.plot([], [], ' ') # To show mean relative error in legend
    #     plt.legend(['Reactions', 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])
    #
    #     if save_fig_path is not None:
    #         plt.savefig(save_fig_path)

    def testO2EthanolProd (self, g_knockout = None, react_id = 'r_2115', cs = 'glucose', range_o2 = list(np.arange(-20, 0, 2))):
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
                res[str(i)] = fluxes.loc[react_id]
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
    case8 = Case8()
    case8.model = case8.loadObjectFromFile('model_yeast_76.sav')
    case8.model.solver = 'optlang-cplex'
    case8.setMedium('MINIMAL')
    case8.dictsForCase8()

    #S. Cerevisiae BY4741 deltas
    # genes = ['HIS3', 'LEU2', 'MET17', 'URA3']
    # genes = list(case8.convertStdToSyst(genes).values())

    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC127579/
    # https://www.yeastgenome.org/strain/BY4741#overview

    # Some of the most commonly applied marker genes are wild-type alleles of yeast genes that encode key enzymes in
    # the metabolic pathways towards essential monomers used in biosynthesis. An example is the URA3 gene, which encodes
    # orotidine-5′-phosphate decarboxylase, an essential enzyme in pyrimidine biosynthesis in Saccharomyces cerevisiae (3).
    # Similarly, the HIS3, LEU2, TRP1, and MET15 (=MET17) marker genes encode essential enzymes for de novo synthesis of the amino
    # acids l-histidine, l-leucine, l-tryptophan, and l-methionine, respectively (4, 8).

    # BY4742 (MAT a his3D1 leu2 D0 lys2 D0 ura3 D0) represents the haploid strain BY4742, which is of the a mating type
    # and auxotrophic for histidine, leucine, lysine, and uracil biosynthetic pathways.

    #General datasets
    exp_dataset, reactions = case8.loadExperimentalRes('Results/Case 8/case8_experimental_fluxes.csv')

    # ====== CS: GLUCOSE ======
    g_exp_df = case8.getColumnWithoutNAs(exp_dataset, 0, 'X')

    # O2 FLUX ESTIMATION - O2 flux of -0.43
    # g_etOH = case8.testO2EthanolProd(g_knockout = None, cs = 'glucose', range_o2 = list(np.arange(-2, 0, 0.1)))
    # case8.saveObjectToFile(g_etOH, 'Results/Case 8/g_dict_etOH_O2_fluxes.sav')
    # g_etOH = case8.loadObjectFromFile('Results/Case 8/g_dict_etOH_O2_fluxes.sav')
    # g_o2_lb = case8.plotO2vsEtOH(g_etOH, real_EtOH_flux = 2.115, legend = 'FY4 - Glucose', fname = 'Results/Case 8/g_etOH_plot.png')
    # plt.close('all')

    #FBA
    g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case8.simulationPipeline(g_exp_df, o2_lb = -0.43, cs = 'glucose', geneko = None, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_glucose_case5.sav')
    g_fba_exp_sim_errors = case8.getDFWithoutExtremeFluxes(g_fba_exp_sim_errors) #without extreme fluxes (for plotting)
    case8.plotExpVsSim(g_fba_exp_sim_errors, save_fig_path = 'Results/Case 8/g_fba_exp_sim_plot.png', title = 'FBA Glucose Carbon Source')
    plt.close('all')

    #pFBA
    g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case8.simulationPipeline(g_exp_df, o2_lb = -0.43, cs = 'glucose', geneko = None, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_glucose_case5.sav')
    case8.plotExpVsSim(g_pfba_exp_sim_errors, save_fig_path = 'Results/Case 8/g_pfba_exp_sim_plot.png', title = 'pFBA Glucose Carbon Source')
    plt.close('all')

    # #LMOMA
    # g_lmoma_res, g_lmoma_exp_sim, g_lmoma_exp_sim_errors = case8.simulationPipeline(g_exp_df, cs = 'glucose', geneko = None, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_glucose_case5.sav')
    # case8.plotExpVsSim(g_lmoma_exp_sim_errors, save_fig_path = 'Results/Case 8/g_lmoma_exp_sim_plot.png', title = 'LMOMA Glucose Carbon Source')
    # plt.close('all')

    #FVA
    g_fva_res, g_fva_exp_sim, _ = case8.simulationPipeline(g_exp_df, cs = 'glucose', o2_lb = -0.43, geneko = None, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_glucose_case5.sav')



    # ====== CS: GALACTOSE ======
    gal_exp_df = case8.getColumnWithoutNAs(exp_dataset, 1, 'X')

    # O2 FLUX ESTIMATION - O2 flux of -0.96
    # gal_etOH = case8.testO2EthanolProd(g_knockout = None, cs = 'galactose', range_o2 = list(np.arange(-2, 0, 0.1)))
    # case8.saveObjectToFile(gal_etOH, 'Results/Case 8/gal_dict_etOH_O2_fluxes.sav')
    # gal_etOH = case8.loadObjectFromFile('Results/Case 8/gal_dict_etOH_O2_fluxes.sav')
    # gal_o2_lb = case8.plotO2vsEtOH(gal_etOH, real_EtOH_flux = 1.545, legend = 'FY4 - Galactose', fname = 'Results/Case 8/gal_etOH_plot.png')
    # plt.close('all')

    #FBA
    gal_fba_res, gal_fba_exp_sim, gal_fba_exp_sim_errors = case8.simulationPipeline(gal_exp_df, cs = 'galactose', o2_lb = -0.96, geneko = None, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_galactose_case5.sav')
    gal_fba_exp_sim_errors = case8.getDFWithoutExtremeFluxes(gal_fba_exp_sim_errors)
    case8.plotExpVsSim(gal_fba_exp_sim_errors, save_fig_path = 'Results/Case 8/gal_fba_exp_sim_plot.png', title = 'FBA Galactose Carbon Source')
    plt.close('all')

    #pFBA
    gal_pfba_res, gal_pfba_exp_sim, gal_pfba_exp_sim_errors = case8.simulationPipeline(gal_exp_df, cs = 'galactose', o2_lb = -0.96, geneko = None, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_galactose_case5.sav')
    case8.plotExpVsSim(gal_pfba_exp_sim_errors, save_fig_path = 'Results/Case 8/gal_pfba_exp_sim_plot.png', title = 'pFBA Galactose Carbon Source')
    plt.close('all')

    # #LMOMA
    # gal_lmoma_res, gal_lmoma_exp_sim, gal_lmoma_exp_sim_errors = case8.simulationPipeline(gal_exp_df, cs = 'galactose', geneko = None, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_galactose_case5.sav')
    # case8.plotExpVsSim(gal_lmoma_exp_sim_errors, save_fig_path = 'Results/Case 8/gal_lmoma_exp_sim_plot.png', title = 'LMOMA Galactose Carbon Source')
    # plt.close('all')

    #FVA
    gal_fva_res, gal_fva_exp_sim, _ = case8.simulationPipeline(gal_exp_df, cs = 'galactose', o2_lb = -0.96, geneko = None, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_galactose_case5.sav')


    # ====== CS: GLYCEROL ======
    gly_exp_df = case8.getColumnWithoutNAs(exp_dataset, 2, 'X')

    # O2 FLUX ESTIMATION - O2 flux of -0.75
    # gly_etOH = case8.testO2EthanolProd(g_knockout = None, cs = 'glycerol', range_o2 = list(np.arange(-2, 0, 0.1)))
    # case8.saveObjectToFile(gly_etOH, 'Results/Case 8/gly_dict_etOH_O2_fluxes.sav')
    # gly_etOH = case8.loadObjectFromFile('Results/Case 8/gly_dict_etOH_O2_fluxes.sav')
    # gly_o2_lb = case8.plotO2vsEtOH(gly_etOH, real_EtOH_flux = 0.175, legend = 'FY4 - Glycerol', fname = 'Results/Case 8/gly_etOH_plot.png')
    # plt.close('all')

    #FBA
    gly_fba_res, gly_fba_exp_sim, gly_fba_exp_sim_errors = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', o2_lb = -0.75, geneko = None, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_glycerol_case5.sav')
    gly_fba_exp_sim_errors = case8.getDFWithoutExtremeFluxes(gly_fba_exp_sim_errors)
    case8.plotExpVsSim(gly_fba_exp_sim_errors, save_fig_path = 'Results/Case 8/gly_fba_exp_sim_plot.png', title = 'FBA Glycerol Carbon Source')
    plt.close('all')

    #pFBA
    gly_pfba_res, gly_pfba_exp_sim, gly_pfba_exp_sim_errors = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', o2_lb = -0.75, geneko = None, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_glycerol_case5.sav')
    case8.plotExpVsSim(gly_pfba_exp_sim_errors, save_fig_path = 'Results/Case 8/gly_pfba_exp_sim_plot.png', title = 'pFBA Glycerol Carbon Source')
    plt.close('all')

    # case8.getListOfMetabolitesSummary(gly_pfba_res)
    # case8.getMetaboliteSummaryWithNames('s_0629', gly_pfba_res) #Glycerol (c)

    # #LMOMA
    # gly_lmoma_res, gly_lmoma_exp_sim, gly_lmoma_exp_sim_errors = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', geneko = None, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_glycerol_case5.sav')
    # case8.plotExpVsSim(gly_lmoma_exp_sim_errors, save_fig_path = 'Results/Case 8/gly_lmoma_exp_sim_plot.png', title = 'LMOMA Glycerol Carbon Source')
    # plt.close('all')

    #FVA
    gly_fva_res, gly_fva_exp_sim, _ = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', o2_lb = -0.75, geneko = None, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_glycerol_case5.sav')


    # ====== CS: ETHANOL ======
    e_exp_df = case8.getColumnWithoutNAs(exp_dataset, 3, 'X')

    #FBA
    e_fba_res, e_fba_exp_sim, e_fba_exp_sim_errors = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = None, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_ethanol_case5.sav')
    e_fba_exp_sim_errors = case8.getDFWithoutExtremeFluxes(e_fba_exp_sim_errors) #without r_0302, for plotting
    case8.plotExpVsSim(e_fba_exp_sim_errors, save_fig_path = 'Results/Case 8/e_fba_exp_sim_plot.png', title = 'FBA Ethanol Carbon Source')
    plt.close('all')

    #pFBA
    e_pfba_res, e_pfba_exp_sim, e_pfba_exp_sim_errors = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = None, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_ethanol_case5.sav')
    case8.plotExpVsSim(e_pfba_exp_sim_errors, save_fig_path = 'Results/Case 8/e_pfba_exp_sim_plot.png', title = 'pFBA Ethanol Carbon Source')
    plt.close('all')

    # #LMOMA
    # e_lmoma_res, e_lmoma_exp_sim, e_lmoma_exp_sim_errors = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = None, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_ethanol_case5.sav')
    # case8.plotExpVsSim(e_lmoma_exp_sim_errors, save_fig_path = 'Results/Case 8/e_lmoma_exp_sim_plot.png', title = 'LMOMA Ethanol Carbon Source')
    # plt.close('all')

    #FVA
    e_fva_res, e_fva_exp_sim, _ = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = None, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_ethanol_case5.sav')




    #TESTS
    # pd.concat([pd.DataFrame(reactions), gly_fba_exp_sim_errors], axis = 1, join = 'inner')
    # result = pd.DataFrame(reactions).join(gly_fba_exp_sim_errors, how = 'inner')

    # R03668 - KEGGID Wrong in cecafdb dataset
    # D-Glyceraldehyde-3-phosphate <==> Glycerol (transport) not present in model
    # Acetyl-CoA <==> Acetyl-CoA-mit	(transport) not present in model


    # case5.checkReaction(case5.convertKeggID('R03668'))
    # case8.checkReaction('r_0454')




