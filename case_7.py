#!/usr/bin/env python

'''
File case 7 simulations

Author: Telma Afonso
'''

from yeastpack.simulation import fba, fva, pfba, lmoma
from yeastpack.data import Media
from types import *
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.backends.backend_pdf import PdfPages

from yeastpack_test import YeastpackSim


class Case7 (YeastpackSim):

    def __int__(self, cobra_model = None):
        super(Case7, self).__int__(self, cobra_model)
        self.fba_res = None
        self.pfba_res = None
        self.fva_res = None
        self.dictsForCase7() #Not working for some reason

    def dictsForCase7 (self):
        # Glucose lbs
        d = {'ADH3': -14.4, 'ALD5': -16.2, 'ALD6': -4.7, 'COX5A': -13.3, 'CTP1': -13.9, 'DAL7': -14.3,
             'FUM1': -9.1, 'GND2': -15.5, 'GCV2': -16, 'GLY1': -13.1, 'GPD1': -16.4, 'ICL1': -15.7,
             'IDP1': -15, 'IDP2': -13.4, 'LSC1': -17.8, 'MAE1': -12.5, 'MDH1': -10, 'MHD2': -13,
             'MDH3': -14.4, 'MSL1': -17.5, 'OAC1': -11.7, 'PCK1': -13.4, 'PDA1': -8, 'PGM1': -14.7,
             'PGM2': -16.1, 'RPE1': -4.6, 'SDH1': -13.4, 'SER33': -14.5, 'SFC1': -12.9, 'SOL1': -14.4,
             'SOL2': -15.6, 'SOL3': -12.2, 'SOL4': -15.9, 'TAL1': -12.2, 'YGR043C': -16, 'ZWF1': -6.5}

        exp_biomass = {'ADH3':0.3716, 'ALD5': 0.4126, 'ALD6': 0.1387, 'COX5A': 0.256, 'CTP1': 0.37, 'DAL7': 0.35,
                       'FUM1': 0.2097, 'GND2': 0.338, 'GCV2': 0.3745, 'GLY1': 0.318, 'GPD1': 0.34, 'ICL1': 0.386,
                       'IDP1': 0.372, 'IDP2': 0.356, 'LSC1': 0.423, 'MAE1': 0.411, 'MDH1': 0.29, 'MHD2': 0.36,
                       'MDH3': 0.405, 'MSL1': 0.4, 'OAC1': 0.29, 'PCK1': 0.4, 'PDA1': 0.17, 'PGM1': 0.37,
                       'PGM2': 0.41, 'RPE1': 0.1316, 'SDH1': 0.293, 'SER33': 0.33, 'SFC1': 0.339, 'SOL1': 0.37,
                       'SOL2': 0.402, 'SOL3': 0.3557, 'SOL4': 0.3843, 'TAL1': 0.361, 'YGR043C': 0.3745, 'ZWF1': 0.16}

        l = self.convertStdToSyst(d.keys())         #Dict with gene : yeastpack gene
        l_inv = {v: k for k, v in l.items()}        #Dict with yeastpack gene : gene
        g_lb = {v: d[k] for k, v in l.items()}      #Dict with yeastpack gene : lb
        self.d, self.l, self.l_inv, self.g_lb, self.exp_biomass = d, l, l_inv, g_lb, exp_biomass

    def case7fba (self, o2_lb = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res = {}
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                m.set_carbon_source('r_1714', lb = lb)
                m.set_environmental_conditions(gene_knockout = g)
                if o2_lb is not None:
                    m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                key = ''.join(self.l_inv[g] + ' (%s)' % g)
                res[key] = fba(m)
                print('Objective value:', res[key].objective_value, '\n')
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7) # GLC WT = -16.7
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res['WildType'] = fba(m)
            print('Objective value:', res['WildType'].objective_value, '\n')
        return res

    def case7pfba (self, o2_lb = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res_pfba = {}
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                m.set_carbon_source('r_1714', lb = lb)
                m.set_environmental_conditions(gene_knockout = g)
                if o2_lb is not None:
                    m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                key = ''.join(self.l_inv[g] + ' (%s)' % g)
                res_pfba[key] = pfba(m)
                print('Objective value:', res_pfba[key].objective_value, '\n')
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res_pfba['WildType'] = pfba(m)
            print('Objective value:', res_pfba['WildType'].objective_value, '\n')
        return res_pfba

    def case7fva (self, o2_lb = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res_fva = {}
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                try:
                    m.set_carbon_source('r_1714', lb = lb)
                    m.set_environmental_conditions(gene_knockout = g)
                    if o2_lb is not None:
                        m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                    key = ''.join(self.l_inv[g] + ' (%s)' % g)
                    res_fva[key] = fva(m, reaction_list = m.reactions, fix_biomass = True)
                    print('Done!', '\n')
                except:
                    print('FVA solution status infeasible for case ' + key + '\n')
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res_fva['WildType'] = fva(m, reaction_list = m.reactions, fix_biomass = True)
            print('Done!', '\n')
        return res_fva

    def readExpFluxes (self, filename, sep = ';'):
        exp_res = pd.read_csv(filename, sep = sep)
        col_names = list(exp_res.columns)

        genes = [gene for gene in col_names if len(gene) <= 5]

        for i in range(0, len(col_names)):
            if 'Real' in col_names[i]:
                col_names[i] = col_names[i-1] + '_real_flux'

        exp_res.columns = col_names

        return exp_res, genes

    def loadExperimentalRes (self, filename_exp):
        exp_fluxes, genes = self.readExpFluxes(filename_exp)
        exp_fluxes.set_index('Yeast7_ID', inplace = True)
        sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if '_real_flux' not in col], 1)
        sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.')).astype('float')

        return sub_exp_fluxes

    def createDatasetExpVsSimul (self, dataset_exp, dataset_sim):
        # For FBA and pFBA
        dsim = dataset_sim.copy()
        dsim.columns = [n.split()[0] + '_sim_flux' for n in list(dsim.columns)]
        df = pd.concat([dataset_exp, dsim], axis = 1, join = 'inner')
        df = df[list(sum(zip(dataset_exp.columns, dsim.columns), ()))]

        return df

    def createDatasetExpVsSimulFVA (self, dataset_exp, dataset_sim):
        dsim = dataset_sim.copy()
        dsim.columns = [(i.split()[0] + '_minimum' if count % 2 == 1 else (i.split()[0] + '_maximum')) for count, i in enumerate(list(dsim.columns))]
        df = pd.concat([dataset_exp, dsim], axis = 1, join = 'inner')
        df = df.reindex_axis(sorted(df.columns), axis = 1)

        return df

    def createDatasetWithAbsRelError (self, dataset_exp_vs_sim):
        df = dataset_exp_vs_sim.copy()
        for ind in range(len(list(dataset_exp_vs_sim.columns)), 0,  -2):
            ae = self.absoluteError(df.ix[:, ind - 2], df.ix[:, ind - 1])
            colname_ae = df.ix[:, ind - 2].name.split('_')[0] + '_abs_error'
            df.insert(loc = ind, column = colname_ae, value = ae)

            re = self.relativeError(df.ix[:, ind - 2], df.ix[:, ind - 1])
            colname_re = df.ix[:, ind - 2].name.split('_')[0] + '_rel_error'
            df.insert(loc = ind + 1, column = colname_re, value = re)

        return df

    def getBiomassObj (self, res_dict):
        biom = {}

        for gene, res in res_dict.items():
            biom[gene] = res.objective_value

        return  biom

    def getEthanolFux (self, dataset_exp, reaction_id = 'r_2115'):
        df = dataset_exp.copy()
        # sub_df = df.filter(regex = '_real_flux')
        df.columns = [col.split('_')[0] for col in list(df.columns)]
        # df.columns = list(sorted(l_inv.keys()))
        return dict(df.ix[reaction_id, :])

    def testO2EthanolProd (self, do_fba = True, g_knockout = None, gluc_lb = -10, range_o2 = list(np.arange(-20, 0, 2))):
        res = {}
        for i in range_o2:
            with self.model as m:
                m.set_carbon_source('r_1714', lb = gluc_lb)                   # glucose
                m.reactions.get_by_id('r_1992').lower_bound = float(i)        # oxygen exchange
                m.set_environmental_conditions(gene_knockout = g_knockout)
                if do_fba: r = fba(m)
                else: r = pfba(m)
                res[str(i)] = r.fluxes['r_2115']
        for key, val in sorted(res.items()): print(key, '\t', val)
        return res

    def testAllO2EthanolProd (self, range_o2 = list(np.arange(-20, 0, 2))):
        res_EtOH_fuxes = {}
        genes = sorted(list(self.l_inv.values()))
        for gene in genes:
            print('Gene ' + gene + ':')
            res_EtOH_fuxes[gene] = self.testO2EthanolProd(g_knockout = self.l[gene], gluc_lb = self.d[gene], range_o2 = range_o2)
            print('Done!')
        print('Wild Type:')
        res_EtOH_fuxes['WildType'] = self.testO2EthanolProd(g_knockout = None, gluc_lb = -16.7, range_o2 = range_o2)
        print('Done!')
        return res_EtOH_fuxes

    def fixEtO2FluxesForPlotting (self, dict_EtOH_fuxes):
        # Remove O2 fluxes that yield no EtOH (for plot slope reasons)
        res = dict_EtOH_fuxes.copy()
        for key in list(dict_EtOH_fuxes.keys()):
            for k in list(dict_EtOH_fuxes[key].keys()):
                if res[key][k] == 0:
                    del res[key][k]
        return res

    def getO2flux (self, dict_EtOH_fluxes, dict_real_EtOH_fluxes):
        res = {}
        for gene in dict_EtOH_fluxes.keys():
            x = sorted([int(x) for x in dict_EtOH_fluxes[gene].keys()])
            y = [dict_EtOH_fluxes[gene][str(key)] for key in x]
            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            real_O2 = lambda x0: (y0 - intercept) / slope
            y0 = dict_real_EtOH_fluxes[gene]
            res[gene] = format(real_O2(y0), '.2f')

        return res

    def correctReversedReactions (self, dataset, reactions = None):
        if reactions is None:
            reactions = ['r_0962', 'r_0300', 'r_1022', 'r_1054', 'r_0452', 'r_0892', 'r_0893', 'r_1049', 'r_1048']

        dat = dataset.copy()
        dat.update(dat.loc[reactions] * -1)

        return dat


    # PLOTS

    def plotO2vsEtOH (self, dict_gene_EtOH, dict_real_EtOH_fulxes, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', gene = 'ADH3'):
        plt.figure(figsize = (10, 10))
        x = sorted([int(x) for x in dict_gene_EtOH.keys()])
        y = [dict_gene_EtOH[str(key)] for key in x]
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        plt.plot(x, y, 'o', x, line)
        plt.axhline(y = dict_real_EtOH_fulxes[gene])
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.legend([gene, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f' % dict_real_EtOH_fulxes[gene]])
        plt.show()

    def multiplePlotO2vsEtOH (self, dict_EtOH_fluxes, dict_real_EtOH_fluxes, n = 3, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', pdf_filename = None):
        plt.rcParams["figure.figsize"] = (20,3)
        g_list = [genes for genes in iter(self.divide_list(sorted(list(dict_EtOH_fluxes.keys())), n))]

        for g_ind in range(1, len(g_list)):
            plt.figure(g_ind, figsize = (10, 10))
            i = 1
            for gene in g_list[g_ind]:
                x = sorted([int(x) for x in dict_EtOH_fluxes[gene].keys()])
                y = [dict_EtOH_fluxes[gene][str(key)] for key in x]
                slope, intercept, r_value, p_value, std_err = linregress(x, y)
                line = [slope * x + intercept for x in x]
                real_O2 = lambda x0: (y0 - intercept) / slope
                y0 = dict_real_EtOH_fluxes[gene]
                plt.subplot(n, 1, i)
                plt.plot(x, y, 'o', x, line)
                plt.axhline(y = dict_real_EtOH_fluxes[gene], ls = 'dashed')
                plt.ylabel(ylab)
                # plt.title(title)
                plt.legend([gene, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f (O2 flux of %.2f)' % (dict_real_EtOH_fluxes[gene], real_O2(y0))])

                i += 1
            plt.xlabel(xlab)

        if pdf_filename is not None:
            pdf = PdfPages(pdf_filename)
            for i in plt.get_fignums():
                pdf.savefig(plt.figure(i))
            pdf.close()

        return [plt.figure(i) for i in plt.get_fignums()]

    def plotGeneExpVsSim (self, absRelErrorDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'ADH3', pdf_filename = None, gene = 'ADH3'):
        plt.rcParams["figure.figsize"] = (10,5)

        x = absRelErrorDataset[gene + '_real_flux']
        y = absRelErrorDataset[gene + '_sim_flux']
        gene_IDs = list(absRelErrorDataset.index)
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        meanRelErr = absRelErrorDataset[gene + '_rel_error'].mean()
        corr = absRelErrorDataset[gene + '_real_flux'].corr(absRelErrorDataset[gene + '_sim_flux'])

        plt.plot(x, y, 'o', x, line)
        for ind, gene_ID in enumerate(gene_IDs):
            plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.plot([], [], ' ') # To show correlation in legend
        plt.plot([], [], ' ') # To show mean relative error in legend
        plt.legend([gene, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

    def multipleGenesPlotExpVsSim (self, absRelErrorDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', pdf_filename = None):
        plt.rcParams["figure.figsize"] = (20,3)

        genes = sorted(set([name.split('_')[0] for name in list(absRelErrorDataset.columns)]))
        g_list = [genes for genes in iter(self.divide_list(genes, n))]

        for g_ind in range(1, len(g_list)):
            plt.figure(g_ind, figsize = (10, 10))
            i = 1
            for gene in g_list[g_ind]:
                x = absRelErrorDataset[gene + '_real_flux']
                y = absRelErrorDataset[gene + '_sim_flux']
                gene_IDs = list(absRelErrorDataset.index)
                slope, intercept, r_value, p_value, std_err = linregress(x, y)
                line = [slope * x + intercept for x in x]
                meanRelErr = absRelErrorDataset[gene + '_rel_error'].mean()
                corr = absRelErrorDataset[gene + '_real_flux'].corr(absRelErrorDataset[gene + '_sim_flux'])

                plt.subplot(n, 1, i)
                plt.plot(x, y, 'o', x, line)
                for ind, gene_ID in enumerate(gene_IDs):
                    plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))

                plt.ylabel(ylab)
                plt.plot([], [], ' ') # To show correlation in legend
                plt.plot([], [], ' ') # To show mean relative error in legend
                plt.legend([gene, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

                i += 1
            plt.xlabel(xlab)

        if pdf_filename is not None:
            pdf = PdfPages(pdf_filename)
            for i in plt.get_fignums():
                pdf.savefig(plt.figure(i))
            pdf.close()

        return [plt.figure(i) for i in plt.get_fignums()]

    def plotReactExpVsSim (self, simVsExpDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'r_0534', reaction = 'r_0534'):
        plt.rcParams["figure.figsize"] = (10,5)
        # Prepare dataset
        df = simVsExpDataset.copy()
        df_sim = df[list(df.columns[1::2])].transpose()
        df_exp = df[list(df.columns[0::2])].transpose()

        x = df_exp[reaction]
        y = df_sim[reaction]
        genes = [gene.split('_')[0] for gene in list(df_exp.index)]
        x.index = y.index = genes
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        meanRelErr = self.relativeError(x, y).mean()
        corr = x.corr(y)

        plt.plot(x, y, 'o', x, line)
        for ind, gene_ID in enumerate(genes):
            plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind], y[ind]))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.plot([], [], ' ') # To show correlation in legend
        plt.plot([], [], ' ') # To show mean relative error in legend
        plt.legend([reaction, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

    def multipleReactsPlotExpVsSim (self, simVsExpDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', pdf_filename = None):
        plt.rcParams["figure.figsize"] = (20,3)
        # Prepare dataset
        df = simVsExpDataset.copy()
        df_sim = df[list(df.columns[1::2])].transpose()
        df_exp = df[list(df.columns[0::2])].transpose()
        r_list = [r for r in iter(self.divide_list(list(df_exp.columns), n))]
        genes = [gene.split('_')[0] for gene in list(df_exp.index)]

        for r_ind in range(1, len(r_list)):
            plt.figure(r_ind, figsize = (10, 10))
            i = 1
            for reaction in r_list[r_ind]:
                x = df_exp[reaction]
                y = df_sim[reaction]
                x.index = y.index = genes
                slope, intercept, r_value, p_value, std_err = linregress(x, y)
                line = [slope * x + intercept for x in x]
                meanRelErr = self.relativeError(x, y).mean()
                corr = x.corr(y)
                plt.subplot(n, 1, i)
                plt.plot(x, y, 'o', x, line)
                for ind, gene_ID in enumerate(genes):
                    plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind], y[ind]))

                plt.ylabel(ylab)
                plt.plot([], [], ' ') # To show correlation in legend
                plt.plot([], [], ' ') # To show mean relative error in legend
                plt.legend([reaction, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

                i += 1
            plt.xlabel(xlab)

        if pdf_filename is not None:
            pdf = PdfPages(pdf_filename)
            for i in plt.get_fignums():
                pdf.savefig(plt.figure(i))
            pdf.close()

        return [plt.figure(i) for i in plt.get_fignums()]

    def plotGeneExpVsSimFVA (self, simVsExpDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'ADH3', pdf_filename = None, gene = 'ADH3'):
        # NOT FINISHED
        plt.rcParams["figure.figsize"] = (10,5)

        x = simVsExpDataset[gene + '_real_flux']
        ymin = simVsExpDataset[gene + '_minimum']
        ymax = simVsExpDataset[gene + '_maximum']

        plt.plot(x, 'o')
        plt.fill_between(range(len(ymin)), ymin, ymax, alpha=.1)

        # react_IDs = list(simVsExpDataset.index)
        #
        # for ind, react_ID in enumerate(react_IDs):
        #     plt.annotate(react_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)



# Pipelines

def case7Pipeline (plot = False):
    # Experimental Fluxes
    exp_dataset = case7.loadExperimentalRes('Results/Case 7/case7_experimental_fluxes.csv')

    # Get real EtOH fluxes
    real_EtOH_fluxes = case7.getEthanolFux(exp_dataset, 'r_2115')
    real_EtOH_fluxes['WildType'] = 23.6318184606019 #From authors file

    # Testing EtOH fluxes with O2 consumption
    # sim_EtOH_O2_fluxes = case7.testAllO2EthanolProd()
    # case7.saveObjectToFile(sim_EtOH_O2_fluxes, 'Results/Case 7/dict_etOH_O2_fluxes.sav')
    sim_EtOH_O2_fluxes = case7.loadObjectFromFile('Results/Case 7/dict_etOH_O2_fluxes.sav')

    # Fix EtOH fluxes with O2 consumption for plotting
    sim_EtOH_O2_fluxes_fixed = case7.fixEtO2FluxesForPlotting(sim_EtOH_O2_fluxes)
    if plot:
        case7.multiplePlotO2vsEtOH (sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes, pdf_filename = 'Results/Case 7/etOH_O2_fluxes_plot.pdf')

    # Get correct O2 flux according to exp EtOH flux
    fluxes_O2 = case7.getO2flux(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes)
    #case7.printDict(fluxes_O2)

    return exp_dataset, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2

def fbaPipeline (plot = None):
    # res_fba = case7.case7fba(fluxes_O2)
    # case7.saveObjectToFile(res_fba, 'Results/Case 7/res_fba_case7.sav')
    res_fba = case7.loadObjectFromFile('Results/Case 7/res_fba_case7.sav')
    res_fba_df, wt_fba_df = case7.createResultsDataset(res_fba)
    res_fba_df = case7.correctReversedReactions(res_fba_df) # Reactions fixed
    wt_fba_df = case7.correctReversedReactions(wt_fba_df)

    # Dataset with experimental vs simulated fluxes
    df_fba_exp_sim = case7.createDatasetExpVsSimul(exp_dataset, res_fba_df)

    # Dataset with absolute and realtive errors
    df_fba_exp_sim_errors = case7.createDatasetWithAbsRelError(df_fba_exp_sim)

    if plot == 'genes':
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_fba_exp_sim_errors, pdf_filename = 'Results/Case 7/fba_genes_exp_vs_sim_plots.pdf')
    elif plot == 'reactions':
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_fba_exp_sim, pdf_filename = 'Results/Case 7/fba_reacts_exp_vs_sim_plots.pdf')

    return res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors

def pfbaPipeline (plot = None):
    # res_pfba = case7.case7pfba(fluxes_O2)
    # case7.saveObjectToFile(res_pfba, 'Results/Case 7/res_pfba_case7.sav')
    res_pfba = case7.loadObjectFromFile('Results/Case 7/res_pfba_case7.sav')
    res_pfba_df, wt_pfba_df = case7.createResultsDataset(res_pfba)
    res_pfba_df = case7.correctReversedReactions(res_pfba_df) # Reactions fixed
    wt_pfba_df = case7.correctReversedReactions(wt_pfba_df)

    # Dataset with experimental vs simulated fluxes
    df_pfba_exp_sim = case7.createDatasetExpVsSimul(exp_dataset, res_pfba_df)

    # Dataset with absolute and realtive errors
    df_pfba_exp_sim_errors = case7.createDatasetWithAbsRelError(df_pfba_exp_sim)

    if plot == 'genes':
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_pfba_exp_sim_errors, pdf_filename = 'Results/Case 7/pfba_genes_exp_vs_sim_plots.pdf')
    elif plot == 'reactions':
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_fba_exp_sim, pdf_filename = 'Results/Case 7/pfba_reacts_exp_vs_sim_plots.pdf')

    return res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors

def fvaPipeline ():
    # res_fva = case7.case7fva(fluxes_O2)
    # case7.saveObjectToFile(res_fva, 'Results/Case 7/res_fva_case7.sav')
    res_fva = case7.loadObjectFromFile('Results/Case 7/res_fva_case7.sav')
    res_fva_df, wt_fva_df = case7.createResultsDatasetFVA(res_fva)
    res_fva_df = case7.correctReversedReactions(res_fva_df) # Reactions fixed
    wt_fva_df = case7.correctReversedReactions(wt_fva_df)

    # Dataset with experimental vs simulated fluxes
    df_fva_exp_sim = case7.createDatasetExpVsSimulFVA(exp_dataset, res_fva_df)

    #NOT FiNISHED - PLOTS NEXT (?)

    # if plot == 'genes':
    #     #Plots w/ errors for genes
    #     # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
    #     case7.multipleGenesPlotExpVsSim(df_pfba_exp_sim_errors, pdf_filename = 'Results/Case 7/pfba_genes_exp_vs_sim_plots.pdf')
    # elif plot == 'reactions':
    #     #Plots w/ errors for reactions
    #     # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
    #     case7.multipleReactsPlotExpVsSim(df_fba_exp_sim, pdf_filename = 'Results/Case 7/pfba_reacts_exp_vs_sim_plots.pdf')

    return res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim





if __name__ == '__main__':

    #Initialization
    case7 = Case7()
    case7.model = case7.loadObjectFromFile('model_yeast_76.sav')
    case7.setMedium('MINIMAL_CASE7')
    case7.dictsForCase7()

    # General datasets
    exp_dataset, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2 = case7Pipeline(plot = False)

    # FBA
    res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors = fbaPipeline(plot = 'genes')

    # pFBA
    res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors = pfbaPipeline()

    # FVA
    res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim = fvaPipeline()


    # Compare biomasses
    biom = case7.getBiomassObj(res_fba)
    case7.printDict(biom)
    case7.printDict(case7.exp_biomass)








    # ========================
    #       OLD STUFF
    # ========================

    #
    # d, l, l_inv, g_lb = dictsForCase7()
    # model = loadObjectFromFile('model_yeast_76.sav')
    # model.set_medium(translate_medium_modified('MINIMAL_CASE7', model))
    #
    # # FBA
    # res_fba = case7fba(model, g_lb, l_inv)
    # saveObjectToFile(res_fba, 'Results/Case 7/res_fba_case7.sav')
    # res_fba = loadObjectFromFile('Results/Case 7/res_fba_case7.sav')
    # res_fba_df, wt_fba_df = createResultsDataset(res_fba)
    # res_fba_df.to_csv('Results/Case 7/res_fba_case7.csv', sep = ';')
    #
    # # pFBA
    # res_pfba = case7pfba(model, g_lb, l_inv)
    # saveObjectToFile(res_pfba, 'Results/Case 7/res_pfba_case7.sav')
    # res_pfba = loadObjectFromFile('Results/Case 7/res_pfba_case7.sav')
    # res_pfba_df, wt_pfba_df  = createResultsDataset(res_pfba)
    # res_pfba_df.to_csv('Results/Case 7/res_pfba_case7.csv', sep = ';')
    #
    # # FVA
    # res_fva = case7fva(model, g_lb, l_inv)
    # saveObjectToFile(res_fva, 'Results/Case 7/res_fva_case7.sav')
    # res_fva = loadObjectFromFile('Results/Case 7/res_fva_case7.sav')
    # res_fva_df, wt_fva_df  = createResultsDatasetFVA(res_fva)
    # res_fva_df.to_csv('Results/Case 7/res_fva_case7.csv', sep = ';')
    #
    #
    # # Datasets with experimental vs simulated fluxes
    # exp_dataset = loadExperimentalRes('Results/Case 7/case7_experimental_fluxes.csv')
    # df_fba = createDatasetExpVsSimul(exp_dataset, res_fba_df)
    # df_fba.to_csv('Results/Case 7/exp_vs_sim_fba_case7.csv', sep = ';')
    # df_pfba = createDatasetExpVsSimul(exp_dataset, res_pfba_df)
    # df_fva = createDatasetExpVsSimulFVA(exp_dataset, res_fva_df)
    #
    #
    # #Plots w/ errors
    # scatterPlot(df_fba['ADH3_real_flux'], df_fba['ADH3_sim_flux'], title = 'ADH3 (FBA)')
    # scatterPlot(df_pfba['ADH3_real_flux'], df_pfba['ADH3_sim_flux'], title = 'ADH3 (pFBA)')
    # plt.close()
    # scatterPlot(df_fba['ADH3_real_flux'], df_fba['ADH3_sim_flux'], title = 'ADH3', abs = True)
    # plt.close()
    #
    # # Datasets with absolute and realtive errors
    # df_fba_errors = createDatasetWithAbsRelError(df_fba)
    # df_pfba_errors = createDatasetWithAbsRelError(df_pfba)
    # list(df_pfba_errors.columns)
    #
    # # Biomass Values
    # fba_biomass = getBiomassObj(res_fba)
    # pfba_biomass = getBiomassObj(res_pfba)
    #
    # # Check reactions
    # checkReaction('r_4041') #Biomass
    # checkReaction('r_2115') #EtOH
    #
    # # Get ethanol experimental values
    # et_fba = getEthanolFux(df_fba)
    # et_pfba = getEthanolFux(df_pfba)
    #
    # # Testing EtOH fluxes with O2 consumption
    # genes = sorted(list(l_inv.values()))
    # range_o2 = list(np.arange(-20, 0, 2))
    #
    # res_EtOH_fuxes = {}
    # for gene in genes:
    #     print('Gene ' + gene + ':')
    #     res_EtOH_fuxes[gene] = testO2EthanolProd(model, g_knockout = l[gene], gluc_lb = d[gene], range_o2 = range_o2)
    #     print('Done!')
    #
    # saveObjectToFile(res_EtOH_fuxes, 'Results/Case 7/res_EtOH_fuxes.csv')
    # res_EtOH_fuxes = loadObjectFromFile('Results/Case 7/res_EtOH_fuxes.csv')
    #
    # # Get real EtOH fluxes
    # real_EtOH_fluxes = getEthanolFux(df_fba, 'r_2115')
    # real_EtOH_fluxes['ADH3']
    #
    # #Plot results (horizontal line for real flux)
    # res_EtOH_fuxes_fix = fixEtO2FluxesForPlotting(res_EtOH_fuxes) # Remove O2 flux with no EtOH yield
    # plotO2vsEtOH(res_EtOH_fuxes_fix['ADH3'], real_EtOH_fluxes, gene = 'ADH3') # For ADH3 gene
    # plots_EtOH = multiplePlotO2vsEtOH(res_EtOH_fuxes_fix, real_EtOH_fluxes, pdf_filename = 'EtOH_fluxes_plot.pdf')
    #
    # #Get correct O2 flux according to exp EtOH flux
    # fluxes_O2 = getO2flux(res_EtOH_fuxes_fix, real_EtOH_fluxes)
    #
    # # FBA with O2 flux fixed
    # res_fba_fix = case7fba(model, g_lb, l_inv, fluxes_O2)
    # saveObjectToFile(res_fba_fix, 'Results/Case 7/res_fba_case7_fixed.sav')
    # res_fba_fix = loadObjectFromFile('Results/Case 7/res_fba_case7_fixed.sav')
    # res_fba_fix_df, wt_fba_fix_df = createResultsDataset(res_fba_fix)
    # res_fba_fix_df = correctReversedReactions(res_fba_fix_df) # Reactions fixed
    #
    # # Dataset with experimental vs simulated fluxes
    # exp_dataset = loadExperimentalRes('Results/Case 7/case7_experimental_fluxes.csv')
    # df_fba_fix = createDatasetExpVsSimul(exp_dataset, res_fba_fix_df)
    #
    # # Datasets with absolute and realtive errors
    # df_fba_fix_errors = createDatasetWithAbsRelError(df_fba_fix)
    # df_fba_fix_errors.to_csv('Results/Case 7/df_fba_fix_errors.csv', sep = ';')
    #
    # #Plots w/ errors for genes
    # plotGeneExpVsSim(df_fba_fix_errors, gene = 'ADH3')
    # multipleGenesPlotExpVsSim(df_fba_fix_errors, pdf_filename = 'fba_exp_vs_sim_plots.pdf')
    #
    # #Plots w/ errors for reactions
    # plotReactExpVsSim(df_fba_fix, reaction = 'r_0534')
    # fba_plots = multipleReactsPlotExpVsSim(df_fba_fix, pdf_filename = 'fba_reacts__exp_vs_sim_plots.pdf')
    #
    #
    # # USING CLASS
    # case7 = Case7(exp_fluxes_file = 'Results/Case 7/case7_experimental_fluxes.csv')
    # case7.model
    #



