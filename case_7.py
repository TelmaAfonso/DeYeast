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

from yeastpack_test import loadObjectFromFile, saveObjectToFile, translate_medium_modified, createResultsDataset, createResultsDatasetFVA, convertStdToSyst


def dictsForCase7 ():
    d = {'ADH3': -14.4, 'ALD5': -16.2, 'ALD6': -4.7, 'COX5A': -13.3, 'CTP1': -13.9, 'DAL7': -14.3,
         'FUM1': -9.1, 'GND2': -15.5, 'GCV2': -16, 'GLY1': -13.1, 'GPD1': -16.4, 'ICL1': -15.7,
         'IDP1': -15, 'IDP2': -13.4, 'LSC1': -17.8, 'MAE1': -12.5, 'MDH1': -10, 'MHD2': -13,
         'MDH3': -14.4, 'MSL1': -17.5, 'OAC1': -11.7, 'PCK1': -13.4, 'PDA1': -8, 'PGM1': -14.7,
         'PGM2': -16.1, 'RPE1': -4.6, 'SDH1': -13.4, 'SER33': -14.5, 'SFC1': -12.9, 'SOL1': -14.4,
         'SOL2': -15.6, 'SOL3': -12.2, 'SOL4': -15.9, 'TAL1': -12.2, 'YGR043C': -16, 'ZWF1': -6.5}

    l = convertStdToSyst(d.keys())              #Dict with gene : yeastpack gene
    l_inv = {v: k for k, v in l.items()}        #Dict with yeastpack gene : gene
    g_lb = {v: d[k] for k, v in l.items()}      #Dict with yeastpack gene : lb

    return d, l, l_inv, g_lb


def case7fba (model, g_lb, l_inv):
    # g_lb - dict with yeastpack gene : lb
    # l_inv - dict with yeastpack gene : gene
    res = {}
    for g, lb in g_lb.items():
        print('======= ' + l_inv[g] + ' (' + g + ')' ' =======')
        with model as m:
            m.set_carbon_source('r_1714', lb = lb)
            m.set_environmental_conditions(gene_knockout = g)
            key = ''.join(l_inv[g] + ' (%s)' % g)
            res[key] = fba(m)
            print('Objective value:', res[key].objective_value, '\n')
    print('======= ' + 'Wild Type' + ' =======')
    with model as m:
        m.set_carbon_source('r_1714', lb = -16.7) # GLC WT = -16.7
        res['WildType'] = fba(m)
        print('Objective value:', res['WildType'].objective_value, '\n')
    return res


def case7pfba (model, g_lb, l_inv):
    # g_lb - dict with yeastpack gene : lb
    # l_inv - dict with yeastpack gene : gene
    res_pfba = {}
    for g, lb in g_lb.items():
        print('======= ' + l_inv[g] + ' (' + g + ')' ' =======')
        with model as m:
            m.set_carbon_source('r_1714', lb = lb)
            m.set_environmental_conditions(gene_knockout = g)
            key = ''.join(l_inv[g] + ' (%s)' % g)
            res_pfba[key] = pfba(m)
            print('Objective value:', res_pfba[key].objective_value, '\n')
    print('======= ' + 'Wild Type' + ' =======')
    with model as m:
        m.set_carbon_source('r_1714', lb = -16.7)
        res_pfba['WildType'] = pfba(m)
        print('Objective value:', res_pfba['WildType'].objective_value, '\n')
    return res_pfba


def case7fva (model, g_lb, l_inv):
    # g_lb - dict with yeastpack gene : lb
    # l_inv - dict with yeastpack gene : gene
    res_fva = {}
    for g, lb in g_lb.items():
        print('======= ' + l_inv[g] + ' (' + g + ')' ' =======')
        with model as m:
            try:
                m.set_carbon_source('r_1714', lb = lb)
                m.set_environmental_conditions(gene_knockout = g)
                key = ''.join(l_inv[g] + ' (%s)' % g)
                res_fva[key] = fva(m, reaction_list = m.reactions, fix_biomass = True)
                print('Done!', '\n')
            except:
                print('FVA solution status infeasible for case ' + key + '\n')
    print('======= ' + 'Wild Type' + ' =======')
    with model as m:
        m.set_carbon_source('r_1714', lb = -16.7)
        res_fva['WildType'] = fva(m, reaction_list = m.reactions, fix_biomass = True)
        print('Done!', '\n')
    return res_fva


def readExpFluxes (filename, sep = ';'):
    exp_res = pd.read_csv(filename, sep = sep)
    col_names = list(exp_res.columns)

    genes = [gene for gene in col_names if len(gene) <= 5]

    for i in range(0, len(col_names)):
        if 'Real' in col_names[i]:
            col_names[i] = col_names[i-1] + '_real_flux'

    exp_res.columns = col_names

    return exp_res, genes


def relativeError (actual, measured):
    return (absoluteError(actual, measured)/np.absolute(actual))*100


def absoluteError (actual, measured):
    return np.absolute(actual - measured)


def scatterPlot (real_flux, sim_flux, xlab = 'real_flux', ylab = 'sim_flux', title = 'Gene', abs = False):
    if abs:
        real_flux = real_flux.abs()
        sim_flux = sim_flux.abs()
    slope, intercept, r_value, p_value, std_err = linregress(real_flux, sim_flux)
    line = slope * real_flux + intercept
    plt.plot(real_flux, sim_flux, 'o', real_flux, line)
    #plt.scatter(real_flux.abs(), sim_flux.abs())
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.legend(['R_squared: ' + str(r_value**2)])
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def loadExperimentalRes (filename_exp):
    exp_fluxes, genes = readExpFluxes(filename_exp)
    exp_fluxes.set_index('Yeast7_ID', inplace = True)
    sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if '_real_flux' not in col], 1)
    sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.')).astype('float')

    return sub_exp_fluxes


def createDatasetExpVsSimul (dataset_exp, dataset_sim):
    # For FBA and pFBA
    dsim = dataset_sim.copy()
    dsim.columns = [n.split()[0] + '_sim_flux' for n in list(dsim.columns)]
    df = pd.concat([dataset_exp, dsim], axis = 1, join = 'inner')
    df = df[list(sum(zip(dataset_exp.columns, dsim.columns), ()))]
    #df = pd.concat([exp_fluxes['Reaction'], df], axis = 1)

    return df


def createDatasetExpVsSimulFVA (dataset_exp, dataset_sim):
    dsim = res_fva_df.copy()
    dsim.columns = [(i.split()[0] + '_minimum' if count % 2 == 1 else (i.split()[0] + '_maximum')) for count, i in enumerate(list(dsim.columns))]
    df = pd.concat([exp_dataset, dsim], axis = 1, join = 'inner')
    df = df.reindex_axis(sorted(df.columns), axis = 1)

    return df


def createDatasetWithAbsRelError (dataset_exp_vs_sim):
    df = dataset_exp_vs_sim.copy()
    for ind in range(len(list(dataset_exp_vs_sim.columns)), 0,  -2):
        ae = absoluteError(df.ix[:, ind - 2], df.ix[:, ind - 1])
        colname_ae = df.ix[:, ind - 2].name.split('_')[0] + '_abs_error'
        df.insert(loc = ind, column = colname_ae, value = ae)

        re = relativeError(df.ix[:, ind - 2], df.ix[:, ind - 1])
        colname_re = df.ix[:, ind - 2].name.split('_')[0] + '_rel_error'
        df.insert(loc = ind + 1, column = colname_re, value = re)

    return df


def getBiomassObj (res_dict):
    biom = {}

    for gene, res in res_dict.items():
        biom[gene] = res.objective_value

    return  biom


def getEthanolFux (res_df, reaction_id = 'r_2115'):
    df = res_df.copy()
    sub_df = df.filter(regex = '_real_flux')
    sub_df.columns = [col.split('_')[0] for col in list(sub_df.columns)]
    return dict(sub_df.ix[reaction_id, :])


def checkReaction (reaction_id):
    print('Reation (' + reaction_id + '): ' + model.reactions.get_by_id(reaction_id).name)
    r = model.reactions.get_by_id(reaction_id).reaction.split()
    r2 = [(model.metabolites.get_by_id(met).name if 's' in met else (met)) for met in r]
    print('\n', ' '.join(r2))


def divide_list(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def testO2EthanolProd (model, g_knockout = None, gluc_lb = -10, range_o2 = list(np.arange(-20, 0, 2))):
    res = {}
    for i in range_o2:
        with model as m:
            m.set_carbon_source('r_1714', lb = gluc_lb)                   # glucose
            m.reactions.get_by_id('r_1992').lower_bound = float(i)        # oxygen exchange
            m.set_environmental_conditions(gene_knockout = g_knockout)
            r = fba(m)
            res[str(i)] = r.fluxes['r_2115']
    for key, val in sorted(res.items()): print(key, '\t', val)
    return res


def plotO2vsEtOH (dict_gene_EtOH, dict_real_EtOH_fulxes, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', gene = 'ADH3'):
    plt.figure(figsize = (10, 10))
    x = sorted([int(x) for x in dict_gene_EtOH.keys()])
    y = [dict_gene_EtOH[str(key)] for key in x]
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    line = [slope * x + intercept for x in x]
    #plt.plot(x, y, 'o')
    plt.plot(x, y, 'o', x, line)
    plt.axhline(y = dict_real_EtOH_fulxes[gene])
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.legend([gene, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f' % dict_real_EtOH_fulxes[gene]])
    plt.show()


def multiplePlotO2vsEtOH (dict_EtOH_fluxes, dict_real_EtOH_fluxes, n = 3, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', pdf_filename = None):
    plt.rcParams["figure.figsize"] = (20,3)
    g_list = [genes for genes in iter(divide_list(sorted(list(dict_EtOH_fluxes.keys())), n))]

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
            #plt.plot(x, y, 'o')
            plt.plot(x, y, 'o', x, line)
            plt.axhline(y = dict_real_EtOH_fluxes[gene], ls = 'dashed')
            # plt.axvline(x = real_O2(y0), ls = 'dashed')
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


def fixEtO2FluxesForPlotting (dict_EtOH_fuxes):
    # Remove O2 fluxes that yield no EtOH (for plot slope reasons)
    res = dict_EtOH_fuxes.copy()
    for key in list(dict_EtOH_fuxes.keys()):
        for k in list(dict_EtOH_fuxes[key].keys()):
            if res[key][k] == 0:
                del res[key][k]
    return res


def getO2flux (dict_EtOH_fluxes, dict_real_EtOH_fluxes):
    res = {}
    for gene in dict_EtOH_fluxes.keys():
        x = sorted([int(x) for x in dict_EtOH_fluxes[gene].keys()])
        y = [dict_EtOH_fluxes[gene][str(key)] for key in x]
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        real_O2 = lambda x0: (y0 - intercept) / slope
        y0 = dict_real_EtOH_fluxes[gene]
        res[gene] = format(real_O2(y0), '.2f')

    return  res






if __name__ == '__main__':

    d, l, l_inv, g_lb = dictsForCase7()
    model = loadObjectFromFile('model_yeast_76.sav')
    model.set_medium(translate_medium_modified('MINIMAL_CASE7', model))

    # FBA
    res_fba = case7fba(model, g_lb, l_inv)
    saveObjectToFile(res_fba, 'Results/Case 7/res_fba_case7.sav')
    res_fba = loadObjectFromFile('Results/Case 7/res_fba_case7.sav')
    res_fba_df, wt_fba_df = createResultsDataset(res_fba)
    res_fba_df.to_csv('Results/Case 7/res_fba_case7.csv', sep = ';')

    # pFBA
    res_pfba = case7pfba(model, g_lb, l_inv)
    saveObjectToFile(res_pfba, 'Results/Case 7/res_pfba_case7.sav')
    res_pfba = loadObjectFromFile('Results/Case 7/res_pfba_case7.sav')
    res_pfba_df, wt_pfba_df  = createResultsDataset(res_pfba)
    res_pfba_df.to_csv('Results/Case 7/res_pfba_case7.csv', sep = ';')

    # FVA
    res_fva = case7fva(model, g_lb, l_inv)
    saveObjectToFile(res_fva, 'Results/Case 7/res_fva_case7.sav')
    res_fva = loadObjectFromFile('Results/Case 7/res_fva_case7.sav')
    res_fva_df, wt_fva_df  = createResultsDatasetFVA(res_fva)
    res_fva_df.to_csv('Results/Case 7/res_fva_case7.csv', sep = ';')


    # Datasets with experimental vs simulated fluxes
    exp_dataset = loadExperimentalRes('Results/Case 7/case7_experimental_fluxes.csv')
    df_fba = createDatasetExpVsSimul(exp_dataset, res_fba_df)
    df_fba.to_csv('Results/Case 7/exp_vs_sim_fba_case7.csv', sep = ';')
    df_pfba = createDatasetExpVsSimul(exp_dataset, res_pfba_df)
    df_fva = createDatasetExpVsSimulFVA(exp_dataset, res_fva_df)


    #Plots
    scatterPlot(df_fba['ADH3_real_flux'], df_fba['ADH3_sim_flux'], title = 'ADH3 (FBA)')
    scatterPlot(df_pfba['ADH3_real_flux'], df_pfba['ADH3_sim_flux'], title = 'ADH3 (pFBA)')
    plt.close()
    scatterPlot(df_fba['ADH3_real_flux'], df_fba['ADH3_sim_flux'], title = 'ADH3', abs = True)
    plt.close()

    # Datasets with absolute and realtive errors
    df_fba_errors = createDatasetWithAbsRelError(df_fba)
    df_pfba_errors = createDatasetWithAbsRelError(df_pfba)
    list(df_pfba_errors.columns)

    # Biomass Values
    fba_biomass = getBiomassObj(res_fba)
    pfba_biomass = getBiomassObj(res_pfba)

    # Check reactions
    checkReaction('r_4041') #Biomass
    checkReaction('r_2115') #EtOH

    # Get ethanol experimental values
    et_fba = getEthanolFux(df_fba)
    et_pfba = getEthanolFux(df_pfba)

    # Testing EtOH fluxes with O2 consumption
    genes = sorted(list(l_inv.values()))
    range_o2 = list(np.arange(-20, 0, 2))

    res_EtOH_fuxes = {}
    for gene in genes:
        print('Gene ' + gene + ':')
        res_EtOH_fuxes[gene] = testO2EthanolProd(model, g_knockout = l[gene], gluc_lb = d[gene], range_o2 = range_o2)
        print('Done!')

    saveObjectToFile(res_EtOH_fuxes, 'Results/Case 7/res_EtOH_fuxes.csv')
    res_EtOH_fuxes = loadObjectFromFile('Results/Case 7/res_EtOH_fuxes.csv')

    # Get real EtOH fluxes
    real_EtOH_fluxes = getEthanolFux(df_fba, 'r_2115')
    real_EtOH_fluxes['ADH3']

    #Plot results (horizontal line for real flux)
    res_EtOH_fuxes_fix = fixEtO2FluxesForPlotting(res_EtOH_fuxes) # Remove O2 flux with no EtOH yield
    plotO2vsEtOH(res_EtOH_fuxes_fix['ADH3'], real_EtOH_fluxes, gene = 'ADH3') # For ADH3 gene
    plots_EtOH = multiplePlotO2vsEtOH(res_EtOH_fuxes_fix, real_EtOH_fluxes, pdf_filename = 'EtOH_fluxes_plot.pdf')

    #Get correct O2 flux according to exp EtOH flux
    fluxes_O2 = getO2flux(res_EtOH_fuxes_fix, real_EtOH_fluxes)





    # Next:
    closest_val = lambda num, list: min(list, key = lambda x: abs(x-num))
    closest_val(5, [4,1,88,44,3])




