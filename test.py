

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
from case_7 import Case7


def case7Pipeline (plot = False, makeFigs = False):
    # Experimental Fluxes
    exp_dataset, reactions = case7.loadExperimentalRes('Results/Case 7/case7_experimental_fluxes.csv')

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

    if makeFigs:
        case7.multiplePlotO2vsEtOHSaveFigs(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes, folder = 'Results/Case 7/EtOH_figs')
    # Get correct O2 flux according to exp EtOH flux
    fluxes_O2 = case7.getO2flux(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes)
    #case7.printDict(fluxes_O2)

    return exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2

def fbaPipeline (plotGenes = False, plotReacts = False, saveGenesPlot = False):
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

    if plotGenes:
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_fba_exp_sim_errors, pdf_filename = 'Results/Case 7/fba_genes_exp_vs_sim_plots.pdf')
    if plotReacts:
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_fba_exp_sim, pdf_filename = 'Results/Case 7/fba_reacts_exp_vs_sim_plots.pdf', folder = 'Results/Case 7/FBA_figs')
    if saveGenesPlot:
        case7.multipleGenesPlotExpVsSimSaveFigs(df_fba_exp_sim_errors, folder = 'Results/Case 7/FBA_figs')

    return res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors

def pfbaPipeline (plotGenes = False, plotReacts = False, saveGenesPlot = False):
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

    if plotGenes:
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_pfba_exp_sim_errors, pdf_filename = 'Results/Case 7/pfba_genes_exp_vs_sim_plots.pdf')
    if plotReacts:
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_pfba_exp_sim, pdf_filename = 'Results/Case 7/pfba_reacts_exp_vs_sim_plots.pdf', folder = 'Results/Case 7/pFBA_figs')
    if saveGenesPlot:
        case7.multipleGenesPlotExpVsSimSaveFigs(df_pfba_exp_sim_errors, folder = 'Results/Case 7/pFBA_figs')

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
    exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2 = case7Pipeline(plot = False, makeFigs = False)

    # FBA
    res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors = fbaPipeline(saveGenesPlot = False)

    # pFBA
    res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors = pfbaPipeline(saveGenesPlot = False, plotReacts = False)

    # Create xlsx with results
    #case7.convertPandasDFToExcel(reactions, df_fba_exp_sim_errors, filename = 'Results/Case 7/fba_results_case7_test.xlsx', imgFolder = 'Results/Case 7/FBA_figs')
    case7.convertPandasDFToExcel(reactions, df_pfba_exp_sim_errors, title = 'pFBA Results Case 7', filename = 'Results/Case 7/pfba_results_case7_test.xlsx', imgFolder = 'Results/Case 7/pFBA_figs')













